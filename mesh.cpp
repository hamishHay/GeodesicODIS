#include "mesh.h"
#include "depth.h"
#include "globals.h"
#include "math.h"
#include "vector"

//Mesh::Mesh():Mesh(2., 2.) {}; //default constructor

Mesh::Mesh(Globals * Globals) //non-default constructor
{

	globals = Globals;

	dLat = globals->dLat.Value();
	dLon = globals->dLon.Value();

	latLength = dLat*2. + 1.;
	dLat = 0.5*180. / dLat;

	std::cout<<"dlat: "<<dLat<<std::endl;


	lonLength = dLon*2.;
	dLon = 0.5*360. / dLon;

	std::cout<<"dlon: "<<dLon<<std::endl;

	//Populate lat vector with co-latitude
	int count = 0;
	int i = 0;
	double currentLat = 90.0;
	while (count < latLength) {
		lat.push_back(currentLat);
		currentLat -= dLat;
		count++;
	}

	//latLength = lat.size(); //EDIT

	for (int j = 0; j < lonLength; j++) lon.push_back(j*dLon);

        CalculateCellArea();

	CalculateDt();
};

int Mesh::CalculateCellArea(void)
{
        cellArea = new double*[ReturnLatLen()];
        for (int i=0; i<ReturnLatLen()-1; i++) {
           cellArea[i] = new double[ReturnLonLen()];
           for (int j=0; j<ReturnLonLen(); j++) {
              cellArea[i][j] = globals->radius.Value()*globals->radius.Value();
              cellArea[i][j] *= sin(radConv*lat[i]) - sin(radConv*lat[i+1]);
              cellArea[i][j] *= radConv*dLon;
              //if (j==0) std::cout<<cellArea[i][j]<<std::endl;
           }
        }
};

int Mesh::ReturnLatLen(void) const
{
	return latLength;
};

int Mesh::ReturnLonLen(void) const
{
	return lonLength;
};

std::vector<double> Mesh::CenterP(int i, int j)
{
	std::vector<double> point = { lat[i], lon[j] };
	return point;
};

std::vector<double> Mesh::NorthP(int i, int j)
{
	std::vector<double> point = { lat[i-2], lon[j] };
	return point;
};

std::vector<double> Mesh::SouthP(int i, int j)
{
	std::vector<double> point;
	if (i == this->ReturnLatLen() - 1) {
		point = { lat[i], lon[this->ReturnLonLen()/2] };
	}
	else {
		point = { lat[i+2], lon[j] };
	}
	return point;
};

std::vector<double> Mesh::EastP(int i, int j)
{
	std::vector<double> point;
	if (j == this->ReturnLonLen() - 1) {
		point = { lat[i], lon[0] };
	}
	else {
		point = { lat[i], lon[j+2] };
	}
	return point;
};

std::vector<double> Mesh::WestP(int i, int j)
{
	std::vector<double> point = { lat[i], lon[j-2] };
	return point;
};

int Mesh::CalculateDt(void){
	double waveSpeed;
	double dt;
	double minDistance;

	outstring << std::endl << "Calculating maximum time step: " << std::endl;

	waveSpeed = sqrt(globals->g.Value()*globals->h.Value());

	outstring << std::endl << "\t\t Average gravitational wave speed = sqrt(gh): " << waveSpeed << " m/s."<<std::endl;

	minDistance = globals->radius.Value() * sin(dLat*radConv) * dLon*radConv;

	dt = minDistance/waveSpeed;

	outstring << "\t\t Time step calculated: " << dt <<" s."<< std::endl;

	if (dt<globals->timeStep.Value()) {
		outstring << std::endl << "\t\t WARNING: Calculated timestep is less than user selected dt = " << globals->timeStep.Value() <<" s." <<std::endl;
		outstring << "\t\t Setting dt to calculated value of " << dt <<" s." <<std::endl;

		globals->timeStep.SetValue(dt);

		globals->Output.Write(OUT_MESSAGE,&outstring);

		outstring << std::endl << "WARNING: Using maximum calculated timestep as user defined dt is too large." << std::endl;
		globals->Output.Write(ERR_MESSAGE,&outstring);
	}
	else {
		outstring << std::endl << "\t\t Calculated timestep is greater than user selected dt = " << globals->timeStep.Value() <<" s." <<std::endl;
		outstring << "\t\t ODIS will keep user defined timestep. You could save time if you chose a larger timestep."<<std::endl;

		globals->Output.Write(OUT_MESSAGE,&outstring);
	}





	return 1;
}
