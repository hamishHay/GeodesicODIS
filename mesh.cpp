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

    // TODO - Make sure dLat is an even number

    lonLength = 4*dLat;
	dLon = 360./lonLength;

	std::cout<<"dlon: "<<dLon<<std::endl;
    std::cout<<"lon nodes: "<<lonLength<<std::endl;

	latLength = 2*dLat + 1;
	dLat = 180. / (latLength - 1);

	std::cout<<"dlat: "<<dLat<<std::endl;
    std::cout<<"lat nodes: "<<latLength<<std::endl;

	//Populate lat vector with co-latitude
	int count = 0;
	int i = 0;
	double currentLat = 90.0;
	while (count < latLength) {
		lat.push_back(currentLat);

		currentLat -= dLat;
		count++;
	}

	for (int j = 0; j < lonLength; j++) {
        lon.push_back(j*dLon);
    }

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
