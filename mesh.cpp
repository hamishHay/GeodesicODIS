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

	//Check dlat and dlon give an integer divison into 360
	//if (fmod(180., dLat) == 0)
	//{
	//	latLength = (180 / dLat)+1;
	//}
	//else dLat = NULL; //Raise exception

	latLength = dLat*2. + 1.;
	dLat = 0.5*180. / dLat;

	std::cout<<"dlat: "<<dLat<<std::endl;

	//dLat = 180. / dLat;

	// if (fmod(360., dLon) == 0)
	// {
	// 	lonLength = (int)(360. / dLon);
	// }
	// else dLon = NULL; //Raise exception

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

	//CalculateDt();
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

void Mesh::CalculateDt(void){
	outstring << std::endl << "Calculating time step: " << std::endl;

	double waveSpeed = sqrt(globals->g.Value()*globals->h.Value());

	outstring << std::endl << "\t Average gravitational wave speed = sqrt(gh): " << waveSpeed <<std::endl;

	double dt = 0.5*globals->radius.Value()/(2*waveSpeed)*cos((90 - dLat)*radConv)*dLon*radConv;

	outstring << "\t Time step calculated: " << dt << std::endl;

	globals->timeStep.SetValue(dt);

	globals->Output.Write(OUT_MESSAGE,&outstring);

}
