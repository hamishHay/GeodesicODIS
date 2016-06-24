#ifndef FIELD_H
#define FIELD_H

#include "vector"
#include "mesh.h"

//Class to populate fields (i.e. velocity etc) at HALF grid spacing
class Field {
private:


public:
	int fieldLatLen;
	int fieldLonLen;

	double * lat;
	double * lon;
	double * cosLat;
	double * sinLat;
	double * cos2Lat;
	double * sin2Lat;
	double * cosLon;
	double * sinLon;
	double * cos2Lon;
	double * sin2Lon;

	double dLat;
	double dLon;
	//int latStaggered = 0;


	std::vector<double> opp; //Pointer to opposite cell if at the pole
	Mesh * grid; //pointer to global Mesh object

	double ** solution;

	Field(Mesh *,int,int); //Constructor

	int ReturnFieldLatLen();
	int ReturnFieldLonLen();

	double SouthWestAvg(int, int);
	double NorthEastAvg(int, int);
};

#endif
