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

	std::vector<double> lat;
	std::vector<double> lon;

	int dLat;
	int dLon;
	//int latStaggered = 0;
	
	std::vector<double> opp; //Pointer to opposite cell if at the pole
	Mesh * grid; //pointer to global Mesh object
	std::vector<std::vector<double>> solution; //2D vector array 

	Field(Mesh *,int,int); //Constructor 

	int ReturnFieldLatLen();
	int ReturnFieldLonLen();

	double CenterP(int, int); //solution value at center point
	double NorthP(int, int); //solution value at point (i+2,j)
	double NorthEastP(int, int); //solution value at point (i+2,j+2)
	double SouthP(int, int); //solution value at point (i-2,j)
	double SouthWestP(int, int);
	double EastP(int, int); //solution value at point (i,j+2)
	double WestP(int, int); //solution value at point (i,j-2)

	double SouthWestAvg(int,int);
	double NorthEastAvg(int,int);
};

#endif