#ifndef MESH_H
#define MESH_H

#include "vector"
#include "globals.h"

//Class stores all coordinate information about the mesh
class Mesh {
private:
	int latLength;
	int lonLength;

	int CalculateDt(void);
        int CalculateCellArea(void);

public:
//	Mesh(); //constructor

	std::ostringstream outstring;

	double dLon;
	double dLat;

	Globals * globals;

	Mesh(Globals *);

	int ReturnLatLen(void) const;
	int ReturnLonLen(void) const;

        double ** cellArea;

	std::vector<double> lat;
	std::vector<double> lon;

	std::vector<double> CenterP(int, int); //Lat and lon at center point
	std::vector<double> NorthP(int, int); //Lat and lon at point (i+2,j)
	std::vector<double> SouthP(int, int); //Lat and lon at point (i-2,j)
	std::vector<double> EastP(int, int); //Lat and lon at point (i,j+2)
	std::vector<double> WestP(int, int); //Lat and lon at point (i,j-2)
};

#endif
