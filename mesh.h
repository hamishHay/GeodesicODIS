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
};

#endif
