#ifndef MESH_H
#define MESH_H

#include "globals.h"

//Class stores all coordinate information about the mesh
class Mesh {
private:
  int ReadMeshFile(void);
  int InitializeArrays(void);
  int CalcMappingCoords(void);

public:
//	Mesh(); //constructor

	std::ostringstream outstring;

	// double dLon;
	// double dLat;

  double ** node_pos_sph;    // Position array for each node in spherical coords
  double ** node_pos_map;    // Position array for each node in mapping coords

  int ** node_friends;

	Globals * globals;

	Mesh(Globals *);

};

#endif
