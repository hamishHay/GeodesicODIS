#ifndef MESH_H
#define MESH_H

#include "globals.h"
#include "array2d.h"

//Class stores all coordinate information about the mesh
class Mesh {
private:
  int ReadMeshFile(void);
  // int InitializeArrays(void);
  int CalcMappingCoords(void);

public:
//	Mesh(); //constructor
  Mesh(Globals *);
  Mesh(Globals *, int);

	std::ostringstream outstring;

  Array2D<double> node_pos_sph;
  Array2D<double> node_pos_map;
  Array2D<int> node_friends;

	Globals * globals;

};

#endif
