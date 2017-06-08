#ifndef MESH_H
#define MESH_H

#include "globals.h"
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include "mathRoutines.h"

//Class stores all coordinate information about the mesh
class Mesh {
private:
  int ReadMeshFile(void);
  // int InitializeArrays(void);
  int CalcMappingCoords(void);
  int CalcVelocityTransformFactors(void);

public:
//	Mesh(); //constructor
  Mesh(Globals *);
  Mesh(Globals *, int);

	std::ostringstream outstring;

  // Each row in node_pos_sph represents the lat,lon coords of node. Node ID 5
  // corresponds to node_pos_sph(5,x).
  Array2D<double> node_pos_sph;

  // Same as node_pos_sph but in mapped coords
  Array3D<double> node_pos_map;

  // Each row in node_friends contains the ID (index)
  // to all surrounding friends. E.g, the ID's of the neighbouring nodes to
  // node 4 would be accessed via node_friends[4][0], node_friends[4][1], etc.
  Array2D<int> node_friends;


  // Each row contains the coordinates of the centroids to the 5 or 6 surrounding
  // cells for a specific node (e.g., node ID 8 corresponds to row 8).
  // Column 0 holds the coordinates for the centroid formed from the parent
  // node, friend 0, and friend 1. Column 1 holds the coordinates for the centroid
  // formed from the parent node, friend 1, and friend 2, etc.
  Array3D<double> centroid_pos_sph;
  Array3D<double> centroid_pos_map;

  Array2D<double> node_m;
  Array2D<double> centroid_m;

  Array3D<double> node_vel_trans;

	Globals * globals;

};

#endif
