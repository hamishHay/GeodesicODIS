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
  int CalcControlVolumeEdgeLengths(void);
  int CalcControlVolumeEdgeCentres(void);
  int CalcControlVolumeEdgeNormals(void);

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

  // Array to hold mapping factor at at each node and its friends
  Array2D<double> node_m;

  // Array to hold mapping factor at at each centroid
  Array2D<double> centroid_m;

  // Array to hold the cos(alpha) and sin(alpha) velocity conversion coefficients
  // from Lee and Macdonald (2009)
  Array3D<double> node_vel_trans;

  // Array to store the edge length of each side to the control volume surrounding
  // the central node
  Array2D<double> control_vol_edge_len;

  // Array to store the midpoint of each side to the control volume surrounding
  // the central node
  Array3D<double> control_vol_edge_centre_pos_map;

  // Array to store the midpoint mapping factor of each side to the control
  // volume surrounding the central node
  Array2D<double> control_vol_edge_centre_m;

  // Array to store the midpoint normal vectors to each side of the control
  // volume surroinding the central node
  Array3D<double> control_vol_edge_normal_map;


	Globals * globals;

};

#endif
