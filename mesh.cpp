#include "mesh.h"
#include "globals.h"
#include "math.h"
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include "mathRoutines.h"

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cstring>

//Mesh::Mesh():Mesh(2., 2.) {}; //default constructor

Mesh::Mesh(Globals * Globals, int N)
  :node_pos_sph(N,2),
   node_pos_map(N,7,2),        // len 7 to include central node coords (even though it's zero)
   node_friends(N,6),
   centroid_pos_sph(N,6,2),
   centroid_pos_map(N,6,2),
   node_m(N,7),                // len 7 to include central (parent) node
   centroid_m(N,6),            // len 6 as there is are only 5 or 6 centroids
   node_vel_trans(N,7,2),      // In the last dimension, element 1 is cos_a, element 2 is sin_a
   control_vol_edge_len(N,6),
   control_vol_edge_centre_pos_map(N,6,2),
   control_vol_edge_centre_m(N,6),
   control_vol_edge_normal_map(N,6,2),
   control_volume_surf_area_map(N)
{

	globals = Globals;                   // define reference to all constants


  // Read in grid file
  ReadMeshFile();

  // Calculate mapping coordinates
  CalcMappingCoords();

  // Calculate velocity transform factors
  CalcVelocityTransformFactors();

  // Find control volume edge lengths in mapping coords
  CalcControlVolumeEdgeLengths();

  // Find control volume edge midpoints in mapping coords
  CalcControlVolumeEdgeCentres();

  // Find control volume edge outward normal unit vectors
  CalcControlVolumeEdgeNormals();

  // Find control volume area for each node
  CalcControlVolumeArea();

};

// Function to calculate the cosine(alpha) and sine(alpha) velocity tranform
// factors from Lee and Macdonald (2009). These factors are constant in time,
// and so are pre-calculated here to be easily accessed during the spatial
// integration.
int Mesh::CalcVelocityTransformFactors(void)
{
  int i, j, node_num, f;
  double * cos_a, * sin_a;
  double lat1, lat2, lon1, lon2;

  node_num = globals->node_num;

  for (i=0; i<node_num; i++)
  {
    lat1 = node_pos_sph(i,0);
    lon1 = node_pos_sph(i,1);

    // Set pointers to address of variables we want to change
    cos_a = &node_vel_trans(i,0,0);
    sin_a = &node_vel_trans(i,0,1);

    // Pass pointers by reference
    velTransform(*cos_a, *sin_a, lat1, lat1, lon1, lon1);

    // Loop through all friends (len is 7 as first element is parent node)
    // Assign node mapping factors and coordinates
    for (j = 1; j<7; j++)
    {
      cos_a = &node_vel_trans(i,j,0);
      sin_a = &node_vel_trans(i,j,1);

      f = node_friends(i,j-1);            // j-1 because len node_friends[i] = 6;
      switch (f) {
        case -1:                          // if node is pentagon center
          *cos_a = -1.0;
          *sin_a = -1.0;
          break;
        default:                          // if node is a hexagon center
          lat2 = node_pos_sph(f, 0);
          lon2 = node_pos_sph(f, 1);

          // assign transform factors values to arrays
          velTransform(*cos_a, *sin_a, lat1, lat2, lon1, lon2);
          break;
      }
    }
  }
}

// Function to convert all spherical coorinate quantities to mapping coordinates
// x and y from Lee and Macdonald (2009). Mapping factors are also calculated
// and stored. These quantities are stored for node positions, centroid positions,
// node mapping factors, and centroid mapping factors.
int Mesh::CalcMappingCoords(void)
{
  int i, j, node_num, f;
  double * x, * y, * m, r;
  double lat1, lat2, lon1, lon2;

  node_num = globals->node_num;
  r = globals->radius.Value();

  for (i=0; i<node_num; i++)
  {
    lat1 = node_pos_sph(i,0);
    lon1 = node_pos_sph(i,1);

    // Set pointers to address of variables we want to change
    m = &node_m(i,0);
    x = &node_pos_map(i,0,0);
    y = &node_pos_map(i,0,1);

    // Pass pointers by reference
    mapAtPoint(*m, *x, *y, lat1, lat1, lon1, lon1, r);

    // Loop through all friends (len is 7 as first element is parent node)
    // Assign node mapping factors and coordinates
    for (j = 1; j<7; j++)
    {
      m = &node_m(i,j);
      x = &node_pos_map(i,j,0);
      y = &node_pos_map(i,j,1);

      f = node_friends(i,j-1);            // j-1 because len node_friends[i] = 6;
      switch (f) {
        case -1:                          // if node is pentagon center
          *m = -1.0;
          *x = -1.0;
          *y = -1.0;
          break;
        default:                          // if node is a hexagon center
          lat2 = node_pos_sph(f, 0);
          lon2 = node_pos_sph(f, 1);

          // assign mapped values to arrays
          mapAtPoint(*m, *x, *y, lat1, lat2, lon1, lon2, r);
          break;
      }
    }
  }

  for (i=0; i<node_num; i++)
  {
    lat1 = node_pos_sph(i,0);
    lon1 = node_pos_sph(i,1);

    // Assign centroid mapping factors and coordinates
    for (j = 0; j<6; j++)
    {
      m = &centroid_m(i,j);
      x = &centroid_pos_map(i,j,0);
      y = &centroid_pos_map(i,j,1);

      f = node_friends(i,j);
      switch (f) {
        case -1:                          // if node is pentagon center
          *m = -1.0;
          *x = -1.0;
          *y = -1.0;
          break;
        default:                          // if node is a hexagon center
          lat2 = centroid_pos_sph(i, j, 0);
          lon2 = centroid_pos_sph(i, j, 1);

          // assign mapped values to arrays
          mapAtPoint(*m, *x, *y, lat1, lat2, lon1, lon2, r);
          break;
      }
    }
  }

  return 1;
};

// Function to find the length of control volume edges for each node, in mapping
// coordinates. Values are stored in control_vol_edge_len 2D array.
int Mesh::CalcControlVolumeEdgeLengths(void)
{
  int i, j, node_num, f, friend_num;
  double * x1, * y1, * x2, * y2, * edge_len;

  node_num = globals->node_num;

  for (i=0; i<node_num; i++)
  {

    f = node_friends(i,5);

    friend_num = 6;                                     // Assume hexagon (6 centroids)
    if (f == -1) {
      friend_num = 5;                                   // Check if pentagon (5 centroids)
      control_vol_edge_len(i,5) = -1.0;
    }

    for (j=0; j<friend_num; j++)                        // Loop through all centroids in the control volume
    {
      f = node_friends(i,j);

      edge_len = &control_vol_edge_len(i,j);            // set pointer to edge length array

      x1 = &centroid_pos_map(i, j, 0);                  // set map coords for first centroid
      y1 = &centroid_pos_map(i, j, 1);

      x2 = &centroid_pos_map(i, (j+1)%friend_num, 0);   // set map coords for second centroid
      y2 = &centroid_pos_map(i, (j+1)%friend_num, 1);   // automatically loops around using %

      distanceBetween(*edge_len, *x1, *x2, *y1, *y2);   // calculate distance between the two centroids.
                                                        // Edge_len automatically assigned the length
    }
  }

  return 1;
};

// Function to find the midpoint of control volume edges for each node, in
// mapping coordinates. Values are stored in control_vol_edge_centre_pos_map 2D
// array.
int Mesh::CalcControlVolumeEdgeCentres(void)
{
  int i, j, node_num, f, friend_num;
  double * x1, * y1, * x2, * y2, * xc, * yc, * m;
  double lat1, lat2, lon1, lon2;

  node_num = globals->node_num;

  for (i=0; i<node_num; i++)
  {

    f = node_friends(i,5);

    friend_num = 6;                                     // Assume hexagon (6 centroids)
    if (f == -1) {
      friend_num = 5;                                   // Check if pentagon (5 centroids)
      control_vol_edge_len(i,5) = -1.0;
      control_vol_edge_centre_m(i,5) = -1.0;
    }

    for (j=0; j<friend_num; j++)                        // Loop through all centroids in the control volume
    {
      f = node_friends(i,j);

      xc = &control_vol_edge_centre_pos_map(i,j,0);     // set pointer to edge length array
      yc = &control_vol_edge_centre_pos_map(i,j,1);

      x1 = &centroid_pos_map(i, j, 0);                  // set map coords for first centroid
      y1 = &centroid_pos_map(i, j, 1);

      x2 = &centroid_pos_map(i, (j+1)%friend_num, 0);   // set map coords for second centroid
      y2 = &centroid_pos_map(i, (j+1)%friend_num, 1);   // automatically loops around using %

      midpointBetween(*xc, *yc, *x1, *x2, *y1, *y2);    // calculate center coords between the two centroids.
                                                        // xc and yc automatically assigned the coords

      m = &control_vol_edge_centre_m(i,j);              // assign pointer to midpoint map factor

      lat1 = centroid_pos_sph(i, j, 0);                 // get first centroid coords
      lon1 = centroid_pos_sph(i, j, 1);

      lat2 = centroid_pos_sph(i, (j+1)%friend_num, 0);  // get second centroid coords
      lon2 = centroid_pos_sph(i, (j+1)%friend_num, 1);

      mapFactorAtPoint(*m, lat1, lat2, lon1, lon2);     // calculate map factor and store
    }
  }

  return 1;
};

// Function to find the midpoint of control volume edges for each node, in
// mapping coordinates. Values are stored in control_vol_edge_centre_pos_map 2D
// array.
int Mesh::CalcControlVolumeEdgeNormals(void)
{
  int i, j, node_num, f, friend_num;
  double * x1, * y1, * x2, * y2, * xn, * yn;

  node_num = globals->node_num;

  for (i=0; i<node_num; i++)
  {

    f = node_friends(i,5);

    friend_num = 6;                                     // Assume hexagon (6 centroids)
    if (f == -1) {
      friend_num = 5;                                   // Check if pentagon (5 centroids)
      control_vol_edge_len(i,5) = -1.0;
      control_vol_edge_normal_map(i,5,0) = -1.0;        // set pointer to edge length array
      control_vol_edge_normal_map(i,5,1) = -1.0;
    }

    for (j=0; j<friend_num; j++)                        // Loop through all centroids in the control volume
    {
      f = node_friends(i,j);

      xn = &control_vol_edge_normal_map(i,j,0);         // set pointer to edge length array
      yn = &control_vol_edge_normal_map(i,j,1);

      x1 = &centroid_pos_map(i, j, 0);                  // set map coords for first centroid
      y1 = &centroid_pos_map(i, j, 1);

      x2 = &centroid_pos_map(i, (j+1)%friend_num, 0);   // set map coords for second centroid
      y2 = &centroid_pos_map(i, (j+1)%friend_num, 1);   // automatically loops around using %

      normalVectorBetween(*xn, *yn, *x1, *x2, *y1, *y2);    // calculate center coords between the two centroids.
                                                            // xc and yc automatically assigned the coords
    }
  }

  return 1;
};

// Function to find the area of each node's control volume. THe function loops
// over each pair of centroids belonging to a node, while summing the area of
// each triangle to find the total area of the hexagon or pentagon.
int Mesh::CalcControlVolumeArea(void)
{
  int i, j, node_num, f, friend_num;
  double * x1, * y1, * x2, * y2, * xc, * yc, *t_area, area;

  node_num = globals->node_num;

  for (i=0; i<node_num; i++)
  {

    f = node_friends(i,5);

    friend_num = 6;                                     // Assume hexagon (6 centroids)
    if (f == -1) {
      friend_num = 5;                                   // Check if pentagon (5 centroids)
    }

    xc = &node_pos_map(i,0,0);                          // set pointer to edge length array
    yc = &node_pos_map(i,0,1);

    t_area = &control_volume_surf_area_map(i);

    area = 0.0;
    for (j=0; j<friend_num; j++)                        // Loop through all centroids in the control volume
    {
      f = node_friends(i,j);

      x1 = &centroid_pos_map(i, j, 0);                  // set map coords for first centroid
      y1 = &centroid_pos_map(i, j, 1);

      x2 = &centroid_pos_map(i, (j+1)%friend_num, 0);   // set map coords for second centroid
      y2 = &centroid_pos_map(i, (j+1)%friend_num, 1);   // automatically loops around using %

      triangularArea(*t_area, *xc, *yc, *x1, *x2, *y1, *y2);     // calculate center coords between the two centroids.
                                                        // xc and yc automatically assigned the coords
      area += *t_area;
    }
    control_volume_surf_area_map(i) = area;
  }

  return 1;
}

// Function to read in text file containing mesh information
// The four pieces of information read and stored are:
//      1. Node ID numbers
//      2. Node positions in spherical coords
//      3. Neighbouring node ID numbers
//      4. All centroid positions in spherical coords
int Mesh::ReadMeshFile(void)
{
  std::string line, val;                 // strings for column and individual number
  std::string file_str;                  // string with path to mesh file.
  int i, node_id;

  file_str = globals->path + SEP + "input_files" + SEP + "grid_l" + std::to_string(globals->geodesic_l.Value()) + ".txt";

  // in stream for input.in file
  std::ifstream gridFile(file_str, std::ifstream::in);

  if (gridFile.is_open())
  {
    outstring << "Found mesh file: " + file_str << std::endl;
    globals->Output.Write(OUT_MESSAGE, &outstring);

    std::getline(gridFile, line);                               // READ HEADER
    while (std::getline(gridFile, line))
    {
      // std::cout<<line<<std::endl;
      std::istringstream line_ss(line);
      std::getline(line_ss >> std::ws,val,' ');                 // COL 0: Node ID number

      node_id = std::stoi(val);

      std::getline(line_ss >> std::ws,val,' ');                 // COL 1: Node Latitude
      node_pos_sph(node_id,0) = std::atof(val.c_str())*radConv;

      std::getline(line_ss >> std::ws,val,' ');                 // COL 2: Node Longitude
      node_pos_sph(node_id,1) = std::atof(val.c_str())*radConv;

      std::getline(line_ss >> std::ws,val,'{');                 // Read up to friends open bracket
      for (i = 0; i<5; i++)
      {
        std::getline(line_ss >> std::ws,val,',');               // Read each friend ID
        node_friends(node_id,i) = std::stoi(val);
      }
      std::getline(line_ss >> std::ws,val,'}');                 // Read last friend ID
      node_friends(node_id,5) = std::stoi(val);

      std::getline(line_ss >> std::ws,val,',');
      std::getline(line_ss >> std::ws,val,'{');                 // Read up to centroid list
      for (i = 0; i<5; i++)
      {
        std::getline(line_ss >> std::ws,val,'(');               // Read up to coord open bracket
        std::getline(line_ss >> std::ws,val,',');               // Read coord lat
        centroid_pos_sph(node_id,i,0) = std::atof(val.c_str())*radConv;

        std::getline(line_ss >> std::ws,val,')');               // Read coord lon
        centroid_pos_sph(node_id,i,1) = std::atof(val.c_str())*radConv;
        std::getline(line_ss >> std::ws,val,',');               // Read end of first coord
      }
      std::getline(line_ss >> std::ws,val,'(');                 // Read up to last coord open bracket
      std::getline(line_ss >> std::ws,val,',');                 // Read last coord lat
      centroid_pos_sph(node_id,5,0) = std::atof(val.c_str())*radConv;

      std::getline(line_ss >> std::ws,val,')');                 // Read last coord lon
      centroid_pos_sph(node_id,5,1) = std::atof(val.c_str())*radConv;

      std::getline(line_ss >> std::ws,val,'}');                 // Finish line


      // UNCOMMENT THIS BLOCK TO VERIFY THAT MESH FILE IS BEING READ CORRECTLY
      // std::cout<<node_id<<'\t';
      // std::cout << node_pos_sph(node_id, 0) <<'\t'<< node_pos_sph(node_id, 1)<<'\t';
      // for (i=0; i<6; i++) {
      //   std::cout<<node_friends(node_id, i)<<'\t';
      // }
      // for (i=0; i<6; i++) {
      //   std::cout<<'('<<centroid_pos_sph(node_id, i, 0)<<", "<<centroid_pos_sph(node_id, i, 1)<<")\t";
      // }
      // std::cout<<std::endl;
    }

    gridFile.close();
  }
  else
  {
    outstring << "ERROR: GRID FILE NOT FOUND AT " + file_str << std::endl;
    globals->Output.Write(ERR_MESSAGE, &outstring);
    globals->Output.TerminateODIS();
  }

  return 1;
};


// int Mesh::CalculateCellArea(void)
// {
//         cellArea = new double*[ReturnLatLen()];
//         for (int i=0; i<ReturnLatLen()-1; i++) {
//            cellArea[i] = new double[ReturnLonLen()];
//            for (int j=0; j<ReturnLonLen(); j++) {
//               cellArea[i][j] = globals->radius.Value()*globals->radius.Value();
//               cellArea[i][j] *= sin(radConv*lat[i]) - sin(radConv*lat[i+1]);
//               cellArea[i][j] *= radConv*dLon;
//            }
//         }
// };

// int Mesh::CalculateDt(void){
// 	double waveSpeed;
// 	double dt;
// 	double minDistance;
//
// 	outstring << std::endl << "Calculating maximum time step: " << std::endl;
//
// 	waveSpeed = sqrt(globals->g.Value()*globals->h.Value());
//
// 	outstring << std::endl << "\t\t Average gravitational wave speed = sqrt(gh): " << waveSpeed << " m/s."<<std::endl;
//
// 	minDistance = globals->radius.Value() * sin(dLat*radConv) * dLon*radConv;
//
// 	dt = minDistance/waveSpeed*1.5;
//
// 	outstring << "\t\t Time step calculated: " << dt <<" s."<< std::endl;
//
// 	if (dt<globals->timeStep.Value()) {
// 		outstring << std::endl << "\t\t WARNING: Calculated timestep is less than user selected dt = " << globals->timeStep.Value() <<" s." <<std::endl;
// 		outstring << "\t\t Setting dt to calculated value of " << dt <<" s." <<std::endl;
//
// 		globals->timeStep.SetValue(dt);
//
// 		globals->Output.Write(OUT_MESSAGE,&outstring);
//
// 		outstring << std::endl << "WARNING: Using maximum calculated timestep as user defined dt is too large." << std::endl;
// 		globals->Output.Write(ERR_MESSAGE,&outstring);
// 	}
// 	else {
// 		outstring << std::endl << "\t\t Calculated timestep is greater than user selected dt = " << globals->timeStep.Value() <<" s." <<std::endl;
// 		outstring << "\t\t ODIS will keep user defined timestep. You could save time if you chose a larger timestep."<<std::endl;
//
// 		globals->Output.Write(OUT_MESSAGE,&outstring);
// 	}
//
//
//
//
//
// 	return 1;
// }

//
// int Mesh::InitializeArrays(void)
// {
//   int node_num, i;
//
//   // geodesic node num expression (Lee and Macdonald, 2008)
//   node_num = 10 * pow(pow(2, globals->geodesic_l.Value() - 1), 2) + 2;
//
//   // Initialize node position array
//
//   // node_pos_sph = new double*[node_num];
//   // node_pos_sph[0] = new double[node_num * 2];
//   // for (i=1; i < node_num; i++) {
//   //   node_pos_sph[i] = &node_pos_sph[0][i*2];
//   // }
//
//   //
//   // node_pos_sph = new double*[node_num];
//   // // node_pos_sph[0] = new double[node_num * 2];
//   // for (i=0; i < node_num; i++) {
//   //   node_pos_sph[i] = new double[2];
//   // }
//
//
//   // Initialize node position array
//   // node_pos_map = new double*[node_num];
//   // node_pos_map[0] = new double[node_num * 2];
//   // for (i=1; i < node_num; i++) {
//   //   node_pos_map[i] = &node_pos_map[0][i*2];
//   // }
//
//   // Initialize node friends array. Each row in node_friends contains the ID (index)
//   // to all surrounding friends. E.g, the ID's of the neighbouring nodes to node 4
//   // would be accessed via node_friends[4][0], node_friends[4][1], etc.
//   // node_friends = new int*[node_num];
//   // node_friends[0] = new int[node_num * 6];
//   // for (i=1; i < node_num; i++) {
//   //   node_friends[i] = &node_friends[0][i*6];
//   // }
//
//   return 1;
// };
