/* File: mesh.cpp

   Contains all member functions for the Mesh class. Mesh exists to store and
   calculate information about the numerical grid. Primarily, quantities relevant
   to the spatial operators (gradient, divergence, laplacian, and curl) are
   calculated here based on the grid that is loaded into memory. Other useful
   properties that are stored are the node map, i.e., which nodes are linked to
   each other. Additionally, analytical functions that vary with space are also
   calculated and stored here, such as cos, sin, and Legendre functions.

*/

#include "mesh.h"
#include "globals.h"
#include "math.h"
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include "mathRoutines.h"

#include <mkl.h>
#include <mkl_spblas.h>
#include <math.h>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cstring>

#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

Mesh::Mesh(Globals * Globals, int N, int N_ll, int l_max)
   :node_pos_sph(N,2),
    node_pos_map(N,7,2),        // len 7 to include central node coords (even though it's zero)
    node_friends(N,6),
    node_dists(N,6),
    centroid_pos_sph(N,6,2),
    centroid_pos_map(N,6,2),
    node_m(N,7),                // len 7 to include central (parent) node
    centroid_m(N,6),            // len 6 as there is are only 5 or 6 centroids
    node_vel_trans(N,7,2),      // In the last dimension, element 1 is cos_a, element 2 is sin_a
    control_vol_edge_len(N,6),
    control_vol_edge_centre_pos_map(N,6,2),
    control_vol_edge_centre_m(N,6),
    control_vol_edge_normal_map(N,6,2),
    control_volume_surf_area_map(N),
    control_volume_mass(N),
    node_friend_element_areas_map(N,6,3),
    centroid_node_dists_map(N,6),
    grad_coeffs(N, 7, 2),
    div_coeffs(N, 7, 2),

    trigLat(N,2),
    trigLon(N,2),
    trig2Lat(N,2),
    trig2Lon(N,2),
    trigSqLat(N,2),
    trigSqLon(N,2),

    Pbar_lm(N, l_max+1, l_max+1),
    Pbar_lm_deriv(N, l_max+1, l_max+1),

    Pbar_cosMLon(N, l_max+1, l_max+1),
    Pbar_sinMLon(N, l_max+1, l_max+1),
    Pbar_deriv_cosMLon(N, l_max+1, l_max+1),
    Pbar_deriv_sinMLon(N, l_max+1, l_max+1),

    sh_matrix(N, (l_max+1)*(l_max+2) - 6), //-6 is to ignore degrees 0 and 1
    sh_matrix_fort(N*((l_max+1)*(l_max+2) - 6)),

    trigMLon(N, l_max+1, 2)
{

    globals = Globals;                   // define reference to all constants

    node_num = globals->node_num;

    // Read in grid file
    ReadMeshFile();

    // Calculate mapping coordinates
    CalcMappingCoords();

    CalcNodeDists();

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

    CalcControlVolumeMass();

    // Find the area of each sub triangle within each element
    CalcElementAreas();

    CalcCentNodeDists();

    // Evaluate trig functions at every node
    CalcTrigFunctions();

    CalcMaxTimeStep();

    CalcLegendreFuncs();

    CalcGradOperatorCoeffs();

    CalcDivOperatorCoeffs();

    GeneratePressureSolver();

    ReadWeightingFile();

};

// Function to calculate the cosine(alpha) and sine(alpha) velocity tranform
// factors from Lee and Macdonald (2009). These factors are constant in time,
// and so are pre-calculated here to be easily accessed during the spatial
// integration.
int Mesh::CalcVelocityTransformFactors(void)
{
    int i, j, f;
    double * cos_a, * sin_a;
    double lat1, lat2, lon1, lon2;

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

    return 1;
};

// Function to convert all spherical coorinate quantities to mapping coordinates
// x and y from Lee and Macdonald (2009). Mapping factors are also calculated
// and stored. These quantities are stored for node positions, centroid positions,
// node mapping factors, and centroid mapping factors.
int Mesh::CalcMappingCoords(void)
{
    int i, j, f;
    double * x, * y, * m, r;
    double lat1, lat2, lon1, lon2;

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
    int i, j, f, friend_num;
    double * x1, * y1, * x2, * y2, * edge_len;

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

// Function to find the spherical distance between a node and each of
// its neighbours.
int Mesh::CalcNodeDists(void)
{
    int i, j, f, friend_num;
    double * lat1, * lat2, * lon1, * lon2, * edge_len;
    double r;

    r = globals->radius.Value();

    for (i=0; i<node_num; i++)
    {

        f = node_friends(i,5);

        friend_num = 6;                                     // Assume hexagon (6 centroids)
        if (f == -1) {
            friend_num = 5;                                 // Check if pentagon (5 centroids)
            node_dists(i,5) = -1.0;
        }

        lat1 = &node_pos_sph(i, 0);                     // set map coords for first
        lon1 = &node_pos_sph(i, 1);
        for (j=0; j<friend_num; j++)                        // Loop through all centroids in the control volume
        {
            edge_len = &node_dists(i,j);                    // set pointer to edge length array

            f = node_friends(i,j);

            lat2 = &node_pos_sph(f, 0);                     // set map coords for second centroid
            lon2 = &node_pos_sph(f, 1);

            distanceBetweenSph(*edge_len, *lat1, *lat2, *lon1, *lon2, r);   // calculate distance between the two centroids.
            // Edge_len automatically assigned the length
        }
    }


    return 1;
};


int Mesh::CalcMaxTimeStep(void)
{
  int i, j, f, friend_num;
  double dist;
  double g, h_max, dt;

  g = globals->g.Value();
  h_max = globals->h.Value();


  dt = globals->timeStep.Value();

  if (globals->surface_type == FREE ||
      globals->surface_type == FREE_LOADING ||
      globals->surface_type == LID_LOVE)
  {
      if (globals->surface_type == LID_LOVE) g *= -globals->shell_factor_beta[globals->l_max.Value()];
      std::cout<<g<<std::endl;
      for (i=0; i<node_num; i++)
      {
          f = node_friends(i,5);
          friend_num = 6;                                     // Assume hexagon (6 centroids)
          if (f == -1) {
              friend_num = 5;                                   // Check if pentagon (5 centroids)
          }
          for (j=0; j<friend_num; j++)                      // Loop through all centroids in the control volume
          {
              dist = node_dists(i,j) * 0.5;                 // consider half node-node distance
              dt = std::min(dt, dist/sqrt(g*h_max));
          }
      }

      dt *= 1.0;         // take some caution
  }

  std::cout<<"DT: "<<dt<<std::endl;

  globals->timeStep.SetValue(dt);



  return 1;
}

// Function to find the length of control volume edges for each node, in mapping
// coordinates. Values are stored in control_vol_edge_len 2D array.
int Mesh::CalcCentNodeDists(void)
{
    int i, j, f, friend_num;
    double * x1, * y1, * x2, * y2, * edge_len;

    for (i=0; i<node_num; i++)
    {

        f = node_friends(i,5);

        friend_num = 6;                                     // Assume hexagon (6 centroids)
        if (f == -1) {
            friend_num = 5;                                   // Check if pentagon (5 centroids)
            centroid_node_dists_map(i,5) = -1.0;
        }

        for (j=0; j<friend_num; j++)                        // Loop through all centroids in the control volume
        {
            f = node_friends(i,j);

            edge_len = &centroid_node_dists_map(i,j);            // set pointer to edge length array

            x1 = &node_pos_map(i, 0, 0);   // set map coords for second centroid
            y1 = &node_pos_map(i, 0, 1);   // automatically loops around using %

            x2 = &centroid_pos_map(i, j%friend_num, 0);   // set map coords for second centroid
            y2 = &centroid_pos_map(i, j%friend_num, 1);   // automatically loops around using %

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
    int i, j, f, friend_num;
    double * x1, * y1, * x2, * y2, * xc, * yc, * m;
    double lat1, lat2, lat3, lon1, lon2, lon3;
    double r = globals->radius.Value();

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

            // lat1 = centroid_pos_sph(i, j, 0);                 // get first centroid coords
            // lon1 = centroid_pos_sph(i, j, 1);

            lat1 = node_pos_sph(i, 0);                 // get first centroid coords
            lon1 = node_pos_sph(i, 1);

            lat2 = centroid_pos_sph(i, j, 0);  // get second centroid coords
            lon2 = centroid_pos_sph(i, j, 1);

            // lat3 = centroid_pos_sph(i, (j+1)%friend_num, 0);  // get second centroid coords
            // lon3 = centroid_pos_sph(i, (j+1)%friend_num, 1);
            //
            // lat2 = 0.5*(lat2 + lat3);
            // lon2 = 0.5*(lon2 + lon3);

            mapFactorAtPoint(*m, lat1, lat2, lon1, lon2);     // calculate map factor and store


            // Calculate here as the midpoint to the cv edge is only known in
            // mapped coordinates
            *m = (4.*pow(r,2.0) + pow(*xc, 2.0) + pow(*yc, 2.0))/(4.*pow(r,2.0));

            // std::cout<<*m<<'\t'<<(4.*r*r + (*xc)*(*xc)+ (*yc)*(*yc))/(4.*r*r)<<std::endl;
        }
    }

    return 1;
};

// Function to find the unit normal vector to the control volume edge, in
// mapping coordinates. Values are stored in control_vol_edge_normal_map 2Dd
// array.
int Mesh::CalcControlVolumeEdgeNormals(void)
{
    int i, j, f, friend_num;
    double * x1, * y1, * x2, * y2, * xn, * yn;

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

            xn = &control_vol_edge_normal_map(i,j,0);         // set pointer to edge length array
            yn = &control_vol_edge_normal_map(i,j,1);

            x1 = &centroid_pos_map(i, j, 0);                  // set map coords for first centroid
            y1 = &centroid_pos_map(i, j, 1);

            x2 = &centroid_pos_map(i, (j+1)%friend_num, 0);   // set map coords for second centroid
            y2 = &centroid_pos_map(i, (j+1)%friend_num, 1);   // automatically loops around using %

            normalVectorBetween(*xn, *yn, *x1, *x2, *y1, *y2);    // calculate center coords between the two centroids.
            // xc and yc automatically assigned the coords

            // std::cout<<i<<"   ("<<*x1<<", "<<*y1<<"), "<<'('<<*x2<<", "<<*y2<<"), "<<'('<<*xn<<", "<<*yn<<"), "<<std::endl;

        }
    }

    return 1;
};

// Function to find the area of each node's control volume. THe function loops
// over each pair of centroids belonging to a node, while summing the area of
// each triangle to find the total area of the hexagon or pentagon.
int Mesh::CalcControlVolumeArea(void)
{
    int i, j, f, friend_num;
    double * x1, * y1, * x2, * y2, * xc, * yc, *t_area, area;

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
};


int Mesh::CalcControlVolumeMass(void)
{
  int i, j, k, as, ae, f, friend_num;
  double * lat1, * lon1, * lat2, * lon2, * lat3, * lon3, *t_area, r;
  double area, avg_area;
  double ax, ay, az;
  double bx, by, bz;
  double cx, cy, cz;
  double * vol;
  double vol1, vol2;

  // r = globals->radius.Value();//# + globals->shell_thickness.Value();
  // r_core = r - globals->h.Value();
  // r_ocean = r;
  // double mass_sum = 0.0;
  //
  // avg_area = 0.0;
  // for (i=0; i<node_num; i++)
  // {
  //
  //     f = node_friends(i,5);
  //
  //     friend_num = 6;                    // Assume hexagon (6 centroids)
  //     if (f == -1) {
  //         friend_num = 5;                // Check if pentagon (5 centroids)
  //     }
  //
  //     lat1 = &node_pos_sph(i,0);
  //     lon1 = &node_pos_sph(i,1);
  //
  //     vol = &control_volume_mass(i);     // set pointer to area of sub element
  //     *vol = 0.0;
  //     for (j=0; j<friend_num; j++)                       // Loop through all centroids in the control volume
  //     {
  //
  //         lat2 = &centroid_pos_sph(i,j%friend_num, 0);
  //         lon2 = &centroid_pos_sph(i,j%friend_num, 1);
  //
  //         lat3 = &centroid_pos_sph(i,(j+1)%friend_num, 0);
  //         lon3 = &centroid_pos_sph(i,(j+1)%friend_num, 1);
  //
  //         sph2cart(ax, ay, az, r_ocean, radConv*0.5 - *lat1, *lon1);
  //         sph2cart(bx, by, bz, r_ocean, radConv*0.5 - *lat2, *lon2);
  //         sph2cart(cx, cy, cz, r_ocean, radConv*0.5 - *lat3, *lon3);
  //
  //         volumeSphericalTriangle(vol1, ax, bx, cx, ay, by, cy, az, bz, cz, r_ocean);
  //
  //         sph2cart(ax, ay, az, r_core, radConv*0.5 - *lat1, *lon1);
  //         sph2cart(bx, by, bz, r_core, radConv*0.5 - *lat2, *lon2);
  //         sph2cart(cx, cy, cz, r_core, radConv*0.5 - *lat3, *lon3);
  //
  //         volumeSphericalTriangle(vol2, ax, bx, cx, ay, by, cy, az, bz, cz, r_core);
  //
  //         *vol += fabs(vol1)-fabs(vol2);
  //     }
  //     *vol *= 1000.0;
  // }
  //
  // for (i=0; i<node_num; i++)
  // {
  //   vol = &control_volume_mass(i);
  //   mass_sum += *vol;
  // }
  // std::cout<<"MASS: "<<mass_sum<<std::endl;
  // return 1;

  r = globals->radius.Value();//# + globals->shell_thickness.Value();

  double mass_sum = 0.0;

  avg_area = 0.0;
  for (i=0; i<node_num; i++)
  {

      f = node_friends(i,5);

      friend_num = 6;                    // Assume hexagon (6 centroids)
      if (f == -1) {
          friend_num = 5;                // Check if pentagon (5 centroids)
      }

      lat1 = &node_pos_sph(i,0);
      lon1 = &node_pos_sph(i,1);

      t_area = &control_volume_mass(i);
      *t_area = 0.0;
      for (j=0; j<friend_num; j++)                       // Loop through all centroids in the control volume
      {

          lat2 = &centroid_pos_sph(i,j%friend_num, 0);
          lon2 = &centroid_pos_sph(i,j%friend_num, 1);

          lat3 = &centroid_pos_sph(i,(j+1)%friend_num, 0);
          lon3 = &centroid_pos_sph(i,(j+1)%friend_num, 1);


          triangularAreaSph(area, *lat1, *lat2, *lat3, *lon1, *lon2, *lon3, r);     // calculate subelement area
          *t_area += area;

      }
      avg_area += *t_area;
  }
  avg_area /= node_num;
  for (i=0; i<node_num; i++)
  {
     t_area = &control_volume_mass(i);
    *t_area *= avg_area/(*t_area) * 1000.0 * globals->h.Value();

    mass_sum += *t_area;
  }

  return 1;
}

// Function to fine the areas of each subelement. Each subelement is a triangle
// subtended by the central node, one friend, and the centroid belonging to the
// friend. Values are stored in node_friend_element_areas_map
int Mesh::CalcElementAreas(void)
{
    int i, j, k, as, ae, f, friend_num;
    double * x1, * y1, * x2, * y2, * xc, * yc, *t_area;

    for (i=0; i<node_num; i++)
    {

        f = node_friends(i,5);

        friend_num = 6;                    // Assume hexagon (6 centroids)
        if (f == -1) {
            friend_num = 5;                // Check if pentagon (5 centroids)
        }
        for (j=0; j<friend_num; j++)                       // Loop through all centroids in the control volume
        {
            xc = &centroid_pos_map(i,j,0);                 // set pointer element centroid
            yc = &centroid_pos_map(i,j,1);

            f = node_friends(i,j);

            for (k=0; k < 3; k++)
            {
                t_area = &node_friend_element_areas_map(i, j, k);     // set pointer to area of sub element


                // messy code TODO - try and tidy
                as = j+k;
                ae = (j+(k+1)%3)%(friend_num+2);
                if (k==0) as = 0;
                else if (k==2) {
                    ae = 0;
                }

                if (j == friend_num-1) {
                    if (k==1) ae = 1;
                    else if (k==2) as = 1;
                }

                x1 = &node_pos_map(i, as, 0);   // set map coords for first vertex of the element
                y1 = &node_pos_map(i, as, 1);

                x2 = &node_pos_map(i, ae, 0);   // set map coords for second vertex on the element
                y2 = &node_pos_map(i, ae, 1);

                triangularArea(*t_area, *xc, *yc, *x1, *x2, *y1, *y2);     // calculate subelement area

            }
        }
    }

    return 1;
};

// Function to evaluate trigonometric functions over latitude and longitude for
// every node.
int Mesh::CalcTrigFunctions(void)
{
    int i;
    double lat, lon;

    for (i=0; i<node_num; i++)
    {
        lat = node_pos_sph(i, 0);
        lon = node_pos_sph(i, 1);

        trigLat(i,0) = cos(lat);
        trigLat(i,1) = sin(lat);

        trigLon(i,0) = cos(lon);
        trigLon(i,1) = sin(lon);

        trig2Lat(i,0) = cos(2.0 * lat);
        trig2Lat(i,1) = sin(2.0 * lat);

        trig2Lon(i,0) = cos(2.0 * lon);
        trig2Lon(i,1) = sin(2.0 * lon);

        trigSqLat(i,0) = pow(cos(lat), 2.0);
        trigSqLat(i,1) = pow(sin(lat), 2.0);

        trigSqLon(i,0) = pow(cos(lon), 2.0);
        trigSqLon(i,1) = pow(sin(lon), 2.0);
    }

    return 1;
};

int Mesh::CalcLegendreFuncs(void)
{
    int i, l, m, l_max;
    double cosCoLat;

    Array2D<double> * temp_legendre;    // temp array for legendre polynomials
    Array2D<double> * temp_dlegendre;   // temp array for legendre polynomial derivs

    //-----------------------------
    l_max = globals->l_max.Value();

    temp_legendre = new Array2D<double>(2 * (l_max+1), 2 * (l_max+1));
    temp_dlegendre = new Array2D<double>(2 * (l_max+1), 2 * (l_max+1));

    for (i=0; i<node_num; i++)
    {
        cosCoLat = cos(pi*0.5 - node_pos_sph(i,0));     // cos(colatitude) of node i

        getLegendreFuncs(cosCoLat, *temp_legendre, l_max);
        getLegendreFuncsDeriv(cosCoLat, *temp_dlegendre, l_max);

        for (l = 0; l < l_max+1; l++)
        {
            for (m = 0; m <= l; m++)
            {

                // assign legendre pol at node i for degree l and order m
                Pbar_lm(i, l, m) = (*temp_legendre)(l, m);

                // assign legendre pol derivative at node i for degree l and order m
                Pbar_lm_deriv(i, l, m) = (*temp_dlegendre)(l, m)*sin(pi*0.5 - node_pos_sph(i, 0));

            }
        }

        // calculate cos(m*longitude) and sin(m*longitude) for node i
        for (m=0; m < l_max+1; m++)
        {
            trigMLon(i, m, 0) = cos( (double)m * node_pos_sph( i, 1 ) );
            trigMLon(i, m, 1) = sin( (double)m * node_pos_sph( i, 1 ) );
        }

        // int count = 0;
        for (l = 0; l < l_max+1; l++)
        {
            for (m = 0; m <= l; m++)
            {
                Pbar_cosMLon(i, l, m) = Pbar_lm(i, l, m)*trigMLon(i, m, 0);
                Pbar_sinMLon(i, l, m) = Pbar_lm(i, l, m)*trigMLon(i, m, 1);

                Pbar_deriv_cosMLon(i, l, m) = Pbar_lm_deriv(i, l, m)*trigMLon(i, m, 0);
                Pbar_deriv_sinMLon(i, l, m) = Pbar_lm_deriv(i, l, m)*trigMLon(i, m, 1);

            }
        }

        int count = 0;
        for (l = 2; l < l_max+1; l++)
        {
            for (m = 0; m <= l; m++)
            {
                sh_matrix(i, count) = Pbar_cosMLon(i, l, m);
                if (globals->surface_type == FREE_LOADING) sh_matrix(i, count) *= globals->loading_factor[l];
                else if (globals->surface_type == LID_MEMBR ||
                         globals->surface_type == LID_LOVE) sh_matrix(i, count) *= globals->shell_factor_beta[l];


                count++;

                sh_matrix(i, count) = Pbar_sinMLon(i, l, m);
                if (globals->surface_type == FREE_LOADING) sh_matrix(i, count) *= globals->loading_factor[l];
                else if (globals->surface_type == LID_MEMBR ||
                         globals->surface_type == LID_LOVE) sh_matrix(i, count) *= globals->shell_factor_beta[l];

                count++;
            }
        }
    }

    int count = 0;
    for (m = 0; m <= (l_max+1)*(l_max+2) - 6; m++)
    {
        for (i=0; i<node_num; i++)
        {
            sh_matrix_fort(count) = sh_matrix(i, m);
            count++;
        }
    }


    // free up memory!
    delete temp_legendre;
    delete temp_dlegendre;

    return 1;
}


int Mesh::CalcGradOperatorCoeffs(void)
{
    int i, j, j2, j3, f, f_num;
    double areaCV, areaElement, *coeff_x, *coeff_y;
    double alpha, beta, gamma, B_x, B_y;

    for (i=0; i<node_num; i++)
    {

        f_num = 6;                                     // Assume hexagon (6 centroids)
        f = node_friends(i,5);
        if (f == -1) {
            f_num = 5;                                   // Check if pentagon (5 centroids)
        }

        // Control volume area
        areaCV = control_volume_surf_area_map(i);

        // Calculate the central node coefficient first
        coeff_x = &grad_coeffs(i, 0, 0);
        coeff_y = &grad_coeffs(i, 0, 1);
        *coeff_x = 0.0;
        *coeff_y = 0.0;
        for (j=0; j<f_num; j++)
        {
            j2 = (j-1)%f_num;
            if (j2 < 0) j2 += f_num;

            alpha = node_friend_element_areas_map(i, j, 1);

            B_x = control_vol_edge_centre_m(i, j) * control_vol_edge_normal_map(i, j, 0) * control_vol_edge_len(i, j);
            B_y = control_vol_edge_centre_m(i, j) * control_vol_edge_normal_map(i, j, 1) * control_vol_edge_len(i, j);

            B_x += control_vol_edge_centre_m(i, j2) * control_vol_edge_normal_map(i, j2, 0) * control_vol_edge_len(i, j2);
            B_y += control_vol_edge_centre_m(i, j2) * control_vol_edge_normal_map(i, j2, 1) * control_vol_edge_len(i, j2);

            // Total element area = alpha + beta + gamma
            areaElement = alpha + node_friend_element_areas_map(i, j, 0) + node_friend_element_areas_map(i, j, 2);

            *coeff_x += B_x/areaElement * alpha;
            *coeff_y += B_y/areaElement * alpha;

        }
        *coeff_x *= 0.5 * 1./areaCV;
        *coeff_y *= 0.5 * 1./areaCV;

        // Now calculate coeffs for neighbouring nodes (matrix diagonals)
        for (j=0; j<f_num; j++)
        {
            // First term ------------------------------------------------------
            coeff_x = &grad_coeffs(i, j+1, 0);
            coeff_y = &grad_coeffs(i, j+1, 1);

            j2 = (j-1)%f_num;
            if (j2 < 0) j2 += f_num;

            beta = node_friend_element_areas_map(i, j, 2);

            B_x = control_vol_edge_centre_m(i, j) * control_vol_edge_normal_map(i, j, 0) * control_vol_edge_len(i, j);
            B_y = control_vol_edge_centre_m(i, j) * control_vol_edge_normal_map(i, j, 1) * control_vol_edge_len(i, j);

            B_x += control_vol_edge_centre_m(i, j2) * control_vol_edge_normal_map(i, j2, 0) * control_vol_edge_len(i, j2);
            B_y += control_vol_edge_centre_m(i, j2) * control_vol_edge_normal_map(i, j2, 1) * control_vol_edge_len(i, j2);

            // Total element area = alpha + beta + gamma
            areaElement = node_friend_element_areas_map(i, j, 0) + beta + node_friend_element_areas_map(i, j, 1);

            *coeff_x = B_x/areaElement * beta;
            *coeff_y = B_y/areaElement * beta;

            // Second term -----------------------------------------------------
            gamma = node_friend_element_areas_map(i, j2, 0);

            j3 = (j-2)%f_num;
            if (j3 < 0) j3 += f_num;

            B_x = control_vol_edge_centre_m(i, j2) * control_vol_edge_normal_map(i, j2, 0) * control_vol_edge_len(i, j2);
            B_y = control_vol_edge_centre_m(i, j2) * control_vol_edge_normal_map(i, j2, 1) * control_vol_edge_len(i, j2);

            B_x += control_vol_edge_centre_m(i, j3) * control_vol_edge_normal_map(i, j3, 0) * control_vol_edge_len(i, j3);
            B_y += control_vol_edge_centre_m(i, j3) * control_vol_edge_normal_map(i, j3, 1) * control_vol_edge_len(i, j3);


            // Total element area = alpha + beta + gamma
            areaElement = node_friend_element_areas_map(i, j2, 2) + node_friend_element_areas_map(i, j2, 1) + gamma;

            *coeff_x += B_x/areaElement * gamma;
            *coeff_y += B_y/areaElement * gamma;

            *coeff_x *= 0.5 * 1./areaCV;
            *coeff_y *= 0.5 * 1./areaCV;


        }

    }

    return 1;
}

int Mesh::CalcDivOperatorCoeffs(void)
{
    int i, j, j2, j3, f, f_num;
    double areaCV, areaElement, *coeff_x, *coeff_y;
    double alpha, beta, gamma, C_x, C_y;

    for (i=0; i<node_num; i++)
    {

        f_num = 6;                                     // Assume hexagon (6 centroids)
        f = node_friends(i,5);
        if (f == -1) {
            f_num = 5;                                   // Check if pentagon (5 centroids)
        }

        // Control volume area
        areaCV = control_volume_surf_area_map(i);

        // Calculate the central node coefficient first
        coeff_x = &div_coeffs(i, 0, 0);
        coeff_y = &div_coeffs(i, 0, 1);
        *coeff_x = 0.0;
        *coeff_y = 0.0;
        for (j=0; j<f_num; j++)
        {
            j2 = (j-1)%f_num;
            if (j2 < 0) j2 += f_num;

            alpha = node_friend_element_areas_map(i, j, 1);

            C_x = control_vol_edge_normal_map(i, j, 0) * control_vol_edge_len(i, j) /  control_vol_edge_centre_m(i, j);
            C_y = control_vol_edge_normal_map(i, j, 1) * control_vol_edge_len(i, j) /  control_vol_edge_centre_m(i, j);

            C_x += control_vol_edge_normal_map(i, j2, 0) * control_vol_edge_len(i, j2) / control_vol_edge_centre_m(i, j2);
            C_y += control_vol_edge_normal_map(i, j2, 1) * control_vol_edge_len(i, j2) / control_vol_edge_centre_m(i, j2);

            // Total element area = alpha + beta + gamma
            areaElement = alpha + node_friend_element_areas_map(i, j, 0) + node_friend_element_areas_map(i, j, 2);

            *coeff_x += C_x/areaElement * alpha;
            *coeff_y += C_y/areaElement * alpha;

        }
        *coeff_x *= 0.5 * 1./areaCV;
        *coeff_y *= 0.5 * 1./areaCV;

        // Now calculate coeffs for neighbouring nodes (matrix diagonals)
        for (j=0; j<f_num; j++)
        {
            // First term ------------------------------------------------------
            coeff_x = &div_coeffs(i, j+1, 0);
            coeff_y = &div_coeffs(i, j+1, 1);

            j2 = (j-1)%f_num;
            if (j2 < 0) j2 += f_num;

            beta = node_friend_element_areas_map(i, j, 2);

            C_x = control_vol_edge_normal_map(i, j, 0) * control_vol_edge_len(i, j) /  control_vol_edge_centre_m(i, j);
            C_y = control_vol_edge_normal_map(i, j, 1) * control_vol_edge_len(i, j) /  control_vol_edge_centre_m(i, j);

            C_x += control_vol_edge_normal_map(i, j2, 0) * control_vol_edge_len(i, j2) / control_vol_edge_centre_m(i, j2);
            C_y += control_vol_edge_normal_map(i, j2, 1) * control_vol_edge_len(i, j2) / control_vol_edge_centre_m(i, j2);

            // Total element area = alpha + beta + gamma
            areaElement = node_friend_element_areas_map(i, j, 0) + beta + node_friend_element_areas_map(i, j, 1);

            *coeff_x = C_x/areaElement * beta;
            *coeff_y = C_y/areaElement * beta;

            // Second term -----------------------------------------------------
            gamma = node_friend_element_areas_map(i, j2, 0);

            j3 = (j-2)%f_num;
            if (j3 < 0) j3 += f_num;

            C_x = control_vol_edge_normal_map(i, j2, 0) * control_vol_edge_len(i, j2) / control_vol_edge_centre_m(i, j2);
            C_y = control_vol_edge_normal_map(i, j2, 1) * control_vol_edge_len(i, j2) / control_vol_edge_centre_m(i, j2);

            C_x += control_vol_edge_normal_map(i, j3, 0) * control_vol_edge_len(i, j3) / control_vol_edge_centre_m(i, j3);
            C_y += control_vol_edge_normal_map(i, j3, 1) * control_vol_edge_len(i, j3) / control_vol_edge_centre_m(i, j3);


            // Total element area = alpha + beta + gamma
            areaElement = node_friend_element_areas_map(i, j2, 2) + node_friend_element_areas_map(i, j2, 1) + gamma;

            *coeff_x += C_x/areaElement * gamma;
            *coeff_y += C_y/areaElement * gamma;

            *coeff_x *= 0.5 * 1./areaCV;
            *coeff_y *= 0.5 * 1./areaCV;
        }
    }

    return 1;
}

// Function to generate a sparse matrix A for the Laplacian operator, convert it
// to CSR format, and create a solver object that is capable of solving the
// linear system A*p = d to find the pressure field, p, give a rhs vector d.
int Mesh::GeneratePressureSolver(void)
{
    int i, j, j2, j3, f, f_num;
    double areaCV, *coeff;
    double D;


    double * hVar = new double[node_num];
    double a_ratio = 0.0;
    for (i=0; i<node_num; i++) hVar[i] = 1.0;//globals->h.Value()*(1. + a_ratio*trigLat(i,0)*pow(trigMLon(i,4,0),4.0));


    // }

    Array2D<double> * pressure_matrix;
    pressure_matrix = new Array2D<double>(node_num, 7);

    // why is this line here? Where would it be better placed?
    mkl_set_num_threads(1);

    //-------------- Calculate the laplacian matrix coefficients ---------------

    for (i=0; i<node_num; i++)
    {

        f_num = 6;                  // Assume hexagon (6 centroids)
        f = node_friends(i,5);
        if (f == -1) f_num = 5;     // Check if pentagon (5 centroids)

        // Control volume area
        areaCV = control_volume_surf_area_map(i);

        // Calculate the central node coefficient first
        coeff = &(*pressure_matrix)(i, 0);
        *coeff = 0.0;
        for (j=0; j<f_num; j++)
        {
            j2 = (j-1)%f_num;
            if (j2 < 0) j2 += f_num;
            f = node_friends(i, j);

            //          this part is the spatially
            //          varying scalar, TODO change
            //          to ocean thickness!!!!!
            //        -----------------------------
            D =  0.5*(hVar[i] + hVar[f]) * control_vol_edge_len(i, j2) /  ( control_vol_edge_centre_m(i, j2) * node_dists(i, j) );
            // 0.5*(trigLon(i, 1) + trigLon(f, 1)) *
            *coeff += D;
        }
        *coeff *= -1./areaCV; // Need the negative sign here!

        // Now calculate coeffs for neighbouring nodes (matrix diagonals)
        for (j=0; j<f_num; j++)
        {
            coeff = &(*pressure_matrix)(i, j+1);

            j2 = (j-1)%f_num;
            if (j2 < 0) j2 += f_num;
            f = node_friends(i,j);


            D = 0.5*(hVar[i] + hVar[f]) * control_vol_edge_len(i, j2) /  ( control_vol_edge_centre_m(i, j2) * node_dists(i, j) );
            // 0.5*(trigLon(i, 1) + trigLon(f, 1)) *

            *coeff = D;
            *coeff *= 1./areaCV; // No negative sign here!
        }
    }


    //--------------- Construct sparse matrix in CSR format --------------------

    //--------------- Define objects to construct matrix -----------------------

    double * nzCoeffs_grad_x;   // non-zero grad operator coeffis in x direction
    double * nzCoeffs_grad_y;   // non-zero grad operator coeffis in y direction
    double * nzCoeffs_div_x;    // non-zero div operator coeffis in x direction
    double * nzCoeffs_div_y;    // non-zero div operator coeffis in y direction

    // Here we assign the number of non-zero matrix elements for each operator.
    // There are 12 pentagons on the geodesic grid, each with 5 neighbours and
    // one central node (total of 6). The other cells are hexagons with 6
    // neighbours and one central node (total of 7). Therefore, the total number
    // of non-zero coefficients in our sparse matrix is given as:
    int nNonZero = (node_num-12)*7 + 12*6;

    int     * colIndx;    // column index for each non-zero element
    int     * rowIndx;    // row indices for CSR matrix format
    double * nzCoeffs;    // array for each non-zero coefficient

    int         error;    // error message from MKL

    int * rowStart, * rowEnd; // row indexes for CSR sparse matrix format


    //---------------- Initialise objects into memory --------------------------

    rowStart = new int[node_num];
    rowEnd   = new int[node_num];

    colIndx  = new int[nNonZero];
    nzCoeffs = new double[nNonZero];
    rowIndx  = new int[node_num + 1];

    nzCoeffs_grad_x = new double[nNonZero];
    nzCoeffs_grad_y = new double[nNonZero];
    nzCoeffs_div_x  = new double[nNonZero];
    nzCoeffs_div_y  = new double[nNonZero];

    // assign the non-zero coefficients and their column indexes in CSR format
    int count=0;
    for (i=0; i<node_num; i++)
    {
        f_num = 6;               // Assume hexagon (6 neighbours)
        f = node_friends(i,5);
        if (f == -1) f_num = 5;  // Check if pentagon (5 neighbours)

        nzCoeffs[count] = (*pressure_matrix)(i, 0);  // assign diagonal element

        nzCoeffs_grad_x[count] = grad_coeffs(i, 0, 0);
        nzCoeffs_grad_y[count] = grad_coeffs(i, 0, 1);

        nzCoeffs_div_x[count]  = hVar[i] * div_coeffs(i, 0, 0);
        nzCoeffs_div_y[count]  = hVar[i] * div_coeffs(i, 0, 1);

        colIndx[count] = i;

        count++;
        double cos_a, sin_a;
        for (j=0; j<f_num; j++)                       // assign off-diagonals
        {
            nzCoeffs[count] = (*pressure_matrix)(i, j+1);
            colIndx[count]  = node_friends(i, j);

            cos_a = node_vel_trans(i, j+1, 0);
            sin_a = node_vel_trans(i, j+1, 1);

            nzCoeffs_grad_x[count] = grad_coeffs(i, j+1, 0);
            nzCoeffs_grad_y[count] = grad_coeffs(i, j+1, 1);

            f = node_friends(i,j);

            // Because the gradient operator returns a vector, we have to transform
            // the vector into mapped coordinates. We do this here by premultiplying
            // the components of the divergence operator by the transpose of the
            // the transform matrix defined in node_vel_trans.
            // nzCoeffs_div_x[count]  =  cos_a*div_coeffs(i, j+1, 0) - sin_a*div_coeffs(i,  j+1, 1);
            // nzCoeffs_div_y[count]  =  sin_a*div_coeffs(i, j+1, 0) + cos_a*div_coeffs(i,  j+1, 1);

            nzCoeffs_div_x[count]  =  hVar[f] * (cos_a*div_coeffs(i, j+1, 0) - sin_a*div_coeffs(i,  j+1, 1));
            nzCoeffs_div_y[count]  =  hVar[f] * (sin_a*div_coeffs(i, j+1, 0) + cos_a*div_coeffs(i,  j+1, 1));

            count++;
        }


    }

    // assign the row indexes for the sparse matrix in CSR format

    rowIndx[0] = 0;         // first element is always zero

    // number of non-zero elements for pentagons
    for (i=0; i<12; i++) rowIndx[i+1] = 6*(i+1);

    // number of non-zero elements for hexagons
    for (i=12; i<node_num; i++) rowIndx[i+1] = (12*6) + 7*(i-11);

    // last element is always the number of non-zero coefficients
    rowIndx[node_num] = nNonZero;

    for (i=0; i<node_num; i++)
    {
      rowStart[i] = rowIndx[i];
      rowEnd[i] = rowIndx[i+1];
    }

    //----------------- Set up sparse matrix for div and grad ------------------


    // define handles to the sparse matrices for each component of the grad and
    // div operators. The Laplacian matrix will be constructed from these.
    sparse_matrix_t * DIV_X, * DIV_Y, * GRAD_X, * GRAD_Y;
    sparse_matrix_t * LAP_X, * LAP_Y;

    DIV_X   = new sparse_matrix_t;
    DIV_Y   = new sparse_matrix_t;
    GRAD_X  = new sparse_matrix_t;
    GRAD_Y  = new sparse_matrix_t;
    LAP_X   = new sparse_matrix_t;
    LAP_Y   = new sparse_matrix_t;

    // Assing memory to the handle for the Laplacian operator matrix
    operatorLaplacian = new sparse_matrix_t;  // defined in mesh.h
    // operatorLaplacian2 = new sparse_matrix_t;

    // Define c (0) based indexing
    sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;

    // Define all operations the original (non-transpose) matrix
    sparse_operation_t operation = SPARSE_OPERATION_NON_TRANSPOSE;

    // Assign non-zero values to matrices for the components of div and grad operators
    error = mkl_sparse_d_create_csr(DIV_X,  indexing, node_num, node_num, rowStart, rowEnd, colIndx, nzCoeffs_div_x);
    error = mkl_sparse_d_create_csr(DIV_Y,  indexing, node_num, node_num, rowStart, rowEnd, colIndx, nzCoeffs_div_y);

    error = mkl_sparse_d_create_csr(GRAD_X, indexing, node_num, node_num, rowStart, rowEnd, colIndx, nzCoeffs_grad_x);
    error = mkl_sparse_d_create_csr(GRAD_Y, indexing, node_num, node_num, rowStart, rowEnd, colIndx, nzCoeffs_grad_y);

    // Compute the multiplication of each component of the div dot grad
    error = mkl_sparse_spmm(operation, *DIV_X, *GRAD_X, LAP_X); // div_x * grad_x
    error = mkl_sparse_spmm(operation, *DIV_Y, *GRAD_Y, LAP_Y); // div_y * grad_y

    // Sum each component to give the total Grad^2 operator (Laplacian)
    // This is equivalent to (div_x * grad_x) + (div_y + grad_y)
    error = mkl_sparse_d_add(operation, *LAP_X, 1.0, *LAP_Y, operatorLaplacian);

    // error = mkl_sparse_d_create_csr(operatorLaplacian, indexing, node_num, node_num, rowStart, rowEnd, colIndx, nzCoeffs);

    // error = mkl_sparse_d_create_csr(operatorLaplacian,  indexing, node_num, node_num, rowStart, rowEnd, colIndx, nzCoeffs_div_x);
    // error = mkl_sparse_d_create_csr(operatorLaplacian2,  indexing, node_num, node_num, rowStart, rowEnd, colIndx, nzCoeffs_div_y);


    // Remove uneeded matrices from memory
    delete DIV_X, DIV_Y, GRAD_X, GRAD_Y, LAP_Y, LAP_X;


    // Create objects for the row and column indices, and non-zero coefficients
    // for the Laplace operator. Note that, although the Laplace operator was
    // just create and assigned above, we do not yet know its structure or the
    // number and values of the coefficients in its matrix.

    int    * rowStartLap;
    int      * rowEndLap;
    int     * colIndxLap;
    double * nzCoeffsLap;

    int nc, nr; //nrows, ncols
    int nNonZeroLap;

    rowStartLap = NULL;
    rowEndLap   = NULL;
    colIndxLap  = NULL;
    nzCoeffsLap = NULL;

    // Now extract the the row and column index arrays, and nonzero coefficients
    // from the Laplacian operator. These will be used to set up the solver for
    // the problem L*p = d below, where L is the Laplacian operator and is to be
    // inverted to find the solution p.
    error = mkl_sparse_d_export_csr (*operatorLaplacian, &indexing, &nc, &nr, &rowStartLap, &rowEndLap, &colIndxLap, &nzCoeffsLap);

    nNonZeroLap = rowEndLap[nr - 1];

    for (i=0; i<node_num; i++)
    {
      rowIndx[i] = rowStartLap[i];
    }
    rowIndx[node_num] = nNonZeroLap;



    // std::cout<<rowStartLap[node_num]<<std::endl;

    //----------------- Set up the intel MKL solver for L*p = d ----------------

    // create the solver object and assign it to a handle
    MKL_INT opt = MKL_DSS_ZERO_BASED_INDEXING;
    error = dss_create(pressureSolverHandle, opt);

    if (error != MKL_DSS_SUCCESS)
    {
        outstring << "ERROR: MKL couldn't create the solver object and exited with error code "<<error<< std::endl;
        globals->Output->Write(ERR_MESSAGE, &outstring);
        globals->Output->TerminateODIS();
    }

    // define the structure of the sparse matrix of size N*N with nNonZero elements
    opt =  MKL_DSS_SYMMETRIC_STRUCTURE;//MKL_DSS_SYMMETRIC_STRUCTURE;
    // error = dss_define_structure(pressureSolverHandle, opt, rowIndx, node_num, node_num, colIndx, nNonZero);
    error = dss_define_structure(pressureSolverHandle, opt, rowIndx, node_num, node_num, colIndxLap, nNonZeroLap);

    if (error != MKL_DSS_SUCCESS)
    {
        outstring << "ERROR: MKL couldn't create the laplacian matrix structure and exited with error code "<<error<< std::endl;
        globals->Output->Write(ERR_MESSAGE, &outstring);
        globals->Output->TerminateODIS();
    }

    // reorder the matrix to optimize the permutation order
    opt = MKL_DSS_AUTO_ORDER;
    error = dss_reorder(pressureSolverHandle, opt, 0);

    if (error != MKL_DSS_SUCCESS)
    {
        outstring << "ERROR: MKL couldn't reorder the laplacian matrix and exited with error code "<<error<< std::endl;
        globals->Output->Write(ERR_MESSAGE, &outstring);
        globals->Output->TerminateODIS();
    }

    // LU factorize the matrix using the nzCoeffs for direct solution
    opt = MKL_DSS_POSITIVE_DEFINITE;
    // error = dss_factor_real(pressureSolverHandle, opt, nzCoeffs);
    error = dss_factor_real(pressureSolverHandle, opt, nzCoeffsLap);

    if (error != MKL_DSS_SUCCESS)
    {
        outstring << "ERROR: MKL couldn't factorize the laplacian matrix and exited with error code "<<error<< std::endl;
        globals->Output->Write(ERR_MESSAGE, &outstring);
        globals->Output->TerminateODIS();
    }

    // delete pressure_matrix;

    delete nzCoeffs_grad_x;
    delete nzCoeffs_grad_y;
    delete nzCoeffs_div_x;
    delete nzCoeffs_div_y;

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

    // file_str = globals->path + SEP + "input_files" + SEP + "grid_l" + std::to_string(globals->geodesic_l.Value()) + "_order.txt";
    file_str = globals->path + SEP + "input_files" + SEP + "grid_l" + std::to_string(globals->geodesic_l.Value()) + ".txt";

    // in stream for input.in file
    std::ifstream gridFile(file_str, std::ifstream::in);

    if (gridFile.is_open())
    {
        outstring << std::endl << "Found mesh file: " + file_str << std::endl;
        globals->Output->Write(OUT_MESSAGE, &outstring);

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
        globals->Output->Write(ERR_MESSAGE, &outstring);
        globals->Output->TerminateODIS();
    }

    return 1;
};

// int Mesh::ReadLatLonFile(void)
// {
//     int i, j, k, N_ll;
//
//     N_ll = (int)globals->dLat.Value();
//
//     const H5std_string DSET_CellID("cell_ID");
//     const H5std_string DSET_VInv("vandermonde_inv");
//     const H5std_string DSET_Rot("rotation");
//
//     std::string file_str;
//     file_str = globals->path + SEP
//                + "input_files" + SEP
//                + "grid_l" + std::to_string(globals->geodesic_l.Value())
//                + '_' + std::to_string(N_ll)
//                + 'x' + std::to_string(N_ll)
//                + "_test.h5";
//
//     // Define file name and open file
//     H5std_string FILE_NAME(file_str);
//     H5File file(FILE_NAME, H5F_ACC_RDONLY);
//
//     // Access dataspaces in latlon file
//     DataSet dset_cellID = file.openDataSet(DSET_CellID);
//     DataSet dset_vInv = file.openDataSet(DSET_VInv);
//     DataSet dset_rot = file.openDataSet(DSET_Rot);
//
//     // Create filespaces for the correct rank and dimensions
//     DataSpace fspace_cellID = dset_cellID.getSpace();
//     DataSpace fspace_vInv = dset_vInv.getSpace();
//     DataSpace fspace_rot = dset_rot.getSpace();
//
//     // Get number of dimensions in the files dataspace
//     int rank_cellID = fspace_cellID.getSimpleExtentNdims();
//     int rank_vInv = fspace_vInv.getSimpleExtentNdims();
//     int rank_rot = fspace_rot.getSimpleExtentNdims();
//
//     hsize_t dims_cellID[2];
//     hsize_t dims_vInv[3];
//     hsize_t dims_rot[1];
//
//     // Get size of each dimension
//     rank_cellID = fspace_cellID.getSimpleExtentDims( dims_cellID );
//     rank_vInv = fspace_vInv.getSimpleExtentDims( dims_vInv );
//     rank_rot = fspace_rot.getSimpleExtentDims( dims_rot );
//
//     // Create memoryspace to read the datasets
//     DataSpace mspace_cellID(2, dims_cellID);
//     DataSpace mspace_vInv(3, dims_vInv);
//     DataSpace mspace_rot(1, dims_rot);
//
//     // Create 1D arrays to store file data
//     int * cellID_1D;
//     double * vInv_1D;
//     double * rot_1D;
//
//     cellID_1D = new int[180/N_ll * 360/N_ll];
//     vInv_1D = new double[node_num * 6 * 6];
//     rot_1D = new double[node_num];
//
//     // Read in the data
//     dset_cellID.read( cellID_1D, PredType::NATIVE_INT, mspace_cellID, fspace_cellID );
//     dset_vInv.read( vInv_1D, PredType::NATIVE_DOUBLE, mspace_vInv, fspace_vInv );
//
//     // Load Array classes with 1D dynamic arrays
//     int count = 0;
//     for (i=0; i<node_num; i++)
//     {
//         for (j=0; j<6; j++)
//         {
//             for (k=0; k<6; k++)
//             {
//                 V_inv(i, j, k) = vInv_1D[count];
//                 count++;
//             }
//         }
//     }
//
//     count = 0;
//
//     int ID;
//     double lat1, lat2, lon1, lon2;
//     double *m, *x, *y;
//     double r;
//
//     m = new double;
//
//     double test_solution_gg[node_num];
//     double test_solution_ll[180/N_ll][360/N_ll];
//     r = 1.0;//globals->radius.Value();
//
//     for (i=0; i<180/N_ll; i++)
//     {
//         for (j=0; j<360/N_ll; j++)
//         {
//             cell_ID(i, j) = cellID_1D[count];
//             count++;
//
//             // get cell ID which contains current lat-lon grid point
//             ID = cell_ID(i, j);
//
//             // get sph coords of cell with current lat-lon node
//             lat1 = node_pos_sph(ID,0);
//             lon1 = node_pos_sph(ID,1);
//
//             // calculate regular lat-lon position in radians
//             lat2 = (90.0 - (double)(i*N_ll))*radConv;
//             lon2 = (double)(j*N_ll)*radConv;
//
//             // std::cout<<lat2/radConv<<' '<<lat1/radConv<<std::endl;
//             // std::cout<<lat1/radConv<<' '<<lon2/radConv<<' '<<lon1/radConv<<std::endl;
//
//             test_solution_ll[i][j] = cos(3. * lat2) * sin(5 * lon2);
//
//             // std::cout<<lon2<<' '<<lat2<<' '<<test_solution_ll[i][j]<<std::endl;
//
//             // set pointers to mapped coords of lat-lon node
//             x = &ll_map_coords(i, j, 0);
//             y = &ll_map_coords(i, j, 1);
//             // *m = 0.0;
//
//             // call mapping function to find x y of current lat-lon node
//             mapAtPoint(*m, *x, *y, lat1, lat2, lon1, lon2, r);
//         }
//     }
//
//     delete m;
//
//
//     delete[] cellID_1D;
//     delete[] vInv_1D;
//     delete[] rot_1D;
//
// };

int Mesh::ReadWeightingFile(void)
{
    int i, N_ll;

    N_ll = (int)globals->dLat.Value();

    const H5std_string DSET_cols("column index");
    const H5std_string DSET_rows("row index");
    const H5std_string DSET_data("weights");

    std::string file_str;
    file_str = globals->path + SEP
               + "input_files" + SEP
               + "grid_l" + std::to_string(globals->geodesic_l.Value())
               + '_' + std::to_string(N_ll)
               + 'x' + std::to_string(N_ll)
               + "_weights.h5";

    // Define file name and open file
    H5std_string FILE_NAME(file_str);
    H5File file(FILE_NAME, H5F_ACC_RDONLY);

    // Access dataspaces in latlon file
    DataSet dset_cols = file.openDataSet(DSET_cols);
    DataSet dset_rows = file.openDataSet(DSET_rows);
    DataSet dset_data = file.openDataSet(DSET_data);

    // Create filespaces for the correct rank and dimensions
    DataSpace fspace_cols = dset_cols.getSpace();
    DataSpace fspace_rows = dset_rows.getSpace();
    DataSpace fspace_data = dset_data.getSpace();

    // Get number of dimensions in the files dataspace
    int rank_cols = fspace_cols.getSimpleExtentNdims();
    int rank_rows = fspace_rows.getSimpleExtentNdims();
    int rank_data = fspace_data.getSimpleExtentNdims();

    hsize_t dims_cols[1];    // length no of non-zero elements
    hsize_t dims_rows[1];
    hsize_t dims_data[1];    // length no of non-zero elements

    // Get size of each dimension
    rank_cols = fspace_cols.getSimpleExtentDims( dims_cols );
    rank_rows = fspace_rows.getSimpleExtentDims( dims_rows );
    rank_data = fspace_data.getSimpleExtentDims( dims_data );

    // Create memoryspace to read the datasets
    DataSpace mspace_cols(1, dims_cols);
    DataSpace mspace_rows(1, dims_rows);
    DataSpace mspace_data(1, dims_data);

    // Create 1D arrays to store file data
    int     * interpRows;
    int     * interpCols;
    double  * interpWeights;

    interpCols =    new    int[ dims_cols[0] ];
    interpRows =    new    int[ dims_rows[0] ];
    interpWeights = new double[ dims_data[0] ];

    // Read in the data
    dset_cols.read( interpCols, PredType::NATIVE_INT, mspace_cols, fspace_cols );
    dset_rows.read( interpRows, PredType::NATIVE_INT, mspace_rows, fspace_rows );
    dset_data.read( interpWeights, PredType::NATIVE_DOUBLE, mspace_data, fspace_data );

    mkl_set_num_threads(1);

    // Create  matrix handle using loaded data
    sparse_index_base_t index_type = SPARSE_INDEX_BASE_ZERO;     // we employ 0-based indexing.
    sparse_status_t err;

    int nrows = (360/N_ll)*(180/N_ll);
    int ncols = 3*node_num;

    interpMatrix = new sparse_matrix_t;
    err = mkl_sparse_d_create_csr(interpMatrix, index_type, nrows, ncols, interpRows, interpRows+1, interpCols, interpWeights);


    // Sparse interpolation matrix successfully created. Now we must provide
    // additional information to the matrix handle for optimization purposes

    // Expected number of calls...?
    // sparse_operation_t operation = SPARSE_OPERATION_NON_TRANSPOSE;
    matrix_descr descrp;
    descrp.type = SPARSE_MATRIX_TYPE_GENERAL;
    descrp.mode = SPARSE_FILL_MODE_LOWER;
    descrp.diag = SPARSE_DIAG_NON_UNIT;

    sparse_operation_t operation = SPARSE_OPERATION_NON_TRANSPOSE;

    int numberOfExpectedCalls = 100000000;  // Expect a big number!
    err = mkl_sparse_set_dotmv_hint(*interpMatrix, operation, descrp, numberOfExpectedCalls);

    if (err>0) {
        std::cout<<"Intel matrix optimization error!"<<std::endl;
    }

    // Now optimize matrix
    err = mkl_sparse_optimize(*interpMatrix);

    delete[] interpRows;
    delete[] interpCols;
    delete[] interpWeights;

    return 1;

};
