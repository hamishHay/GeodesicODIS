#include "mesh.h"
#include "globals.h"
#include "array3d.h"
#include "array2d.h"
#include "array1d.h"

void coriolisForce(Mesh * mesh, Array2D<double> & dvdt, Array2D<double> &velocity)
{
    int node_num, friend_num;
    int i, j, i1, i2, j1, j2;

    Array2D<int> * friend_list;
    Array2D<double> * cent_dists;
    Array2D<double> * cent_map_factor;
    Array3D<double> * vel_transform;
    Array3D<double> * element_areas;


    double omega;
    double * sinLat;
    double dist_sum;
    double dist;
    double u_avg, v_avg;
    double u_sph, v_sph;
    double cos_a, sin_a;
    double u0, u1, u2, u_centre;
    double v0, v1, v2, v_centre;
    double a0, a1, a2;
    double m;

    node_num = mesh->node_num;
    friend_list = &(mesh->node_friends);
    cent_map_factor = &(mesh->control_vol_edge_centre_m);
    cent_dists = &(mesh->centroid_node_dists_map);
    vel_transform = &(mesh->node_vel_trans);
    element_areas = &(mesh->node_friend_element_areas_map);


    omega = mesh->globals->angVel.Value();

    sinLat = &(mesh->trigLat(0,1));

    for (i=0; i<node_num; i++)
    {
        dvdt(i,0) += 2.0 * omega * sinLat[i*2] * velocity(i,1);
        dvdt(i,1) -= 2.0 * omega * sinLat[i*2] * velocity(i,0);
    }

};
