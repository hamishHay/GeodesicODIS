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
    //
    // for (i=0; i<node_num; i++)
    // {
    //
    //     friend_num = 6;
    //     if ((*friend_list)(i, 5) == -1) friend_num = 5;
    //
    //     u_sph = velocity(i,0);
    //     v_sph = velocity(i,1);
    //
    //     cos_a = (*vel_transform)(i,0,0);
    //     sin_a = (*vel_transform)(i,0,1);
    //
    //     u0 = u_sph * cos_a + v_sph * sin_a;
    //     v0 = -u_sph * sin_a + v_sph * cos_a;
    //
    //     dist_sum = 0.0;
    //     u_avg = 0.0;
    //     v_avg = 0.0;
    //     for (j=0; j<friend_num; j++)
    //     {
    //
    //         j1 = j%friend_num;
    //         j2 = (j+1)%friend_num;
    //         i1 = (*friend_list)(i,j1);
    //         i2 = (*friend_list)(i,j2);
    //
    //
    //         a0 = (*element_areas)(i,j1,0);
    //         a1 = (*element_areas)(i,j1,1);
    //         a2 = (*element_areas)(i,j1,2);
    //
    //         // FIND VELOCITY AT FIRST FRIEND
    //
    //         u_sph = velocity(i1,0);
    //         v_sph = velocity(i1,1);
    //
    //
    //         // CONVERT TO MAPPED VELOCITIES
    //
    //         cos_a = (*vel_transform)(i, j1+1, 0);
    //         sin_a = (*vel_transform)(i, j1+1, 1);
    //         u1 = u_sph * cos_a + v_sph * sin_a;
    //         v1 = -u_sph * sin_a + v_sph * cos_a;
    //
    //
    //         // FIND VELOCITY AT SECOND FREIND
    //
    //         u_sph = velocity(i2,0);
    //         v_sph = velocity(i2,1);
    //
    //
    //         // CONVERT TO MAPPED VELOCITIES
    //
    //         cos_a = (*vel_transform)(i, j2+1, 0);
    //         sin_a = (*vel_transform)(i, j2+1, 1);
    //         u2 = u_sph * cos_a + v_sph * sin_a;
    //         v2 = -u_sph * sin_a + v_sph * cos_a;
    //
    //         // FIND AVERAGE VELOCITY AT FIRST ELEMENT CENTRE
    //
    //         u_centre = (u0 * a1 + u1 * a2 + u2 * a0) / (a0 + a1 + a2);
    //         v_centre = (v0 * a1 + v1 * a2 + v2 * a0) / (a0 + a1 + a2);
    //
    //
    //         i1 = (*friend_list)(i,j);
    //
    //         m = (*cent_map_factor)(i, j1);
    //
    //         dist = 1.0/((*cent_dists)(i,j)*m);
    //         dist_sum += dist;
    //
    //         u_avg += u_centre * dist;
    //         v_avg += v_centre * dist;
    //
    //     }
    //
    //     // divide sum by total distance, and apply.
    //
    //     u_avg /= dist_sum;
    //     v_avg /= dist_sum;
    //
    //     dvdt(i,0) += 2.0 * omega * sinLat[i*2] * v_avg;
    //     dvdt(i,1) -= 2.0 * omega * sinLat[i*2] * u_avg;
    //
    //     // if (i==5120) dvdt(i,0) -= 2.0 * omega * sinLat[i*2] * v_avg;
    //
    // }


};
