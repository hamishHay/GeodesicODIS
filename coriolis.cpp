#include "mesh.h"
#include "globals.h"
#include "array3d.h"
#include "array2d.h"
#include "array1d.h"

void coriolisForce(Mesh * mesh, Array2D<double> & dvdt, Array2D<double> &velocity)
{
    int node_num, friend_num;
    int i, j, i1, i2, j1, j2;
    int method = 1;

    Array2D<int> * friend_list;
    Array2D<double> * cent_map_factor;
    Array3D<double> * vel_transform;
    Array3D<double> * element_areas;


    double omega;
    double * sinLat;
    double u_avg, v_avg;
    double cos_a, sin_a;
    double u0, u1, u2, u0_cent, u1_cent;
    double v0, v1, v2, v0_cent, v1_cent;
    double u_temp, v_temp;
    double a0, a1, a2;
    double m;
    double total_avg_u, total_avg_v;

    node_num = mesh->node_num;
    friend_list = &(mesh->node_friends);
    cent_map_factor = &(mesh->control_vol_edge_centre_m);
    vel_transform = &(mesh->node_vel_trans);
    element_areas = &(mesh->node_friend_element_areas_map);


    omega = mesh->globals->angVel.Value();

    sinLat = &(mesh->trigLat(0,1));

    switch (method)
    {
    case 0:
        for (i=0; i<node_num; i++)
        {
            dvdt(i,0) += 2.0 * omega * sinLat[i*2] * velocity(i,1);
            dvdt(i,1) -= 2.0 * omega * sinLat[i*2] * velocity(i,0);
        }
        break;

    case 1:
        for (i=0; i<node_num; i++)
        {
            friend_num = 6;
            if ((*friend_list)(i, 5) == -1) friend_num = 5;

            u0 = velocity(i,0);
            v0 = velocity(i,1);

            total_avg_v = 0.0;
            total_avg_u = 0.0;
            // dmass = 0.0;

            for (j=0; j<friend_num; j++)
            {
                // find avg pressure in element j
                j1 = j%friend_num;
                j2 = (j+1)%friend_num;
                i1 = (*friend_list)(i,j1);
                i2 = (*friend_list)(i,j2);


                a0 = (*element_areas)(i,j1,0);
                a1 = (*element_areas)(i,j1,1);
                a2 = (*element_areas)(i,j1,2);

                // FIND VELOCITY AT FIRST FRIEND
                u_temp = velocity(i1,0);
                v_temp = velocity(i1,1);

                // CONVERT TO MAPPED VELOCITIES

                cos_a = (*vel_transform)(i, j1+1, 0);
                sin_a = (*vel_transform)(i, j1+1, 1);
                u1 = u_temp * cos_a + v_temp * sin_a;
                v1 = -u_temp * sin_a + v_temp * cos_a;

                // FIND VELOCITY AT SECOND FRIEND
                u_temp = velocity(i2,0);
                v_temp = velocity(i2,1);

                // CONVERT TO MAPPED VELOCITIES

                cos_a = (*vel_transform)(i, j2+1, 0);
                sin_a = (*vel_transform)(i, j2+1, 1);
                u2 = u_temp * cos_a + v_temp * sin_a;
                v2 = -u_temp * sin_a + v_temp * cos_a;

                // FIND AVERAGE VELOCITY AT FIRST ELEMENT CENTRE
                u0_cent = (u0 * a1 + u1 * a2 + u2 * a0) / (a0 + a1 + a2);
                v0_cent = (v0 * a1 + v1 * a2 + v2 * a0) / (a0 + a1 + a2);

                // find avg pressure in element j+1
                j1 = j2;
                j2 = (j+2)%friend_num;
                i1 = (*friend_list)(i,j1);
                i2 = (*friend_list)(i,j2);

                a0 = (*element_areas)(i,j1,0);
                a1 = (*element_areas)(i,j1,1);
                a2 = (*element_areas)(i,j1,2);

                // FIND VELOCITY AT FIRST FRIEND
                u_temp = velocity(i1,0);
                v_temp = velocity(i1,1);

                // CONVERT TO MAPPED VELOCITIES
                cos_a = (*vel_transform)(i, j1+1, 0);
                sin_a = (*vel_transform)(i, j1+1, 1);
                u1 = u_temp * cos_a + v_temp * sin_a;
                v1 = -u_temp * sin_a + v_temp * cos_a;

                // FIND VELOCITY AT SECOND FREIND
                u_temp = velocity(i2,0);
                v_temp = velocity(i2,1);

                // CONVERT TO MAPPED VELOCITIES
                cos_a = (*vel_transform)(i, j2+1, 0);
                sin_a = (*vel_transform)(i, j2+1, 1);
                u2 = u_temp * cos_a + v_temp * sin_a;
                v2 = -u_temp * sin_a + v_temp * cos_a;

                // FIND AVERAGE VELOCITY AT FIRST ELEMENT CENTRE
                u1_cent = (u0 * a1 + u1 * a2 + u2 * a0) / (a0 + a1 + a2);
                v1_cent = (v0 * a1 + v1 * a2 + v2 * a0) / (a0 + a1 + a2);

                // Find average p at the center of the control volume edge
                v_avg = 0.5*(v0_cent + v1_cent);
                u_avg = 0.5*(u0_cent + u1_cent);

                j1 = j%friend_num;

                // get mapping factor for edge i,j2
                m = (*cent_map_factor)(i, j1);

                // calculate control volume divergence
                total_avg_u += u_avg*m;
                total_avg_v += v_avg*m;

            }

            total_avg_u /= friend_num;
            total_avg_v /= friend_num;

            dvdt(i,0) += 2.0 * omega * sinLat[i*2] * total_avg_v;
            dvdt(i,1) -= 2.0 * omega * sinLat[i*2] * total_avg_u;

        }
        break;
    }



};
