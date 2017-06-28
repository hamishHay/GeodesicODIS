#include "mesh.h"
#include "globals.h"
#include "array3d.h"
#include "array2d.h"
#include "array1d.h"

void pressureGradient(Mesh * mesh, Array2D<double> & dvdt, Array1D<double> & pressure)
{
    int node_num, friend_num;
    int i, j, j1, j2;

    Array2D<int> * friend_list;
    Array2D<double> * cent_map_factor;
    Array3D<double> * element_areas;
    Array3D<double> * normal_vecs;
    Array2D<double> * edge_lens;
    Array1D<double> * cv_areas;

    double m;                          // mapping factor at current cv edge
    double a0, a1, a2;
    double p0, p1, p2, p0_cent, p1_cent, p_avg;
    double x_grad, y_grad;
    double nx, ny;
    double edge_len;
    double g;

    node_num = mesh->node_num;
    friend_list = &(mesh->node_friends);
    cent_map_factor = &(mesh->control_vol_edge_centre_m);
    element_areas = &(mesh->node_friend_element_areas_map);
    normal_vecs = &(mesh->control_vol_edge_normal_map);
    edge_lens = &(mesh->control_vol_edge_len);
    cv_areas = &(mesh->control_volume_surf_area_map);
    g = mesh->globals->g.Value();

    for (i=0; i<node_num; i++)
    {
        friend_num = 6;
        if ((*friend_list)(i, 5) == -1) friend_num = 5;

        p0 = pressure(i);

        x_grad = 0.0;
        y_grad = 0.0;

        for (j=0; j<friend_num; j++)
        {
            // find avg pressure in element j
            j1 = j%friend_num;
            j2 = (j+1)%friend_num;

            a0 = (*element_areas)(i,j1,0);
            a1 = (*element_areas)(i,j1,1);
            a2 = (*element_areas)(i,j1,2);

            p1 = pressure((*friend_list)(i,j1));
            p2 = pressure((*friend_list)(i,j2));

            p0_cent = (p0 * a1 + p1 * a2 + p2 * a0) / (a0 + a1 + a2);

            // find avg pressure in element j+1
            j1 = j2;
            j2 = (j+2)%friend_num;

            a0 = (*element_areas)(i,j1,0);
            a1 = (*element_areas)(i,j1,1);
            a2 = (*element_areas)(i,j1,2);

            p1 = pressure((*friend_list)(i,j1));
            p2 = pressure((*friend_list)(i,j2));

            p1_cent = (p0 * a1 + p1 * a2 + p2 * a0) / (a0 + a1 + a2);

            // Find average p at the center of the control volume edge
            p_avg = 0.5*(p0_cent + p1_cent);


            j1 = j%friend_num;

            // get mapping factor for edge i,j2
            m = (*cent_map_factor)(i, j1);

            // get components of the edge normal vector
            nx = (*normal_vecs)(i, j1, 0);
            ny = (*normal_vecs)(i, j1, 1);

            // get edge length of current edge
            edge_len = (*edge_lens)(i,j1);

            // calculate x gradient
            x_grad += m * p_avg * nx * edge_len;

            // calculate y gradient
            y_grad += m * p_avg * ny * edge_len;


        }
        x_grad = g * x_grad / (*cv_areas)(i);
        y_grad = g * y_grad / (*cv_areas)(i);

        // if (i== 1280) std::cout<<x_grad<<'\t'<<y_grad<<std::endl;
        dvdt(i,0) -= x_grad;
        dvdt(i,1) -= y_grad;

    }
};

void velocityDivergence(Mesh * mesh, Array1D<double> & dpdt, Array2D<double> & velocity)
{
    int node_num, friend_num;
    int i, j, j1, j2, i1, i2;

    Array2D<int> * friend_list;
    Array2D<double> * cent_map_factor;
    Array3D<double> * element_areas;
    Array3D<double> * normal_vecs;
    Array2D<double> * edge_lens;
    Array1D<double> * cv_areas;
    Array3D<double> * vel_transform;

    double m;                          // mapping factor at current cv edge
    double a0, a1, a2;
    double u0, u1, u2, u0_cent, u1_cent, u_avg, u_temp;
    double v0, v1, v2, v0_cent, v1_cent, v_avg, v_temp;
    double div;
    double nx, ny;
    double edge_len;
    double h;
    double cos_a, sin_a;

    node_num = mesh->node_num;
    friend_list = &(mesh->node_friends);
    cent_map_factor = &(mesh->control_vol_edge_centre_m);
    element_areas = &(mesh->node_friend_element_areas_map);
    normal_vecs = &(mesh->control_vol_edge_normal_map);
    edge_lens = &(mesh->control_vol_edge_len);
    cv_areas = &(mesh->control_volume_surf_area_map);
    h = mesh->globals->h.Value();
    vel_transform = &(mesh->node_vel_trans);

    for (i=0; i<node_num; i++)
    {
        friend_num = 6;
        if ((*friend_list)(i, 5) == -1) friend_num = 5;

        u0 = velocity(i,0);
        v0 = velocity(i,1);

        // cos_a = (*vel_transform)(i,0,0);
        // sin_a = (*vel_transform)(i,0,1);
        //
        // u0 = u_temp * cos_a + v_temp * sin_a;
        // v0 = -u_temp * sin_a + v_temp * cos_a;

        div = 0.0;

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

            // FIND VELOCITY AT SECOND FREIND

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
            // u_avg = 0.5*(u0_cent + u1_cent);
            v_avg = 0.5*(v0_cent + v1_cent);
            u_avg = 0.5*(u0_cent + u1_cent);


            j1 = j%friend_num;

            // get mapping factor for edge i,j2
            m = (*cent_map_factor)(i, j1);

            // get components of the edge normal vector
            nx = (*normal_vecs)(i, j1, 0);
            ny = (*normal_vecs)(i, j1, 1);

            // get edge length of current edge
            edge_len = (*edge_lens)(i,j1);

            // calculate control volume divergence
            div += ((u_avg * nx) + (v_avg * ny)) * edge_len / m;

        }
        div /= (*cv_areas)(i);

        dpdt(i) = -h*div;

    }
};



// Function to calculate the Laplacian and diffusion for the velocity field.
// The discretized operator used here is the second order accurate finite
// difference operator given as equation 3.2 in Heikes et al (2013).
void velocityDiffusion(Mesh * mesh, Array2D<double> & dvdt, Array2D<double> & velocity)
{
    int node_num, friend_num;
    int i, j, i1, j1;

    Array2D<int> * friend_list;
    Array2D<double> * cent_map_factor;
    Array2D<double> * edge_lens;
    Array1D<double> * cv_areas;
    Array3D<double> * vel_transform;
    Array2D<double> * node_friend_dists;

    double m;                          // mapping factor at current cv edge
    double u0, u1, u_temp;
    double v0, v1, v_temp;
    double edge_len;
    double cos_a, sin_a;
    double lap_u, lap_v;
    double node_friend_dist;

    node_num = mesh->node_num;
    friend_list = &(mesh->node_friends);
    cent_map_factor = &(mesh->control_vol_edge_centre_m);
    edge_lens = &(mesh->control_vol_edge_len);
    cv_areas = &(mesh->control_volume_surf_area_map);
    vel_transform = &(mesh->node_vel_trans);
    node_friend_dists = &(mesh->node_dists);

    for (i=0; i<node_num; i++)
    {
        friend_num = 6;
        if ((*friend_list)(i, 5) == -1) friend_num = 5;

        u0 = velocity(i,0);
        v0 = velocity(i,1);

        lap_u = 0.0;
        lap_v = 0.0;

        for (j=0; j<friend_num; j++)
        {
            j1 = j%friend_num;
            i1 = (*friend_list)(i,j1);

            u_temp = velocity(i1,0);
            v_temp = velocity(i1,1);

            cos_a = (*vel_transform)(i, j1+1, 0);
            sin_a = (*vel_transform)(i, j1+1, 1);

            u1 = u_temp * cos_a + v_temp * sin_a;
            v1 = -u_temp * sin_a + v_temp * cos_a;

            node_friend_dist = (*node_friend_dists)(i,(j1+1)%friend_num);
            edge_len = (*edge_lens)(i,j1);

            lap_u += (u1 - u0)/node_friend_dist * edge_len;
            lap_v += (v1 - v0)/node_friend_dist * edge_len;

          }

          dvdt(i,0) += 2.0e3 * lap_u/(*cv_areas)(i);
          dvdt(i,1) += 2.0e3 * lap_v/(*cv_areas)(i);
    }

};
