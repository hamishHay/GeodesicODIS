#include "mesh.h"
#include "globals.h"
#include "array3d.h"
#include "array2d.h"
#include "array1d.h"
#include "sphericalHarmonics.h"
#include "interpolation.h"

// static double dummy_sum = 0;

// #pragma omp for
void pressureGradient(Mesh * mesh, Array2D<double> & dvdt, Array1D<double> & pressure, double g = 1.0)
{
    int node_num, friend_num;
    int i, j, j1, j2;
    int end_i;

    Array2D<int> * friend_list;
    Array2D<double> * cent_map_factor;
    Array3D<double> * element_areas;
    Array3D<double> * normal_vecs;
    Array2D<double> * edge_lens;
    Array1D<double> * cv_areas;
    Array1D<int> * mask;

    double m;                          // mapping factor at current cv edge
    double a0, a1, a2;
    double p0, p1, p2, p0_cent, p1_cent, p_avg;
    double x_grad, y_grad;
    double nx, ny;
    double edge_len;
    // double g;

    node_num = mesh->node_num;
    friend_list = &(mesh->node_friends);
    cent_map_factor = &(mesh->control_vol_edge_centre_m);
    element_areas = &(mesh->node_friend_element_areas_map);
    normal_vecs = &(mesh->control_vol_edge_normal_map);
    edge_lens = &(mesh->control_vol_edge_len);
    cv_areas = &(mesh->control_volume_surf_area_map);
    // g = mesh->globals->g.Value();
    mask = &(mesh->land_mask);

    end_i = node_num;
    if (mesh->globals->surface_type == FREE_LOADING ||
        mesh->globals->surface_type == LID_MEMBR ||
        mesh->globals->surface_type == LID_INF)
    {
        end_i = 2;
    }


    for (i=0; i<end_i; i++)
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

            p1_cent = (p0*a1 + p1*a2 + p2*a0) / (a0 + a1 + a2);

            // Find average p at the center of the control volume edge
            p_avg = 0.5*(p0_cent + p1_cent);


            j1 = j%friend_num;

            // get mapping factor for edge i,j2
            m = (*cent_map_factor)(i, j1);

            // get components of the edge normal vector
            nx = (*normal_vecs)(i, j1, 0);
            ny = (*normal_vecs)(i, j1, 1);

            // get edge length of current edge
            edge_len = (*edge_lens)(i, j1);

            // calculate x gradient
            x_grad += m * p_avg * nx * edge_len;

            // calculate y gradient
            y_grad += m * p_avg * ny * edge_len;

        }

        x_grad = g * x_grad / (*cv_areas)(i);
        y_grad = g * y_grad / (*cv_areas)(i);

        dvdt(i,0) -= x_grad;
        dvdt(i,1) -= y_grad;

        // dvdt(i,0) -= y_grad;
        // dvdt(i,1) -= x_grad;
    }
};

void pressureGradientSH(Globals * globals, Mesh * mesh, Array2D<double> & dvdt, Array1D<double> & gg_scalar, double factor)
{
    int node_num, l_max, N_ll;
    int i, j, f, l, m;
    int start_l;
    double lat_factor, lon_factor;
    double r;
    double dvdt_x, dvdt_y;
    double dvdt_x_total, dvdt_y_total;
    double dummy_val, dummy_val_total;
    double * loading_factor;
    double * shell_factor_beta;
    double * cosLat;

    Array3D<double> * scalar_lm;
    Array3D<double> * scalar_lm_dummy;

    Array1D<double> * scalar_dummy;

    Array3D<double> * Pbar_lm;
    Array3D<double> * Pbar_lm_deriv;
    Array3D<double> * trigMLon;
    Array2D<double> * trigLat;

    Array2D<double> * ll_scalar;

    Array2D<int> * friend_list;
    friend_list = &(mesh->node_friends);

    node_num = globals->node_num;
    l_max = globals->l_max.Value();
    N_ll = (int)globals->dLat.Value();
    r = globals->radius.Value();
    lat_factor = factor * 1.0/r;

    loading_factor = globals->loading_factor;
    shell_factor_beta = globals->shell_factor_beta;

    scalar_lm = new Array3D<double>(2.*(l_max+1), 2.*(l_max+1), 2);
    scalar_lm_dummy = new Array3D<double>(2.*(l_max+1), 2*(l_max+1), 2);
    ll_scalar = new Array2D<double>(180/N_ll, 360/N_ll);
    scalar_dummy = new Array1D<double>(node_num);

    Pbar_lm = &(mesh->Pbar_lm);             // 4-pi Normalised associated leg funcs
    Pbar_lm_deriv = &(mesh->Pbar_lm_deriv); // 4-pi Normalised associated leg funcs derivs
    trigMLon = &(mesh->trigMLon);           // cos and sin of m*longitude
    trigLat = &(mesh->trigLat);
    cosLat = &(mesh->trigLat(0,0));

    start_l = 2;
    if (globals->surface_type == LID_LOVE) start_l = 0;
    if (globals->surface_type == FREE_LOADING) start_l = 2;

    // INTERPOLATE GEODESIC GRID DATA TO LAT-LON GRID
    interpolateGG2LL(globals,
                       mesh,
                       *ll_scalar,
                       gg_scalar,
                       mesh->ll_map_coords,
                       mesh->V_inv,
                       mesh->cell_ID);

    // FIND SPHERICAL HARMONIC EXPANSION COEFFICIENTS
    // OF THE PRESSURE FIELD
    getSHCoeffsLL(*ll_scalar, *scalar_lm, N_ll, l_max);

    // getSHCoeffsGG((mesh->node_pos_sph), gg_scalar, *scalar_lm, node_num, l_max);

    //for (l=0; l<l_max+1; l++)
    //{
    //    for (m=0; m<=l; m++)
    //    {
        //std::cout<<l<<' '<<m<<' '<<(*scalar_lm)(l, m, 1)<<'\t'<<(*scalar_lm_dummy)(l, m, 1)<<std::endl;
    //    }
    //}



    int method = 1;
    switch (method)
    {
    case 1:
        // LOOP THROUGH EVERY GRID POINT, REBUILDING SCALAR
        // GRADIENT FROM SPHERICAL HARMONIC COEFFICIENTS
        for (i=0; i<node_num; i++)
        {
            // lon_factor = factor * 1.0/(r*(*trigLat)(i,0));
            lon_factor = factor * 1.0/(r*cosLat[i*2]);


            dvdt_x_total = 0.0;
            dvdt_y_total = 0.0;
            dummy_val_total = 0.0;

            (*scalar_dummy)(i) = 0.0;
            for (l=start_l; l<l_max+1; l++)
            {
                dvdt_x = 0.0;
                dvdt_y = 0.0;
                dummy_val = 0.0;
                for (m=0; m<=l; m++)
                {
                    dvdt_x += lon_factor * (*Pbar_lm)(i, l, m)
                              * (-(*scalar_lm)(l, m, 0) * (double)m * (*trigMLon)(i, m, 1)
                              + (*scalar_lm)(l, m, 1) * (double)m * (*trigMLon)(i, m, 0));

                    dvdt_y += lat_factor * (*Pbar_lm_deriv)(i, l, m)
                              * ((*scalar_lm)(l, m, 0) * (*trigMLon)(i, m, 0)
                              + (*scalar_lm)(l, m, 1) * (*trigMLon)(i, m, 1));


                    dummy_val += (*Pbar_lm)(i, l, m)
                              * ((*scalar_lm)(l, m, 0) * (*trigMLon)(i, m, 0)
                              + (*scalar_lm)(l, m, 1) * (*trigMLon)(i, m, 1));
                }

                if (globals->surface_type == FREE_LOADING)
                {
                    dvdt_x *= loading_factor[l];
                    dvdt_y *= loading_factor[l];
                    dummy_val *= loading_factor[l];
                }
                else if (globals->surface_type == LID_LOVE
                         || globals->surface_type == LID_MEMBR)
                {
                    dvdt_x *= shell_factor_beta[l];
                    dvdt_y *= shell_factor_beta[l];
                    dummy_val *= shell_factor_beta[l];
                }


                dvdt_x_total += dvdt_x;
                dvdt_y_total += dvdt_y;
                dummy_val_total += dummy_val;
            }

            (*scalar_dummy)(i) = dummy_val_total;

            if (i > 1) {
                dvdt(i,0) += dvdt_x_total;
                dvdt(i,1) += dvdt_y_total;
            }
        }

        pressureGradient(mesh, dvdt, *scalar_dummy, -factor);

        break;

    case 2:
        // LOOP THROUGH EVERY GRID POINT, REBUILDING SCALAR
        // GRADIENT FROM SPHERICAL HARMONIC COEFFICIENTS
        for (i=0; i<node_num; i++)
        {
            // lon_factor = factor * 1.0/(r*(*trigLat)(i,0));
            // std::cout<<gg_scalar(i)<<'\t';
            // gg_scalar(i) = 0.0;     // THIS IS CHEATING
            dvdt_x_total = 0.0;

            for (l=0.0; l<l_max+1; l++)
            {
                dvdt_x = 0.0;
                for (m=0; m<=l; m++)
                {
                    dvdt_x += (*Pbar_lm)(i, l, m)
                              * ((*scalar_lm)(l, m, 0) * (*trigMLon)(i, m, 0)
                              + (*scalar_lm)(l, m, 1) * (*trigMLon)(i, m, 1));
                }

                if (globals->surface_type == FREE_LOADING)
                {
                    dvdt_x *= loading_factor[l];
                }
                else if (globals->surface_type == LID_LOVE ||
                         globals->surface_type == LID_MEMBR)
                {
                    dvdt_x *= shell_factor_beta[l];
                }

                dvdt_x_total += dvdt_x;
            }

            // std::cout<<gg_scalar(i)<<std::endl;
            (*scalar_dummy)(i) += dvdt_x_total;

        }

        pressureGradient(mesh, dvdt, *scalar_dummy, -factor);

        break;

    }

    delete scalar_dummy;
    delete scalar_lm;
    delete ll_scalar;
    delete scalar_lm_dummy;

}

void velocityDivergence(Mesh * mesh, Array1D<double> & dpdt, Array2D<double> & velocity, double & sum, double h = 1.0)
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
    // Array1D<double> * mass;
    Array1D<int> * mask;

    double m;                          // mapping factor at current cv edge
    double a0, a1, a2;
    double u0, u1, u2, u0_cent, u1_cent, u_avg, u_temp;
    double v0, v1, v2, v0_cent, v1_cent, v_avg, v_temp;
    double div;
    double nx, ny;
    double edge_len;
    // double total_dmass, dmass, dt;
    // double h;
    double cos_a, sin_a;

    node_num = mesh->node_num;
    friend_list = &(mesh->node_friends);
    cent_map_factor = &(mesh->control_vol_edge_centre_m);
    element_areas = &(mesh->node_friend_element_areas_map);
    normal_vecs = &(mesh->control_vol_edge_normal_map);
    edge_lens = &(mesh->control_vol_edge_len);
    cv_areas = &(mesh->control_volume_surf_area_map);
    // mass = &(mesh->control_volume_mass);
    // dt = mesh->globals->timeStep.Value();
    vel_transform = &(mesh->node_vel_trans);
    mask = &(mesh->land_mask);

    // total_dmass = 0.0;
    for (i=0; i<node_num; i++)
    {
        friend_num = 6;
        if ((*friend_list)(i, 5) == -1) friend_num = 5;

        u0 = velocity(i,0);
        v0 = velocity(i,1);

        div = 0.0;
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

            // get components of the edge normal vector
            nx = (*normal_vecs)(i, j1, 0);
            ny = (*normal_vecs)(i, j1, 1);

            // get edge length of current edge
            edge_len = (*edge_lens)(i,j1);

            // calculate control volume divergence
            div += ((u_avg * nx) + (v_avg * ny)) * edge_len / m;

        }

        // dmass = dt * div * 1000.0 * h;


        div /= (*cv_areas)(i);

        dpdt(i) = -h*div;


        if (sum >= 0.0) {
            sum += fabs(div);
            // total_dmass += dmass/(*mass)(i);
        }
    }

    // std::cout<<std::scientific<<total_dmass<<std::endl;
};

// Function to calculate the Laplacian and diffusion for the velocity field.
// The discretized operator used here is the second order accurate finite
// difference operator given as equation 3.2 in Heikes et al (2013).
void velocityDiffusion(Mesh * mesh, Array2D<double> & dvdt, Array2D<double> & velocity, double viscosity)
{
    int node_num, friend_num;
    int i, j, i1, j1;

    Array2D<int> * friend_list;
    Array2D<double> * cent_map_factor;
    Array2D<double> * edge_lens;
    Array1D<double> * cv_areas;
    Array3D<double> * vel_transform;
    Array2D<double> * node_friend_dists;
    Array1D<int> * mask;

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
    mask = &(mesh->land_mask);

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

        dvdt(i,0) += viscosity * lap_u/(*cv_areas)(i);
        dvdt(i,1) += viscosity * lap_v/(*cv_areas)(i);

    }
};

// Function to calculate the Laplacian and diffusion for the velocity field.
// The discretized operator used here is the second order accurate finite
// difference operator given as equation 3.2 in Heikes et al (2013).
void scalarDiffusion(Mesh * mesh, Array1D<double> & d2s, Array1D<double> & s, double viscosity)
{
    int node_num, friend_num;
    int i, j, i1, j1;

    Array2D<int> * friend_list;
    Array2D<double> * cent_map_factor;
    Array2D<double> * edge_lens;
    Array1D<double> * cv_areas;
    Array3D<double> * vel_transform;
    Array2D<double> * node_friend_dists;
    Array1D<int> * mask;

    double m;                          // mapping factor at current cv edge
    double s0, s1, s_temp;
    double edge_len;
    double cos_a, sin_a;
    double lap_s, lap_v;
    double node_friend_dist;

    node_num = mesh->node_num;
    friend_list = &(mesh->node_friends);
    cent_map_factor = &(mesh->control_vol_edge_centre_m);
    edge_lens = &(mesh->control_vol_edge_len);
    cv_areas = &(mesh->control_volume_surf_area_map);
    vel_transform = &(mesh->node_vel_trans);
    node_friend_dists = &(mesh->node_dists);
    mask = &(mesh->land_mask);

    for (i=0; i<node_num; i++)
    {
        friend_num = 6;
        if ((*friend_list)(i, 5) == -1) friend_num = 5;

        s0 = s(i);

        lap_s = 0.0;

        for (j=0; j<friend_num; j++)
        {
              j1 = j%friend_num;
              i1 = (*friend_list)(i,j1);

              s1 = s(i1);

              node_friend_dist = (*node_friend_dists)(i,j1);
              edge_len = (*edge_lens)(i,(j1+1)%friend_num);

              lap_s += (s1 - s0)/node_friend_dist * edge_len;

        }

        d2s(i) = viscosity * lap_s/(*cv_areas)(i);

    }
};

void smoothingSH(Globals * globals, Mesh * mesh, Array1D<double> & gg_scalar)
{
    int node_num, l_max, N_ll;
    int i, j, f, l, m;
    double r;
    double interp_val;
    double * cosLat;

    Array3D<double> * scalar_lm;

    Array3D<double> * Pbar_lm;
    Array3D<double> * Pbar_lm_deriv;
    Array3D<double> * trigMLon;
    Array2D<double> * trigLat;

    Array2D<double> * ll_scalar;

    node_num = globals->node_num;
    l_max = globals->l_max.Value();
    N_ll = (int)globals->dLat.Value();
    r = globals->radius.Value();

    scalar_lm = new Array3D<double>(2.*(l_max+1), 2.*(l_max+1), 2);
    ll_scalar = new Array2D<double>(180/N_ll, 360/N_ll);
    // scalar_dummy = new Array1D<double>(node_num);

    Pbar_lm = &(mesh->Pbar_lm);             // 4-pi Normalised associated leg funcs
    Pbar_lm_deriv = &(mesh->Pbar_lm_deriv); // 4-pi Normalised associated leg funcs derivs
    trigMLon = &(mesh->trigMLon);           // cos and sin of m*longitude
    trigLat = &(mesh->trigLat);
    cosLat = &(mesh->trigLat(0,0));

    // INTERPOLATE GEODESIC GRID DATA TO LAT-LON GRID
    interpolateGG2LL(globals,
                        mesh,
                        *ll_scalar,
                        gg_scalar,
                        mesh->ll_map_coords,
                        mesh->V_inv,
                        mesh->cell_ID);

    // FIND SPHERICAL HARMONIC EXPANSION COEFFICIENTS
    // OF THE PRESSURE FIELD
    getSHCoeffsLL(*ll_scalar, *scalar_lm, N_ll, l_max);

    for (i=0; i<node_num; i++)
    {
        interp_val = 0.0;
        for (l=0; l<l_max+1; l++)
        {
            for (m=0; m<=l; m++)
            {
                interp_val += (*Pbar_lm)(i, l, m)
                          * ((*scalar_lm)(l, m, 0) * (*trigMLon)(i, m, 0)
                          + (*scalar_lm)(l, m, 1) * (*trigMLon)(i, m, 1));
            }
        }

        gg_scalar(i) = interp_val;
    }

    delete scalar_lm;
    delete ll_scalar;

};

void avgAtPoles(Mesh * mesh, Array2D<double> & velocity)
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
    // Array1D<double> * mass;
    Array1D<int> * mask;

    double m;                          // mapping factor at current cv edge
    double a0, a1, a2;
    double u0, u1, u2, u0_cent, u1_cent, u_avg, u_temp;
    double v0, v1, v2, v0_cent, v1_cent, v_avg, v_temp;
    double div;
    double nx, ny;
    double edge_len;
    // double total_dmass, dmass, dt;
    // double h;
    double cos_a, sin_a;

    node_num = mesh->node_num;
    friend_list = &(mesh->node_friends);
    cent_map_factor = &(mesh->control_vol_edge_centre_m);
    element_areas = &(mesh->node_friend_element_areas_map);
    normal_vecs = &(mesh->control_vol_edge_normal_map);
    edge_lens = &(mesh->control_vol_edge_len);
    cv_areas = &(mesh->control_volume_surf_area_map);
    // mass = &(mesh->control_volume_mass);
    // dt = mesh->globals->timeStep.Value();
    vel_transform = &(mesh->node_vel_trans);
    mask = &(mesh->land_mask);

    // total_dmass = 0.0;
    for (i=0; i<2; i++)
    {
        friend_num = 5;

        u_avg = 0.0;
        v_avg = 0.0;
        for (j=0; j<friend_num; j++)
        {
            // FIND VELOCITY AT FRIEND
            i1 = (*friend_list)(i,j);

            u_temp = velocity(i1,0);
            v_temp = velocity(i1,1);

            // CONVERT TO MAPPED VELOCITIES
            cos_a = (*vel_transform)(i, j+1, 0);
            sin_a = (*vel_transform)(i, j+1, 1);
            u1 = u_temp * cos_a + v_temp * sin_a;
            v1 = -u_temp * sin_a + v_temp * cos_a;

            // ADD TO AVERAGE
            u_avg += u1;
            v_avg += v1;
        }

        // WEIGHT BY NUMBER OF SURROUNDING POINTS
        u_avg /= 5.0;
        v_avg /= 5.0;

        velocity(i,0) = u_avg;
        velocity(i,1) = v_avg;

    }

};
