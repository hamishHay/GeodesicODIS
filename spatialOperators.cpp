#include "mesh.h"
#include "globals.h"
#include "array3d.h"
#include "array2d.h"
#include "array1d.h"
#include "sphericalHarmonics.h"
#include "interpolation.h"

extern "C"{
// FORTRAN adds _ after all the function names
// and all variables are called by reference
double dgemv_( const char * TRANS,
               const int * m,
               const int * n,
               const double * alpha,
               const double * V,
               const int * ldv,
               const double * x,
               const int * incx,
               const double * beta,
               const double * y,
               const int * incy);
}

// static double dummy_sum = 0;

// #pragma omp for
void pressureGradient(Mesh * mesh, Array2D<double> & dvdt, Array1D<double> & pressure, int N, double g = 1.0)
{
    int node_num, friend_num;
    int i, j, end_i;

    Array2D<int> * friend_list;
    Array3D<double> * grad_coeffs;

    double p0, p1;
    double x_grad, y_grad;

    node_num = mesh->node_num;
    friend_list = &(mesh->node_friends);
    grad_coeffs = &(mesh->grad_coeffs);

    // end_i = node_num;
    // if (mesh->globals->surface_type == FREE_LOADING ||
    //     mesh->globals->surface_type == LID_MEMBR ||
    //     mesh->globals->surface_type == LID_INF)
    // {
    //     end_i = 2;
    // }

    // end_i = node_num;

    for (i=0; i<N; i++)
    {
        friend_num = 6;
        if ((*friend_list)(i, 5) == -1) friend_num = 5;

        p0 = pressure(i);

        x_grad = p0 * (*grad_coeffs)(i, 0, 0);
        y_grad = p0 * (*grad_coeffs)(i, 0, 1);

        for (j=0; j<friend_num; j++)
        {
            p1 = pressure( (*friend_list)(i, j) );

            x_grad += p1 * (*grad_coeffs)(i, j+1, 0);
            y_grad += p1 * (*grad_coeffs)(i, j+1, 1);
        }

        dvdt(i,0) -= g*x_grad;
        dvdt(i,1) -= g*y_grad;
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
    double scalar_lm_cos, scalar_lm_sin;
    double cosMLon, sinMLon;
    double * sh_coeffs;
    double * soln_vec;


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

    double * sh_matrix;
    sh_matrix = &(mesh->sh_matrix_fort(0));
    soln_vec = &((*scalar_dummy)(0));

    start_l = 2;
    if (globals->surface_type == LID_LOVE) start_l = 2;
    if (globals->surface_type == FREE_LOADING) start_l = 2;

    sh_coeffs = new double[(l_max+1)*(l_max+2) - 6];

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

    char trans = 'n';
    int M = node_num;
    int N = (l_max+1)*(l_max+2) - 6;
    double alpha = 1.0;
    int lda = M;
    int incx = 1;
    int incy = 1;
    double beta = 0.0;


    int count = 0;

    int method = 3;
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
                // for (m=0; ; m++)
                {
                    scalar_lm_cos = (*scalar_lm)(l, m, 0);
                    scalar_lm_sin = (*scalar_lm)(l, m, 1);

                    cosMLon = (*trigMLon)(i, m, 0);
                    sinMLon = (*trigMLon)(i, m, 1);

                    dvdt_x += lon_factor * (*Pbar_lm)(i, l, m)
                              * (-scalar_lm_cos * (double)m * sinMLon
                              + scalar_lm_sin * (double)m * cosMLon);

                    dvdt_y += lat_factor * (*Pbar_lm_deriv)(i, l, m)
                              * (scalar_lm_cos * cosMLon
                              + scalar_lm_sin * sinMLon);


                    dummy_val += (*Pbar_lm)(i, l, m)
                              * (scalar_lm_cos * cosMLon
                              + scalar_lm_sin * sinMLon);

                    // if (m==l) break;

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

        pressureGradient(mesh, dvdt, *scalar_dummy, 2, -factor);

        break;

    case 2:
        // LOOP THROUGH EVERY GRID POINT, REBUILDING SCALAR
        // GRADIENT FROM SPHERICAL HARMONIC COEFFICIENTS
        for (i=0; i<node_num; i++)
        {
            dvdt_x_total = 0.0;

            for (l=0; l<l_max+1; l++)
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

            (*scalar_dummy)(i) += dvdt_x_total;

        }

        pressureGradient(mesh, dvdt, *scalar_dummy, node_num, -factor);

        break;


    case 3:
        // Using blas to perform a large matrix-vector calculation

        int count = 0;
        for (l=2; l<l_max+1; l++)
        {
            for (m=0; m<=l; m++)
            {
                sh_coeffs[count] = (*scalar_lm)(l, m, 0);
                count++;

                sh_coeffs[count] = (*scalar_lm)(l, m, 1);
                count++;
            }
        }

        // Perform sh_matrix * sh_coeffs and write the solution to soln_vec
        dgemv_(&trans, &M, &N, &alpha, sh_matrix, &lda, sh_coeffs, &incx, &beta, soln_vec, &incy);

        pressureGradient(mesh, dvdt, *scalar_dummy, node_num, -factor);

        break;

    }


    delete scalar_dummy;
    delete scalar_lm;
    delete ll_scalar;
    delete scalar_lm_dummy;
    delete[] sh_coeffs;
    // delete[] soln_vec;

}

void velocityDivergence(Mesh * mesh, Array1D<double> & dpdt, Array2D<double> & velocity, double & sum, double h = 1.0)
{
    int node_num, friend_num;
    int i, j, j1;

    Array2D<int> * friend_list;
    Array3D<double> * vel_transform;
    Array3D<double> * div_coeffs;

    double u0, u1, u_temp;
    double v0, v1, v_temp;
    double div_v;
    double cos_a, sin_a;

    node_num = mesh->node_num;
    friend_list = &(mesh->node_friends);
    vel_transform = &(mesh->node_vel_trans);
    div_coeffs = &(mesh->div_coeffs);

    for (i=0; i<node_num; i++)
    {
        friend_num = 6;
        if ((*friend_list)(i, 5) == -1) friend_num = 5;

        u0 = velocity(i,0);
        v0 = velocity(i,1);

        div_v = (*div_coeffs)(i, 0, 0) * u0 + (*div_coeffs)(i, 0, 1) * v0;

        for (j=0; j<friend_num; j++)
        {
            j1 = (*friend_list)(i, j);

            // FIND VELOCITY AT FRIEND j1
            u_temp = velocity(j1, 0);
            v_temp = velocity(j1, 1);

            // CONVERT TO MAPPED VELOCITIES
            cos_a = (*vel_transform)(i, j+1, 0);
            sin_a = (*vel_transform)(i, j+1, 1);

            u1 = u_temp * cos_a + v_temp * sin_a;
            v1 = -u_temp * sin_a + v_temp * cos_a;

            div_v += (*div_coeffs)(i, j+1, 0) * u1 + (*div_coeffs)(i, j+1, 1) * v1;
        }

        dpdt(i) = -h*div_v;

        if (sum >= 0.0) sum += fabs(div_v);
    }
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
