#include "mesh.h"
#include "globals.h"
#include "array3d.h"
#include "array2d.h"
#include "array1d.h"
#include "sphericalHarmonics.h"
#include "interpolation.h"
#include "mkl.h"

extern "C"{
// FORTRAN adds _ after all the function names
// and all variables are called by reference

// define FORTRAN matrix-vector multiplication function
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

void pressureGradientN(Mesh * mesh, Array1D<double> & dvdt, Array1D<double> & pressure, int x=1, double g = 1.0)
{
    int node_num, friend_num;
    int i, j, j1, j2,f1,f2, f3;

    int face_num = mesh->face_num;


    for (i=0; i<mesh->face_num; i++) {
        double inner, outer;
        int node_in, node_out;
        int dir_in, dir_out;
        double dist;

        node_in = mesh->face_nodes(i, 0);
        node_out = mesh->face_nodes(i, 1);
        dist = mesh->face_node_dist(i);

        inner = pressure(node_in);
        outer = pressure(node_out);

        double grad = -(inner*mesh->face_centre_m(i, 0) - outer*mesh->face_centre_m(i, 1))/dist;
        //
        // x_grad = g * x_grad / (*cv_areas)(i);
        // y_grad = g * y_grad / (*cv_areas)(i);

        dvdt(i) -= g*grad;
        // dvdt(i,1) -= y_grad;
    }
};

void velocityDivergenceN(Mesh * mesh, Array1D<double> & dpdt, Array1D<double> & velocity_n, double & sum, double h = 1.0)
{
    int node_num, friend_num;
    int i, j, j1, j2, i1, i2, f1, f2, f3;

    Array1D<double> * cv_areas;
    Array2D<int> * friend_list;

    double m;                          // mapping factor at current cv edge
    double div;
    double edge_len;

    node_num = mesh->node_num;
    cv_areas = &(mesh->control_volume_surf_area_map);
    friend_list = &(mesh->node_friends);


    // total_dmass = 0.0;
    for (i=0; i<node_num; i++)
    {
        friend_num = 6;
        if ((*friend_list)(i, 5) == -1) friend_num = 5;

        div = 0.0;

        for (j=0; j<friend_num; j++)
        {

            int face_id = mesh->faces(i,j);

            // get edge length of current edge
            edge_len = mesh->face_len(face_id);
            double dir = mesh->node_face_dir(i, j);

            // calculate control volume divergence
            div += dir * velocity_n(face_id) * edge_len;

        }

        div /= (*cv_areas)(i);
        // if (mesh->cell_is_boundary(i)==1) std::cout<<i<<div<<std::endl;

        dpdt(i) = -h*div;


        if (sum >= 0.0) {
            sum += div;
            // total_dmass += dmass/(*mass)(i);
        }
    }

    // std::cout<<"TOTAL: "<<sum<<std::endl;
    sum = 0.0;
    // std::cout<<std::scientific<<total_dmass<<std::endl;
};

// void velocityDivergence(Mesh * mesh, Array1D<double> & dpdt, Array2D<double> & velocity, double & sum, double h = 1.0)
// {
//     int node_num, friend_num;
//     int i, j, j1, j2, i1, i2, f1, f2, f3;
//
//     Array2D<int> * friend_list;
//     Array2D<double> * cent_map_factor;
//     Array3D<double> * element_areas;
//     Array3D<double> * normal_vecs;
//     Array2D<double> * edge_lens;
//     Array1D<double> * cv_areas;
//     Array3D<double> * vel_transform;
//     // Array1D<double> * mass;
//     Array1D<int> * mask;
//
//     double m;                          // mapping factor at current cv edge
//     double a0, a1, a2;
//     double u0, u1, u2, u0_cent, u1_cent, u_avg, u_temp;
//     double v0, v1, v2, v0_cent, v1_cent, v_avg, v_temp;
//     double div;
//     double nx, ny;
//     double edge_len;
//     // double total_dmass, dmass, dt;
//     // double h;
//     double cos_a, sin_a;
//
//     node_num = mesh->node_num;
//     friend_list = &(mesh->node_friends);
//     cent_map_factor = &(mesh->control_vol_edge_centre_m);
//     element_areas = &(mesh->node_friend_element_areas_map);
//     normal_vecs = &(mesh->control_vol_edge_normal_map);
//     edge_lens = &(mesh->control_vol_edge_len);
//     cv_areas = &(mesh->control_volume_surf_area_map);
//     // mass = &(mesh->control_volume_mass);
//     // dt = mesh->globals->timeStep.Value();
//     vel_transform = &(mesh->node_vel_trans);
//     // mask = &(mesh->land_mask);
//
//     // total_dmass = 0.0;
//     for (i=0; i<node_num; i++)
//     {
//         friend_num = 6;
//         if ((*friend_list)(i, 5) == -1) friend_num = 5;
//
//         u0 = velocity(i,0);
//         v0 = velocity(i,1);
//
//         div = 0.0;
//         // dmass = 0.0;
//
//         for (j=0; j<friend_num; j++)
//         {
//             // find avg pressure in element j
//             j1 = j%friend_num;
//             j2 = (j+1)%friend_num;
//             i1 = (*friend_list)(i,j1);
//             i2 = (*friend_list)(i,j2);
//
//             f1 = (*friend_list)(i, j%friend_num);
//             f2 = (*friend_list)(i, (j+1)%friend_num);
//             f3 = (*friend_list)(i, (j+2)%friend_num);
//
//
//             a0 = (*element_areas)(i,j1,0);
//             a1 = (*element_areas)(i,j1,1);
//             a2 = (*element_areas)(i,j1,2);
//
//             // FIND VELOCITY AT FIRST FRIEND
//
//             u_temp = velocity(i1,0);
//             v_temp = velocity(i1,1);
//
//             // CONVERT TO MAPPED VELOCITIES
//
//             cos_a = (*vel_transform)(i, j1+1, 0);
//             sin_a = (*vel_transform)(i, j1+1, 1);
//             u1 = u_temp * cos_a + v_temp * sin_a;
//             v1 = -u_temp * sin_a + v_temp * cos_a;
//
//             // FIND VELOCITY AT SECOND FRIEND
//
//             u_temp = velocity(i2,0);
//             v_temp = velocity(i2,1);
//
//             // CONVERT TO MAPPED VELOCITIES
//
//             cos_a = (*vel_transform)(i, j2+1, 0);
//             sin_a = (*vel_transform)(i, j2+1, 1);
//             u2 = u_temp * cos_a + v_temp * sin_a;
//             v2 = -u_temp * sin_a + v_temp * cos_a;
//
//             // FIND AVERAGE VELOCITY AT FIRST ELEMENT CENTRE
//             u0_cent = (u0 * a1 + u1 * a2 + u2 * a0) / (a0 + a1 + a2);
//             v0_cent = (v0 * a1 + v1 * a2 + v2 * a0) / (a0 + a1 + a2);
//
//             // find avg pressure in element j+1
//             j1 = j2;
//             j2 = (j+2)%friend_num;
//             i1 = (*friend_list)(i,j1);
//             i2 = (*friend_list)(i,j2);
//
//             a0 = (*element_areas)(i,j1,0);
//             a1 = (*element_areas)(i,j1,1);
//             a2 = (*element_areas)(i,j1,2);
//
//             // FIND VELOCITY AT FIRST FRIEND
//             u_temp = velocity(i1,0);
//             v_temp = velocity(i1,1);
//
//             // CONVERT TO MAPPED VELOCITIES
//             cos_a = (*vel_transform)(i, j1+1, 0);
//             sin_a = (*vel_transform)(i, j1+1, 1);
//             u1 = u_temp * cos_a + v_temp * sin_a;
//             v1 = -u_temp * sin_a + v_temp * cos_a;
//
//             // FIND VELOCITY AT SECOND FREIND
//             u_temp = velocity(i2,0);
//             v_temp = velocity(i2,1);
//
//             // CONVERT TO MAPPED VELOCITIES
//             cos_a = (*vel_transform)(i, j2+1, 0);
//             sin_a = (*vel_transform)(i, j2+1, 1);
//             u2 = u_temp * cos_a + v_temp * sin_a;
//             v2 = -u_temp * sin_a + v_temp * cos_a;
//
//             // FIND AVERAGE VELOCITY AT FIRST ELEMENT CENTRE
//             u1_cent = (u0 * a1 + u1 * a2 + u2 * a0) / (a0 + a1 + a2);
//             v1_cent = (v0 * a1 + v1 * a2 + v2 * a0) / (a0 + a1 + a2);
//
//             u_temp = velocity(f2,0);
//             v_temp = velocity(f2,1);
//
//
//             // Find average p at the center of the control volume edge
//             v_avg = 0.5*(v0_cent + v1_cent);
//             u_avg = 0.5*(u0_cent + u1_cent);
//
//             if (mesh->cell_is_boundary(i)==1)
//             {
//               // std::cout<<i<<' '<<f1<<' '<<f2<<' '<<f3<<std::endl;
//               if (mesh->cell_is_boundary(f1)==2 || mesh->cell_is_boundary(f3)==2)
//               {
//                 // u_avg = 0.5 * (u0 + velocity(f2, 0));
//                 // v_avg = 0.5 * (v0 + velocity(f2, 1));
//                 u_avg = 0.5 * (velocity(f2, 0) - u0) + u0;
//                 v_avg = 0.5 * (velocity(f2, 1) - v0) + v0;
//                 // p_avg = 0.5 * (pressure(f2) - p0) + p0;
//               }
//
//             }
//
//             i1 = (*friend_list)(i,(j+1)%friend_num);
//             if (mesh->cell_is_boundary(i)==1 && mesh->cell_is_boundary(i1)==2)
//             {
//               u_avg = 0.0;
//               v_avg = 0.0;
//             }
//
//             j1 = j%friend_num;
//
//             // get mapping factor for edge i,j2
//             m = (*cent_map_factor)(i, j1);
//
//             // get components of the edge normal vector
//             nx = (*normal_vecs)(i, j1, 0);
//             ny = (*normal_vecs)(i, j1, 1);
//
//             // get edge length of current edge
//             edge_len = (*edge_lens)(i,j1);
//
//             // if (i==973 || i==974) std::cout<<i<<' '<<f2<<' '<<nx<<' '<<ny<<' '<<(*edge_lens)(i,j%(friend_num))<<std::endl;
//
//
//             // double flux = ((u_avg * nx) + (v_avg * ny)) * edge_len / m;
//
//             // if (i==137 && f1==1427 && f2==1453 && f3 ==1457)
//             // {
//             //   std::cout<<i<<' '<<flux<<' '<<f1<<' '<<f2<<' '<<f3<<std::endl;
//             //   // for (k=0; k<3; k++)
//             //   // {
//             //   //   a0 = (*element_areas)(i,j,0);
//             //   //   a1 = (*element_areas)(i,j,1);
//             //   //   a2 = (*element_areas)(i,j,2);
//             //   // }
//             //   // std::cout<<i<<' '<<a0<<' '<<a1<<' '<<a2<<std::endl;
//             //   // std::cout<<i<<' '<<nx<<' '<<ny<<std::endl;
//             //   // std::cout<<i<<' '<<m<<std::endl;
//             //   // std::cout<<i<<' '<<mesh->node_dists(i, (j+1)%friend_num)<<std::endl;
//             //   // for (int k=0; k<3; k++) std::cout<<' '<<f1<<' '<<(*element_areas)(i,(j)%friend_num,k);
//             //   // for (int k=0; k<3; k++) std::cout<<' '<<f2<<' '<<(*element_areas)(i,(j+1)%friend_num,k);
//             //   // for (int k=0; k<3; k++) std::cout<<' '<<f2<<' '<<(*element_areas)(i,(j+2)%friend_num,k);
//             // }
//             // if (i==1453 && f1==1457 && f2==137 && f3 ==1427)
//             // {
//             //   std::cout<<i<<' '<<flux<<' '<<f1<<' '<<f2<<' '<<f3<<std::endl;
//             //   // std::cout<<i<<' '<<a0<<' '<<a1<<' '<<a2<<std::endl;
//             //   // std::cout<<i<<' '<<nx<<' '<<ny<<std::endl;
//             //   // std::cout<<i<<' '<<m<<std::endl;
//             //   // std::cout<<i<<' '<<mesh->node_dists(i, (j+1)%friend_num)<<std::endl;
//             //   // for (int k=0; k<3; k++) std::cout<<' '<<f1<<' '<<(*element_areas)(i,(j)%friend_num,k);
//             //   // for (int k=0; k<3; k++) std::cout<<' '<<f2<<' '<<(*element_areas)(i,(j+1)%friend_num,k);
//             //   // for (int k=0; k<3; k++) std::cout<<' '<<f2<<' '<<(*element_areas)(i,(j+2)%friend_num,k);
//             // }
//             // calculate control volume divergence
//             div += ((u_avg * nx) + (v_avg * ny)) * edge_len / m;
//
//         }
//
//
//         // dmass = dt * div * 1000.0 * h;
//
//
//         div /= (*cv_areas)(i);
//         // if (mesh->cell_is_boundary(i)==1) std::cout<<i<<div<<std::endl;
//
//         dpdt(i) = -h*div;
//
//
//         if (sum >= 0.0) {
//             sum += div;
//             // total_dmass += dmass/(*mass)(i);
//         }
//     }
//
//     std::cout<<"TOTAL: "<<sum<<std::endl;
//     sum = 0.0;
//     // std::cout<<std::scientific<<total_dmass<<std::endl;
// };

// void pressureGradient(Mesh * mesh, Array2D<double> & dvdt, Array1D<double> & pressure, int nodeNum, double g = 1.0)
// {
//     int i, j;
//
//     sparse_operation_t operation = SPARSE_OPERATION_NON_TRANSPOSE;
//     matrix_descr descript;
//
//     descript.type = SPARSE_MATRIX_TYPE_GENERAL;
//     double alpham = 1.0;
//     double betam = 1.0;
//
//     double * vec = new double[3*nodeNum];
//
//     #pragma omp parallel for
//     for (i=0; i<nodeNum; i++)
//     {
//         vec[3*i] = 0.0;
//         vec[3*i+1] = 0.0;
//         vec[3*i+2] = pressure(i);
//     }
//
//     mkl_sparse_d_mv(operation, g, *(mesh->operatorGradient), descript, vec, betam, &(dvdt(0,0)));
//
//     // for (i=0; i<nodeNum; i++)
//     // {
//     //     friend_num = 6;
//     //     if ((*friend_list)(i, 5) == -1) friend_num = 5;
//     //
//     //     p0 = pressure(i);
//     //
//     //     x_grad = p0 * (*grad_coeffs)(i, 0, 0);
//     //     y_grad = p0 * (*grad_coeffs)(i, 0, 1);
//     //
//     //     for (j=0; j<friend_num; j++)
//     //     {
//     //         p1 = pressure( (*friend_list)(i, j) );
//     //
//     //         x_grad += p1 * (*grad_coeffs)(i, j+1, 0);
//     //         y_grad += p1 * (*grad_coeffs)(i, j+1, 1);
//     //     }
//     //
//     //     dvdt(i,0) -= g*x_grad;
//     //     dvdt(i,1) -= g*y_grad;
//     // }
//
//     delete[] vec;
// };


// TODO - move this function to a new file -- maybe selfGravity.cpp ?
void pressureGradientSH(Globals * globals, Mesh * mesh, Array1D<double> & dvdt, Array1D<double> & gg_scalar, Array1D<double> & gg_soln, double factor)
{

    int node_num, l_max, N_ll;
    int i, j, f, l, m;
    int start_l;
    double * sh_coeffs;

    Array3D<double> * scalar_lm;
    Array2D<double> * ll_scalar;

    node_num = globals->node_num;
    l_max = globals->l_max.Value();
    N_ll = (int)globals->dLat.Value();

    scalar_lm = new Array3D<double>(2.*(l_max+1), 2.*(l_max+1), 2);

    ll_scalar = new Array2D<double>(180/N_ll, 360/N_ll);

    if (globals->surface_type == LID_LOVE) start_l = 2;
    if (globals->surface_type == FREE_LOADING) start_l = 2;
    start_l = 0;

    sh_coeffs = new double[(l_max+1)*(l_max+2) - 6];

    // getSHCoeffsGG(mesh->node_pos_sph, gg_scalar, *scalar_lm, node_num, l_max);

    interpolateGG2LLConservative(globals,
                                mesh,
                                *ll_scalar,
                                gg_scalar);

    // FIND SPHERICAL HARMONIC EXPANSION COEFFICIENTS
    // OF THE PRESSURE FIELD
    getSHCoeffsLL(*ll_scalar, *scalar_lm, N_ll, l_max);

    double * sh_matrix;
    double * soln_vec;
    sh_matrix = &(mesh->sh_matrix_fort(0));
    soln_vec = &(gg_soln(0));

    char trans = 'n';
    int M = node_num;
    int N = (l_max+1)*(l_max+2) - 6;
    double alpha = factor;
    int lda = M;
    int incx = 1;
    int incy = 1;
    double beta = 1.0;

    int count = 0;
    for (l=2; l<l_max+1; l++)
    {
        for (m=0; m<=l; m++)
        {
            sh_coeffs[count] = (*scalar_lm)(l, m, 0);
            count++;

            sh_coeffs[count] = (*scalar_lm)(l, m, 1);
            count++;

            // std::cout<<l<<' '<<m<<' '<<(*scalar_lm)(l, m, 0)<<' '<<(*scalar_lm)(l, m, 1)<<std::endl;
        }
    }

    // Perform sh_matrix * sh_coeffs and write the solution to soln_vec
    dgemv_(&trans, &M, &N, &alpha, sh_matrix, &lda, sh_coeffs, &incx, &beta, soln_vec, &incy);

    delete scalar_lm;
    delete ll_scalar;

    delete[] sh_coeffs;

};

// void velocityDivergence(Mesh * mesh, Array1D<double> & dpdt, Array2D<double> & velocity, double & sum, double h = 1.0)
// {
//     int node_num, friend_num;
//     int i, j, j1;
//
//
//     node_num = mesh->node_num;
//
//     sparse_operation_t operation = SPARSE_OPERATION_NON_TRANSPOSE;
//     matrix_descr descript;
//
//     descript.type = SPARSE_MATRIX_TYPE_GENERAL;
//     double alpham = 1.0;
//     double betam = 0.0;
//
//     double * vec = new double[3*node_num];
//
//     #pragma omp parallel for
//     for (i=0; i<node_num; i++)
//     {
//         vec[3*i] = velocity(i,0);
//         vec[3*i+1] = velocity(i,1);
//         vec[3*i+2] = 0.0;
//     }
//
//     mkl_sparse_d_mv(operation, -h, *(mesh->operatorDivergence), descript, vec, betam, &(dpdt(0)));
//
//     delete[] vec;
//     //
//     // for (i=0; i<node_num; i++)
//     // {
//     //     friend_num = 6;
//     //     if ((*friend_list)(i, 5) == -1) friend_num = 5;
//     //
//     //     u0 = velocity(i,0);
//     //     v0 = velocity(i,1);
//     //
//     //     div_v = (*div_coeffs)(i, 0, 0) * u0 + (*div_coeffs)(i, 0, 1) * v0;
//     //
//     //     for (j=0; j<friend_num; j++)
//     //     {
//     //         j1 = (*friend_list)(i, j);
//     //
//     //         // FIND VELOCITY AT FRIEND j1
//     //         u_temp = velocity(j1, 0);
//     //         v_temp = velocity(j1, 1);
//     //
//     //         // CONVERT TO MAPPED VELOCITIES
//     //         cos_a = (*vel_transform)(i, j+1, 0);
//     //         sin_a = (*vel_transform)(i, j+1, 1);
//     //
//     //         u1 = u_temp * cos_a + v_temp * sin_a;
//     //         v1 = -u_temp * sin_a + v_temp * cos_a;
//     //
//     //         div_v += (*div_coeffs)(i, j+1, 0) * u1 + (*div_coeffs)(i, j+1, 1) * v1;
//     //     }
//     //
//     //     // dpdt(i) = -h*div_v;
//     //
//     //     if (sum >= 0.0) sum += fabs(div_v);
//     //     dpdt(i) += test[i];
//     //
//     //     if (fabs(-h*div_v-test[i])>1e-12) std::cout<<i<<' '<<-h*div_v<<' '<<test[i]<<' '<<-h*div_v-test[i]<<std::endl;
//     // }
//
//     // mesh->globals->Output->TerminateODIS();
//
//     // delete[] test;
// };

// Function to calculate the Laplacian and diffusion for the velocity field.
// The discretized operator used here is the second order accurate finite
// difference operator given as equation 3.2 in Heikes et al (2013).
void velocityDiffusion(Mesh * mesh, Array2D<double> & dvdt, Array2D<double> & velocity, double viscosity)
{
    int i, node_num;
    // int node_num, friend_num;
    // int i, j, i1, j1, j2;
    //
    // Array2D<int> * friend_list;
    // Array2D<double> * edge_lens;
    // Array1D<double> * cv_areas;
    // Array3D<double> * vel_transform;
    // Array2D<double> * node_friend_dists;
    //
    // double m;                          // mapping factor at current cv edge
    // double u0, u1, u_temp;
    // double v0, v1, v_temp;
    // double edge_len;
    // double cos_a, sin_a;
    // double lap_u, lap_v;
    // double node_friend_dist;
    //
    node_num = mesh->node_num;
    // friend_list = &(mesh->node_friends);
    // edge_lens = &(mesh->control_vol_edge_len);
    // cv_areas = &(mesh->control_volume_surf_area_map);
    // vel_transform = &(mesh->node_vel_trans);
    // node_friend_dists = &(mesh->node_dists);

    sparse_operation_t operation = SPARSE_OPERATION_NON_TRANSPOSE;
    matrix_descr descript;
    // struct matrix_descr {
    //     sparse_matrix_type_t type;
    // } descript;
    descript.type = SPARSE_MATRIX_TYPE_GENERAL;
    double alpham = 1.0;
    double betam = 1.0;

    double * vec = new double[3*node_num];
    for (i=0; i<node_num; i++)
    {
        vec[3*i] = velocity(i,0);
        vec[3*i+1] = velocity(i,1);
        vec[3*i+2] = 0.0;
    }

    mkl_sparse_d_mv(operation, viscosity, *(mesh->operatorLaplacian2), descript, vec, betam, &(dvdt(0,0)));

    delete[] vec;

    // for (i=0; i<node_num; i++)
    // {
    //     friend_num = 6;
    //     if ((*friend_list)(i, 5) == -1) friend_num = 5;
    //
    //     u0 = velocity(i,0);
    //     v0 = velocity(i,1);
    //
    //     lap_u = 0.0;
    //     lap_v = 0.0;
    //
    //     for (j=0; j<friend_num; j++)
    //     {
    //           i1 = (*friend_list)(i,j);
    //
    //           j1 = (j-1)%friend_num;
    //           if (j1 < 0) j1 += friend_num;
    //
    //           u_temp = velocity(i1,0);
    //           v_temp = velocity(i1,1);
    //
    //           cos_a = (*vel_transform)(i, j+1, 0);
    //           sin_a = (*vel_transform)(i, j+1, 1);
    //
    //           u1 = u_temp * cos_a + v_temp * sin_a;
    //           v1 = -u_temp * sin_a + v_temp * cos_a;
    //
    //           // CHANGED - index of j is now less by 1. The distance between node
    //           //           and friend is stored in the same index system as the
    //           //           friends themselves. The control volume edge lengths are
    //           //           stored in the j-1 index system (i.e., edge s6 lines up
    //           //           with node-friend length l1)
    //
    //           node_friend_dist = (*node_friend_dists)(i, j);  // same j index as the friend node
    //           edge_len = (*edge_lens)(i, j1);                 // index j-1 due to clockwise storage
    //
    //           lap_u += (u1 - u0)/node_friend_dist * edge_len;
    //           lap_v += (v1 - v0)/node_friend_dist * edge_len;
    //
    //     }
    //
    //     dvdt(i,0) += viscosity * lap_u/(*cv_areas)(i);
    //     dvdt(i,1) += viscosity * lap_v/(*cv_areas)(i);
    //
    // }
};

// Function to calculate the Laplacian and diffusion for the velocity field.
// The discretized operator used here is the second order accurate finite
// difference operator given as equation 3.2 in Heikes et al (2013).
void scalarDiffusion(Mesh * mesh, Array1D<double> & d2s, Array1D<double> & s, double viscosity)
{
    int node_num, friend_num;
    int i, j, i1, j1;

    Array2D<int> * friend_list;
    Array2D<double> * edge_lens;
    Array1D<double> * cv_areas;
    Array2D<double> * node_friend_dists;

    double m;                          // mapping factor at current cv edge
    double s0, s1, s_temp;
    double edge_len;
    double cos_a, sin_a;
    double lap_s, lap_v;
    double node_friend_dist;

    node_num = mesh->node_num;
    friend_list = &(mesh->node_friends);
    edge_lens = &(mesh->control_vol_edge_len);
    cv_areas = &(mesh->control_volume_surf_area_map);
    node_friend_dists = &(mesh->node_dists);

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

// void smoothingSH(Globals * globals, Mesh * mesh, Array1D<double> & gg_scalar)
// {
//     int node_num, l_max, N_ll;
//     int i, j, f, l, m;
//     double r;
//     double interp_val;
//     double * cosLat;
//
//     Array3D<double> * scalar_lm;
//
//     Array3D<double> * Pbar_lm;
//     Array3D<double> * Pbar_lm_deriv;
//     Array3D<double> * trigMLon;
//     Array2D<double> * trigLat;
//
//     Array2D<double> * ll_scalar;
//
//     node_num = globals->node_num;
//     l_max = globals->l_max.Value();
//     N_ll = (int)globals->dLat.Value();
//     r = globals->radius.Value();
//
//     scalar_lm = new Array3D<double>(2.*(l_max+1), 2.*(l_max+1), 2);
//     ll_scalar = new Array2D<double>(180/N_ll, 360/N_ll);
//     // scalar_dummy = new Array1D<double>(node_num);
//
//     Pbar_lm = &(mesh->Pbar_lm);             // 4-pi Normalised associated leg funcs
//     // Pbar_lm_deriv = &(mesh->Pbar_lm_deriv); // 4-pi Normalised associated leg funcs derivs
//     trigMLon = &(mesh->trigMLon);           // cos and sin of m*longitude
//     trigLat = &(mesh->trigLat);
//     cosLat = &(mesh->trigLat(0,0));
//
//     // INTERPOLATE GEODESIC GRID DATA TO LAT-LON GRID
//     // interpolateGG2LL(globals,
//     //                     mesh,
//     //                     *ll_scalar,
//     //                     gg_scalar,
//     //                     mesh->ll_map_coords,
//     //                     mesh->V_inv,
//     //                     mesh->cell_ID);
//
//     interpolateGG2LLConservative(globals,
//                             mesh,
//                             *ll_scalar,
//                             gg_scalar);
//
//     // FIND SPHERICAL HARMONIC EXPANSION COEFFICIENTS
//     // OF THE PRESSURE FIELD
//     getSHCoeffsLL(*ll_scalar, *scalar_lm, N_ll, l_max);
//
//     for (i=0; i<node_num; i++)
//     {
//         interp_val = 0.0;
//         for (l=0; l<l_max+1; l++)
//         {
//             for (m=0; m<=l; m++)
//             {
//                 interp_val += (*Pbar_lm)(i, l, m)
//                           * ((*scalar_lm)(l, m, 0) * (*trigMLon)(i, m, 0)
//                           + (*scalar_lm)(l, m, 1) * (*trigMLon)(i, m, 1));
//             }
//         }
//
//         gg_scalar(i) = interp_val;
//     }
//
//     delete scalar_lm;
//     delete ll_scalar;
//
// };

// void smoothingSHVector(Globals * globals, Mesh * mesh, Array2D<double> & gg_vector)
// {
//     int node_num, l_max, N_ll;
//     int i, j, f, l, m;
//     double r;
//     double interp_val;
//     double * cosLat;
//
//     Array3D<double> * scalar_lm;
//
//     Array3D<double> * Pbar_lm;
//     Array3D<double> * Pbar_lm_deriv;
//     Array3D<double> * trigMLon;
//     Array2D<double> * trigLat;
//
//     Array2D<double> * ll_scalar;
//     Array1D<double> * gg_scalar;
//
//
//     node_num = globals->node_num;
//     l_max = globals->l_max.Value();
//     N_ll = (int)globals->dLat.Value();
//     r = globals->radius.Value();
//
//     scalar_lm = new Array3D<double>(2.*(l_max+1), 2.*(l_max+1), 2);
//     ll_scalar = new Array2D<double>(180/N_ll, 360/N_ll);
//     gg_scalar = new Array1D<double>(node_num);
//
//     Pbar_lm = &(mesh->Pbar_lm);             // 4-pi Normalised associated leg funcs
//     Pbar_lm_deriv = &(mesh->Pbar_lm_deriv); // 4-pi Normalised associated leg funcs derivs
//     trigMLon = &(mesh->trigMLon);           // cos and sin of m*longitude
//     trigLat = &(mesh->trigLat);
//     cosLat = &(mesh->trigLat(0,0));
//
//     for (j=0; j<2; j++)
//     {
//         for (i=0; i<node_num; i++)
//         {
//             (*gg_scalar)(i) = gg_vector(i, j);
//         }
//
//         // INTERPOLATE GEODESIC GRID DATA TO LAT-LON GRID
//         interpolateGG2LL(globals,
//                         mesh,
//                         *ll_scalar,
//                         *gg_scalar,
//                         mesh->ll_map_coords,
//                         mesh->V_inv,
//                         mesh->cell_ID);
//
//         // FIND SPHERICAL HARMONIC EXPANSION COEFFICIENTS
//         // OF THE PRESSURE FIELD
//         getSHCoeffsLL(*ll_scalar, *scalar_lm, N_ll, l_max);
//
//         for (i=0; i<node_num; i++)
//         {
//             interp_val = 0.0;
//             for (l=0; l<l_max+1; l++)
//             {
//                 for (m=0; m<=l; m++)
//                 {
//                     interp_val += (*Pbar_lm)(i, l, m)
//                     * ((*scalar_lm)(l, m, 0) * (*trigMLon)(i, m, 0)
//                     + (*scalar_lm)(l, m, 1) * (*trigMLon)(i, m, 1));
//                 }
//             }
//
//             gg_vector(i, j) = interp_val;
//         }
//     }
//
//
//     delete scalar_lm;
//     delete ll_scalar;
//     delete gg_scalar;
//
// };

// TODO - move function to interpolation.cpp
void avgAtPoles(Mesh * mesh, Array2D<double> & velocity)
{
    int friend_num;
    int i, j, j1, j2, i1, i2;

    Array2D<int> * friend_list;
    Array3D<double> * vel_transform;

    double m;                          // mapping factor at current cv edge
    double a0, a1, a2;
    double u0, u1, u2, u0_cent, u1_cent, u_avg, u_temp;
    double v0, v1, v2, v0_cent, v1_cent, v_avg, v_temp;
    double div;
    double nx, ny;
    double edge_len;

    double cos_a, sin_a;

    friend_list = &(mesh->node_friends);

    vel_transform = &(mesh->node_vel_trans);

    for (i=0; i<2; i++)
    {
        friend_num = 5;

        u_avg = 0.0;
        v_avg = 0.0;
        for (j=0; j<friend_num; j++)
        {
            // FIND VELOCITY AT FRIEND
            i1 = (*friend_list)(i, j);

            u_temp = velocity(i1, 0);
            v_temp = velocity(i1, 1);

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

// void pressureGradient(Mesh * mesh, Array2D<double> & dvdt, Array1D<double> & pressure, int x=1, double g = 1.0)
// {
//     int node_num, friend_num;
//     int i, j, j1, j2,f1,f2, f3;
//
//     Array2D<int> * friend_list;
//     Array2D<double> * cent_map_factor;
//     Array3D<double> * element_areas;
//     Array3D<double> * normal_vecs;
//     Array2D<double> * edge_lens;
//     Array1D<double> * cv_areas;
//     Array1D<int> * mask;
//
//     double m;                          // mapping factor at current cv edge
//     double a0, a1, a2;
//     double p0, p1, p2, p0_cent, p1_cent, p_avg;
//     double x_grad, y_grad;
//     double nx, ny;
//     double edge_len;
//     // double g;
//
//     node_num = mesh->node_num;
//     friend_list = &(mesh->node_friends);
//     cent_map_factor = &(mesh->control_vol_edge_centre_m);
//     element_areas = &(mesh->node_friend_element_areas_map);
//     normal_vecs = &(mesh->control_vol_edge_normal_map);
//     edge_lens = &(mesh->control_vol_edge_len);
//     cv_areas = &(mesh->control_volume_surf_area_map);
//     // g = mesh->globals->g.Value();
//     // mask = &(mesh->land_mask);
//
//     for (i=0; i<node_num; i++)
//     {
//         friend_num = 6;
//         if ((*friend_list)(i, 5) == -1) friend_num = 5;
//
//         p0 = pressure(i);
//
//         x_grad = 0.0;
//         y_grad = 0.0;
//
//         for (j=0; j<friend_num; j++)
//         {
//             // find avg pressure in element j
//             j1 = j%friend_num;
//             j2 = (j+1)%friend_num;
//
//             f1 = (*friend_list)(i, j%friend_num);
//             f2 = (*friend_list)(i, (j+1)%friend_num);
//             f3 = (*friend_list)(i, (j+2)%friend_num);
//
//             a0 = (*element_areas)(i,j1,0);
//             a1 = (*element_areas)(i,j1,1);
//             a2 = (*element_areas)(i,j1,2);
//
//             p1 = pressure((*friend_list)(i,j1));
//             p2 = pressure((*friend_list)(i,j2));
//
//             // if (mesh->cell_is_boundary(i)==1 && mesh->cell_is_boundary(j1)==2)
//             // {
//             //   p1 = 1/a2 * ((a0+a1+a2)*0.5*(p0+p2) -a1*p0 - a0*p2);
//             // }
//             //
//             // if (mesh->cell_is_boundary(i)==1 && mesh->cell_is_boundary(j2)==2)
//             // {
//             //   p2 = 1/a0 * ((a0+a1+a2)*0.5*(p0+p1) -a1*p0 - a2*p1);
//             // }
//
//             p0_cent = (p0 * a1 + p1 * a2 + p2 * a0) / (a0 + a1 + a2);
//
//
//
//             // find avg pressure in element j+1
//             j1 = j2;
//             j2 = (j+2)%friend_num;
//
//             a0 = (*element_areas)(i,j1,0);
//             a1 = (*element_areas)(i,j1,1);
//             a2 = (*element_areas)(i,j1,2);
//
//             p1 = pressure((*friend_list)(i,j1));
//             p2 = pressure((*friend_list)(i,j2));
//
//             // if (mesh->cell_is_boundary(i)==1 && mesh->cell_is_boundary(j1)==2)
//             // {
//             //   p1 = 1/a2 * ((a0+a1+a2)*0.5*(p0+p2) -a1*p0 - a0*p2);
//             // }
//             //
//             // if (mesh->cell_is_boundary(i)==1 && mesh->cell_is_boundary(j2)==2)
//             // {
//             //   p2 = 1/a0 * ((a0+a1+a2)*0.5*(p0+p1) -a1*p0 - a2*p1);
//             // }
//
//             p1_cent = (p0*a1 + p1*a2 + p2*a0) / (a0 + a1 + a2);
//
//             // Find average p at the center of the control volume edge
//             p_avg = 0.5*(p0_cent + p1_cent);
//
//             if (mesh->cell_is_boundary(i)==1)
//             {
//               if (mesh->cell_is_boundary(f1)==2 || mesh->cell_is_boundary(f3)==2)
//               {
//                 p_avg = 0.5 * (pressure(f2) - p0) + p0;
//               }
//               else if (mesh->cell_is_boundary(f2)==2)
//               {
//                 p_avg = p0;
//               }
//
//             }
//
//             // if (mesh->cell_is_boundary(i)==1 && mesh->cell_is_boundary(j2)==2)
//             // {
//             //   p2 = 1/a0 * ((a0+a1+a2)*0.5*(p0+p1) -a1*p0 - a2*p1);
//             // }
//
//             // if (mesh->cell_is_boundary(i)==1 && mesh->cell_is_boundary(j1)==2)
//             // {
//             //   p_avg = 0.5*(p0 + p2);
//             // }
//
//             j1 = j%friend_num;
//
//             // get mapping factor for edge i,j2
//             m = (*cent_map_factor)(i, j1);
//
//             // get components of the edge normal vector
//             nx = (*normal_vecs)(i, j1, 0);
//             ny = (*normal_vecs)(i, j1, 1);
//
//             // get edge length of current edge
//             edge_len = (*edge_lens)(i,j1);
//
//             // calculate x gradient
//             x_grad += m * p_avg * nx * edge_len;
//
//             // calculate y gradient
//             y_grad += m * p_avg * ny * edge_len;
//
//             j1 = (*friend_list)(i,(j+5)%friend_num);
//             // if (mesh->cell_is_boundary(i)==1 && mesh->cell_is_boundary(j1)==2)
//             // {
//             //   x_grad = 0.0;
//             //   y_grad = 0.0;
//             // }
//
//         }
//
//         x_grad = g * x_grad / (*cv_areas)(i);
//         y_grad = g * y_grad / (*cv_areas)(i);
//
//         dvdt(i,0) -= x_grad;
//         dvdt(i,1) -= y_grad;
//     }
// };
