#include "mesh.h"
#include "globals.h"
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include <math.h>
#include <iostream>
#include "spatialOperators.h"

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/QR>
// #include <mkl_spblas.h>
// #include <mkl.h>

int interpolateVelocity(Globals * globals,
                        Mesh * mesh,
                        Array2D<double> & interp_vel,
                        Array1D<double> & normal_vel)
{
    int face_num = globals->face_num;

    for (int i=0; i<face_num; ++i)
    {
        double v_tang = 0.0;
        double v_norm = 0.0;
        int friend_num = 10;

        int n1, n2;
        n1 = mesh->face_nodes(i,0);
        n2 = mesh->face_nodes(i,1);
        if (mesh->node_friends(n1,5) < 0) friend_num--;
        if (mesh->node_friends(n2,5) < 0) friend_num--;

        for (int j=0; j<friend_num; j++) {
            int f_ID = mesh->face_interp_friends(i, j);
            v_tang += normal_vel(f_ID)*mesh->face_interp_weights(i,j)*mesh->face_len(f_ID);
        }
        v_tang /= mesh->face_node_dist(i);
        v_norm = normal_vel(i);

        double nx, ny, tx, ty;
        nx = mesh->face_normal_vec_map(i, 0);
        ny = mesh->face_normal_vec_map(i, 1);

        ty = -nx;
        tx = ny;

        interp_vel(i,0) = nx*v_norm + tx*v_tang;
        interp_vel(i,1) = ny*v_norm + ty*v_tang;
    }
}


/**
*   @purpose    Function interpolates data from the geodesic grid to the regular
*               lat-lon grid. Interpolation is biquadratic and takes the form
*               Vc = d, where V is a 6x6 vandermonde matrix, c is a 6x1 column vector
*               of the interpolating function coefficients, and d is a 6x1 column
*               vector of the geodesic grid data to be interpolated. See Lee and
*               Macdonald (2008).
*
*   @params
*
*   *globals        pointer to the Globals object in ODIS
*   *mesh           pointer to the Mesh object in ODIS
*   ll_data         reference to the initialised 2d array for data on the
*                   lat-lon grid
*   gg_data         reference to the initialised geodesic grid data
*   ll_map_coords   reference to the coords of lat-lon nodes in mapping space
*   V_inv           reference to the inverse Vandermonde matrix for every cell
*   cell_ID         reference to the geodesic cell ID for each lat-lon node
*
*   @returns
*
*   ll_data         array with data that has been interpolated from the geodesic
*                   grid to the 2D lat-lon grid
*
*   @author         Hamish Hay
*
*   @history        15/09/2017 - created
*                   05/01/2018 - switched to conservative interpolation
*
*/

int interpolateGG2LLConservative(Globals * globals,
                                 Mesh * mesh,
                                 Array2D<double> & ll_data,
                                 Array1D<double> & gg_data)
{
    double * gg_data_1D;
    double * ll_data_1D;

    int dLat;
    dLat = globals->dLat.Value();

    double r;
    r = globals->radius.Value();

    // No. of nodes on the geodesic and lat-lon grids.
    int node_num_gg = mesh->node_num;
    int node_num_ll = (int)(360/dLat)*(int)(180/dLat);

    // geodesic grid data (0: val, 1: lat gradient, 2: lon gradient)
    gg_data_1D = new double[node_num_gg * 3];

    // -------------------------------------------------------------------------
    // Compute gradient of the solution to be interpolated
    // -------------------------------------------------------------------------

    Array2D<double> * gradient;
    gradient = new Array2D<double>(node_num_gg, 2);

    // pressureGradient(mesh, *gradient, gg_data, node_num_gg, -1.0);

    // -------------------------------------------------------------------------
    // Fill the data vector with the geodesic grid solution and its gradients
    // -------------------------------------------------------------------------

    double * cosLat = &(mesh->trigLat(0,0));
    #pragma omp parallel for
    for (int i=0; i<node_num_gg; i++)
    {
        gg_data_1D[3*i]     = gg_data(i);
        gg_data_1D[3*i + 1] = r*(*gradient)(i, 1);
        gg_data_1D[3*i + 2] = r*(*gradient)(i, 0);
    }


    // -------------------------------------------------------------------------
    // Define variables for matrix-vector multiplication
    // -------------------------------------------------------------------------

    double alpha = 1.0;
    double beta = 1.0;

    // matrix_descr descrp;
    // descrp.type = SPARSE_MATRIX_TYPE_GENERAL;

    // sparse_operation_t operation = SPARSE_OPERATION_NON_TRANSPOSE;

    // sparse_status_t err;

    // sparse_matrix_t dest;
    // sparse_index_base_t index_type = SPARSE_INDEX_BASE_ZERO;
    // int nrows = (360/dLat)*(180/dLat);
    // int ncols = 3*node_num_gg;

    // // -------------------------------------------------------------------------
    // // Multiply the sparse interpolation matrix by the geodesic grid data vector
    // // and return the interpolated solution in ll_data_1D
    // // -------------------------------------------------------------------------

    // err = mkl_sparse_d_mv (operation, alpha, *(mesh->interpMatrix), descrp, gg_data_1D, beta, &(ll_data(0,0)));
    Eigen::Map<Eigen::VectorXd> gg_data_eigen(gg_data_1D, 3*node_num_gg);
    Eigen::Map<Eigen::VectorXd> ll_data_eigen(&ll_data(0,0), node_num_ll);

    // std::cout<<mesh->interpMatrix.cols()<<' '<<3*node_num_gg<<std::endl;
    ll_data_eigen = mesh->interpMatrix*gg_data_eigen;

    delete[] gg_data_1D;
    delete gradient;

    return 1;
}

template <typename T>
int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}





int interpolateLSQFlux(Globals *globals,
                        Mesh *mesh,
                        Array1D<double> & flux,
                        Array1D<double> & edge_vel,
                        Array1D<double> & node_scalar,
                        double beta)
{
    double r = globals->radius.Value();

    int N_face = mesh->face_num;
    int N_node = mesh->node_num;


    Array2D<double> d2sdx2(N_node, 3);
    
    // Loop through every node and solve the least squares
    // inverse problem to find coefficients of the second 
    // order polynomial for the scalar node_salar(i) over 
    // the control volume 
    // Vc = s --> c = V^-1 s, where V is the Vandermonde
    // matrix
    for (int i=0; i<N_node; ++i) {
        
        int f_num = 6;
        if ( mesh->node_friends(i, 5) == -1 ) f_num--;

        // Build solution vector
        Array1D<double> s(f_num+1);
        Array1D<double> c(6);

        s(0) = node_scalar(i);
        for (int j=0; j<f_num; ++j) { 
            int f_ID = mesh->node_friends(i, j);

            s(j+1) = node_scalar(f_ID);    
        }

        Eigen::Map<Eigen::VectorXd> b(&s(0), f_num+1);
        Eigen::Map<Eigen::VectorXd> cVec(&c(0), 6);

        cVec = mesh->interpSolvers[i].solve(b);

        d2sdx2(i, 0) = c(3);
        d2sdx2(i, 1) = c(4);
        d2sdx2(i, 2) = c(5);
    }


    // Then loop through each face, and calculate flux across face 

    for (int i=0; i<N_face; ++i) {
        double d2_inner, d2_outer;
        double dx = mesh->face_node_dist(i);

        int inner_ID = mesh->face_nodes(i, 0);
        int outer_ID = mesh->face_nodes(i, 1);

        double nx = mesh->face_normal_vec_map(i, 0);
        double ny = mesh->face_normal_vec_map(i, 1);

        double cosa = mesh->face_node_vel_trans(i, 0, 0);
        double sina = mesh->face_node_vel_trans(i, 0, 1);

        // Get second derivative in direction normal to face
        double c3, c4, c5;
        c3 = d2sdx2(inner_ID, 0);
        c4 = d2sdx2(inner_ID, 1);
        c5 = d2sdx2(inner_ID, 2);

        // c3, c4, c5 are the second derivatives in the xx, xy, and yy 
        // directions at the node(inner_ID)

        // Now rotate the 3 components of the second derivative 
        // to be in map coordinates centered on the face.
        // This is derived by writing R^T H R, where R is the 
        // mapping rotation matrix and H is the hermitian matrix 
        // at the node centre
        double fxx, fyy, fxy;

        fxx = 2*(c3*cosa*cosa - c4*cosa*sina + c5*sina*sina);
        fyy = 2*(c3*sina*sina + c4*cosa*sina + c5*cosa*cosa);
        fxy = 2*(c3 - c5)*cosa*sina + c4*(cosa*cosa - sina*sina);

        // Now take the directional derivative of the second derivative
        // to get the second derivative in the direction of the face 
        // normal vector
        d2_inner = nx *nx *fxx + 2 *nx *ny * fxy + ny *ny * fyy;


        // Now repeat the above for the outer node 
        cosa = mesh->face_node_vel_trans(i, 1, 0);
        sina = mesh->face_node_vel_trans(i, 1, 1);
        
        c3 = d2sdx2(outer_ID, 0);
        c4 = d2sdx2(outer_ID, 1);
        c5 = d2sdx2(outer_ID, 2);

        fxx = 2*(c3*cosa*cosa - c4*cosa*sina + c5*sina*sina);
        fyy = 2*(c3*sina*sina + c4*cosa*sina + c5*cosa*cosa);
        fxy = 2*(c3 - c5)*cosa*sina + c4*(cosa*cosa - sina*sina);

        d2_outer = nx *nx *fxx + 2 *nx *ny * fxy + ny *ny * fyy;

        // beta = 0.75;


        // Calculate the high order flux (Eq. 11, Skamarock and Gassmann 2011)
        // beta = 1.0 is third order, beta = 0.0 is fourth order
        flux(i) = 0.5 * ( node_scalar(inner_ID) + node_scalar(outer_ID) )
                  - dx*dx/12.0 * (d2_outer + d2_inner);
                   + dx*dx*beta/12.0 * sgn(edge_vel(i)) * (d2_outer - d2_inner);

        flux(i) *= edge_vel(i);
    }

    // Return face flu array

    // globals->Output->TerminateODIS();

    return 1;
}
