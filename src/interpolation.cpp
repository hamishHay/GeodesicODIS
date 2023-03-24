#include "mesh.h"
#include "globals.h"
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include <math.h>
#include <iostream>
#include "spatialOperators.h"
#include "gridConstants.h"

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/QR>

#include<Eigen/SparseCholesky>
#include<Eigen/SparseQR>
#include <Eigen/OrderingMethods>

// #include <mkl_spblas.h>
// #include <mkl.h>


#include <sys/time.h>


int interpolateVelocity(Globals * globals,
                        Mesh * mesh,
                        Array2D<double> & interp_vel,
                        Array1D<double> & normal_vel)
{
    for (int i=0; i<FACE_NUM; ++i)
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

    return 1;
}

int interpolateVelocityCartRBF(Globals * globals,
                           Mesh * mesh,
                           Array2D<double> & interp_vel_xyz,
                           Array1D<double> & normal_vel)
{
    int friend_num;
    double r = globals->radius.Value();

    double sphvj1[2], sphvj2[2], sphfj[2];
    double sphn[2], sphf[2];

    double nxj, nyj, nzj;

    double njx, njy, xj, yj, mj, v1_xj, v2_xj, v1_yj, v2_yj;
    double ue, ve, ue2, ve2, cosa, sina, ut, vt;
    // double x, y, z, zfact;
    ue = 0.0; ve = 0.0;
    njx = 0.0;
    njy = 0.0;
    cosa = 0.0; sina=0.0;

    // Cartesian components of the edge velocity
    // Array2D<double> ue_xyz(face_num, 3);

    // for (int f=0; f<face_num; f++)
    // {
    //     // Get cartesian coordinates of the face
    //     sphf[0] = mesh->face_centre_pos_sph(f, 0);
    //     sphf[1] = mesh->face_centre_pos_sph(f, 1);

    //     sph2cart(x, y, z, 1.0, pi*0.5 - sphf[0], sphf[1]);

    //     // Caluclate velocity normal to the face in cartesian coordinates

    //     ue = normal_vel(f);

    //     zfact = 1.0 / sqrt(1 - pow(z, 2.0));
    //     ue_yxz(f, 0) = ue*zfact * (-y - z*x);
    //     ue_yxz(f, 1) = ue*zfact * (x - z*y);
    //     ue_yxz(f, 2) = ue*zfact * (1 - pow(z,2.0));
    // }

    // Map Eigen object to the input velocity field, defined normal to faces
    Eigen::Map<Eigen::VectorXd> eig_normal_vel(&normal_vel(0), FACE_NUM);

    // Map Eigen object to the output velocity field components (3), defined at 
    // node centres
    Eigen::Map<Eigen::VectorXd> eig_soln(&interp_vel_xyz(0,0), 3*NODE_NUM);

    // Perform Radial Basis Function interpolation of normal velocities 
    // to cartesian components.

    eig_soln = mesh->operatorRBFinterp * eig_normal_vel;

    //     // Use coefficients to obtain vector components at node i
    //     double vx = 0;
    //     double vy = 0.0;
    //     double vz = 0.0;
    //     double zi = 0;
    //     double xi = 0;
    //     double yi = 0;
    // for (int i=0; i<NODE_NUM; i++) 
    // {

        // xi = mesh->control_vol_edge_centre_pos_map(i, 2, 0)/r;
        // yi = mesh->control_vol_edge_centre_pos_map(i, 2, 1)/r;
    

        //  vx = interp_vel_xyz(i,0);
        //  vy = interp_vel_xyz(i,1);
        //  vz = interp_vel_xyz(i,2);
        // Get cartesian coords of node position
        // sph2cart(xi, yi, zi, 1.0, pi*0.5 - sphn[0], sphn[1]);

        // double lonx, lony, lonz;
        // double latx, laty, latz;
        
        // // make sure that zfact is not infinite
        // if (fabs(zi)-1.0 <= 1e-10) zi = sgn(zi)*(1.0 - 1e-15);
        // double zfact = 1.0/sqrt(1.0-pow(zi,2.0)); 
        // // if (i==0 || i==1) std::cout<<zfact<<' '<<zi<<' '<<1-zi*zi<<std::endl;

        

        // // longitude unit vector, in cartesian components
        // lonx = zfact * -yi;
        // lony = zfact * xi;
        // lonz = 0.0;

        // // latitude unit vector, in cartesian components
        // latx = zfact * -zi*xi;
        // laty = zfact * -zi*yi;
        // latz = zfact * (1-pow(zi,2.0));

        

        // interp_vel(i, 0) = lonx*vx + lony*vy + lonz*vz;
        // interp_vel(i, 1) = latx*vx + laty*vy + latz*vz;

        // if (i<2) std::cout<<0.5*(vx*vx + vy*vy + vz*vz)<<std::endl;

        // interp_KE(i) = 0.5*(vx*vx + vy*vy + vz*vz);



    // }

    

    return 1;
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




/**
*   @purpose    Function to find the flux of a scalar across each edge in 
*               the grid using Radial Basis Function interpolation
*
*   @params
*
*   *globals        pointer to the Globals object in ODIS
*   *mesh           pointer to the Mesh object in ODIS
*   flux            reference to the array to store the calculated flux 
*   edge_vel        reference to the velocity normal to each face
*   node_scalar     reference to the array to interpolate
*   beta            param to control accuracy of interpolation. 1=3rd order,
*                   0 = 4th order.
*
*   @returns
*
*   flux            calculated flux across each face
*
*/
int interpolateLSQFlux(Globals *globals,
                        Mesh *mesh,
                        Array1D<double> & flux,
                        Array1D<double> & edge_vel,
                        Array1D<double> & node_scalar,
                        double beta = 1.0)
{
    double r = globals->radius.Value();
    double r_recip_sq = 1.0/(r*r);
    
    double d2_inner, d2_outer;
    double dx2, dx;
    double fact = 1./12.0;
    double vel;

    int inner_ID;
    int outer_ID;

    // Wrap eigen vector around the scalar vector
    Eigen::Map<Eigen::VectorXd> eig_node_scalar(&node_scalar(0), NODE_NUM);

    // Set beta to 1.0 for third order accurate flux calculation
    beta = 1.0;

    // Find the second direction derivative of s2, normal to the face
    // This is the slowest operation in this function.
    Eigen::VectorXd d2(2*FACE_NUM);
    d2 = mesh->operatorDirectionalSecondDeriv * eig_node_scalar; 

  

    double diff=0.0;
    for (int i=0; i<FACE_NUM; ++i) {
        dx = mesh->face_node_dist(i);
        dx2 = dx*dx*fact;
        vel = edge_vel(i);
        
        inner_ID = mesh->face_nodes(i, 0);
        outer_ID = mesh->face_nodes(i, 1);

        d2_outer = d2(2*i+1);
        d2_inner = d2(2*i);
        
        // Calculate the high order flux (Eq. 11, Skamarock and Gassmann 2011)
        // beta = 1.0 is third order, beta = 0.0 is fourth order
        flux(i) = vel*0.5 * ( node_scalar(inner_ID) + node_scalar(outer_ID) )
                - dx2 * (d2_outer + d2_inner)*vel
                + dx2*beta * fabs(vel) * (d2_outer - d2_inner);


    }
    return 1;

}


      // Loop through every node and solve the least squares
    // inverse problem to find coefficients of the second 
    // order polynomial for the scalar node_salar(i) over 
    // the control volume 
    // Vc = s --> c = V^-1 s, where V is the Vandermonde
    // matrix


    // Array1D<double> c(6*12 + 7*(NODE_NUM-12));
    // for (int i=0; i<6*12 + 7*(NODE_NUM-12); i++) c(i) = 0.0;

    // // for (int i=0; i<NODE_NUM; i++) node_scalar(i) = cos(mesh->node_pos_sph(i,0))*cos(mesh->node_pos_sph(i,0))*cos(2*mesh->node_pos_sph(i,1));
    
    // Eigen::Map<Eigen::VectorXd> cVec(&c(0), 6*12 + 7*(NODE_NUM-12));
    
    // Array2D<double> d2sdx2(NODE_NUM, 3);
    // Array2D<double> d2sdx2_test(NODE_NUM, 3);
    // Eigen::Map<Eigen::VectorXd> d2dx2(&d2sdx2(0,0), 3*NODE_NUM);
    // Eigen::Map<Eigen::VectorXd> d2dx2_test(&d2sdx2_test(0,0), 3*NODE_NUM);

    // d2dx2_test = mesh->operatorSecondDeriv*eig_node_scalar;


    // int count=0;
    // for (int i=0; i<NODE_NUM; ++i) {
        
        
    //     int f_num = 6;
    //     if ( mesh->node_friends(i, 5) == -1 ) f_num--;

    //     Array1D<double> c2(f_num+1);
    //     Eigen::Map<Eigen::VectorXd> cVec2(&c2(0), f_num+1);

    //     // // Build solution vector
    //     Array1D<double> s(f_num+1);

        
    //     s(0) = node_scalar(i);
    //     for (int j=0; j<f_num; ++j) { 
    //         int f_ID = mesh->node_friends(i, j);

    //         s(j+1) = node_scalar(f_ID);   
    //     }

    //     Eigen::Map<Eigen::VectorXd> b(&s(0), f_num+1);
        
    //     Eigen::VectorXd Vec2 = r_recip_sq * mesh->interpSolvers[i].solve(b);
        
    //     d2sdx2(i, 0) = Vec2(3);
    //     d2sdx2(i, 1) = Vec2(4);
    //     d2sdx2(i, 2) = Vec2(5);

    //     // std::cout<<"1: "<<d2sdx2_test(i,0)<<' '<<d2sdx2(i, 0)<<std::endl;
    //     // std::cout<<"2: "<<d2sdx2_test(i,1)<<' '<<d2sdx2(i, 1)<<std::endl;
    //     // std::cout<<"3: "<<d2sdx2_test(i,2)<<' '<<d2sdx2(i, 2)<<std::endl;

    // }

// ------------------------------------------

        // double nx;
        // double ny;

        // double cosa;
        // double sina;

        // double fxx, fyy, fxy;

        // // Get second derivative in direction normal to face
        // double c3, c4, c5;

        // nx = mesh->face_normal_vec_map(i, 0);
        // ny = mesh->face_normal_vec_map(i, 1);

        // cosa = mesh->face_node_vel_trans(i, 0, 0);
        // sina = mesh->face_node_vel_trans(i, 0, 1);

        // // Get second derivative in direction normal to face
        // c3 = d2sdx2(inner_ID, 0);
        // c4 = d2sdx2(inner_ID, 1);
        // c5 = d2sdx2(inner_ID, 2);

        // // c3, c4, c5 are the second derivatives in the xx, xy, and yy 
        // // directions at the node(inner_ID)

        // // Now rotate the 3 components of the second derivative 
        // // to be in map coordinates centered on the face.
        // // This is derived by writing R^T H R, where R is the 
        // // mapping rotation matrix and H is the hermitian matrix 
        // // at the node centre

        // fxx = 2*(c3*cosa*cosa - c4*cosa*sina + c5*sina*sina);
        // fyy = 2*(c3*sina*sina + c4*cosa*sina + c5*cosa*cosa);
        // fxy = 2*(c3 - c5)*cosa*sina + c4*(cosa*cosa - sina*sina);

        // // Now take the directional derivative of the second derivative
        // // to get the second derivative in the direction of the face 
        // // normal vector
        // d2_inner = nx *nx *fxx + 2 *nx *ny * fxy + ny *ny * fyy;

        // // Now repeat the above for the outer node 
        // cosa = mesh->face_node_vel_trans(i, 1, 0);
        // sina = mesh->face_node_vel_trans(i, 1, 1);
        
        // c3 = d2sdx2(outer_ID, 0);
        // c4 = d2sdx2(outer_ID, 1);
        // c5 = d2sdx2(outer_ID, 2);

        // fxx = 2*(c3*cosa*cosa - c4*cosa*sina + c5*sina*sina);
        // fyy = 2*(c3*sina*sina + c4*cosa*sina + c5*cosa*cosa);
        // fxy = 2*(c3 - c5)*cosa*sina + c4*(cosa*cosa - sina*sina);

        // d2_outer = nx *nx *fxx + 2 *nx *ny * fxy + ny *ny * fyy;

        // if (fabs(d2_outer) > 1e-10) {
        //     diff += (d2_outer-d2(2*i+1))/d2_outer + (d2_inner-d2(2*i))/d2_inner;
        //     std::cout<<"1: "<<d2_outer<<' '<<d2(2*i+1)<<std::endl;  
        //     std::cout<<"2: "<<d2_inner<<' '<<d2(2*i)<<std::endl;
        // }
        

        // ----------------------------------------------