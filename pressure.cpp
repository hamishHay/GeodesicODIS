
#include "mesh.h"
#include "globals.h"
#include "outFiles.h"
#include "solver.h"
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include "pressure.h"
#include "spatialOperators.h"
#include "sphericalHarmonics.h"
#include "interpolation.h"
#include <math.h>
#include <iostream>
#include <sstream>
#include <iomanip>


// #include <Eigen/Sparse>
#include <mkl.h>

// MKL_NUM_THREADS=1;
/**
    @purpose    Function will iteratively compute the pressure and pressure
                gradient force that exists to force the velocity solution to
                conserve mass. The pressure is found by solving poission's
                equation using the SIMPLE algorithm. The pressure gradient
                is then applied to the velocity solution to update it. Three
                options are available for the pressure solver: 1. simple
                numerical algorithm, 2. spherical harmonic solver on the
                geodesic grid, and 3. spherical harmonic solver on an
                interpolated lat-lon grid.

    @params

    globals     pointer to the Globals object in ODIS
    mesh        pointer to the Mesh object in ODIS
    p           reference to the initialised pressure field on the geodesic grid
    v           reference to the initialised velocity field on the geodesic grid

    @returns

    p           pressure array is updated based on mass conservation
    v           velocity array is also updated such that the flow is divergence-free
    int         if convergence is achieved, returns 0

    @author     Hamish Hay

    @history    01/07/2017  - created
                15/09/2017  - spherical harmonic lat-lon pressure solver added
                03/05/2018  - added pressure solution with linear equation solver

*/
int updatePressure(Globals * globals,
                   Mesh * grid,
                   Array1D<double> & p,
                   Array2D<double> & v,
                   Array2D<double> & dvdt)
{
    int node_num, N_ll;             // node number on geodesic and lat-lon grid
    int i, j, l, m;                 // counters: nodes (i,j), SH deg/order (l, m)
    int iter, max_iter;             // no. of and max iterations
    int l_max;                      // maximum SH degree
    double r, rr, factor;           // radius, radius^2, and laplacian factor
    double epsilon;                 // divergence tolerance
    double total_div;               // sum of divergence field
    double dt;                      // time-step
    double u_factor, v_factor;      // 1/r factors for spherical coord gradients
    double u_corr, v_corr;          // east-west and north-south velocity corrections
    double ll_num;                  // double for SH degree factor

    node_num =  globals->node_num;
    N_ll =      (int)globals->dLat.Value();
    l_max =     globals->l_max.Value();
    r =         globals->radius.Value();
    rr =        pow(r, 2.0);

    iter =      0;
    max_iter =  1;
    epsilon =   1e-8;
    dt =        globals->timeStep.Value();

    u_factor =  1.0/r;
    v_factor =  1.0/r;


    //-------------- Declare and allocate new arrays to be stored --------------
    Array1D<double> * v_div;            // divergence of the velocity field
    Array1D<double> * p_corr;           // pressure correction field
    Array2D<double> * ll_v_div;         // divergence of the vel field on a lat-lon grid
    Array2D<double> * v_temp;

    Array3D<double> * p_lm;             // SH coeffs of the pressure field
    Array3D<double> * div_lm;           // SH coeffs of the divergence field

    v_div =     new Array1D<double>(node_num);
    v_temp =    new Array2D<double>(node_num,2);
    ll_v_div =  new Array2D<double>(180/N_ll, 360/N_ll);
    p_corr =    new Array1D<double>(node_num);
    div_lm =    new Array3D<double>(2.*(l_max+1), 2.*(l_max+1), 2);
    p_lm =      new Array3D<double>(2.*(l_max+1), 2.*(l_max+1), 2);


    //-------------- Declare and point arrays to existing grid arrays ----------
    Array3D<double> * Pbar_lm;          // Associated legendre functions at every node
    Array3D<double> * Pbar_lm_deriv;    // Associated legendre function derivatives at every node
    Array3D<double> * trigMLon;         // cos and sin of m*longitude at every node
    Array2D<double> * trigLat;          // cos and sin of latitude at every node

    Pbar_lm =       &(grid->Pbar_lm);
    Pbar_lm_deriv = &(grid->Pbar_lm_deriv);
    trigMLon =      &(grid->trigMLon);
    trigLat =       &(grid->trigLat);

    //------------------------ Elliptical solver set up ------------------------
    _MKL_DSS_HANDLE_t * solverHandle;   // handle to the linear equation solver
    solverHandle = &(grid->pressureSolverHandle); // point to solver set up in mesh.cpp


    // initialise key arrays to zero
    for (i=0; i<node_num; i++)
    {
        (*v_div)(i) = 0.0;
        (*p_corr)(i) = 0.0;
        p(i) = 0.0;
    }

    double * hVar = new double[node_num];
    double a_ratio = 0.0;
    for (i=0; i<node_num; i++) hVar[i] = globals->h.Value()*(1. + a_ratio*(*trigLat)(i,0)*pow((*trigMLon)(i,4,0),4.0));

    for (i=0; i<node_num; i++)
    {
        (*v_temp)(i,0) = v(i,0) * hVar[i];
        (*v_temp)(i,1) = v(i,1) * hVar[i];
    }


    //------------------------- BEGIN PRESSURE SOLVER --------------------------
    // Below, ODIS attempts to correct the current velocity solution v to obey
    // the continuity constraint (either div(v) = 0 or div(h*v) = 0). A pressure
    // correction is found by solving an elliptical equation, either using
    // Laplacian spherical harmonic identities or a direct numerical solution
    // to A*p = d, where A is a stored coefficient matrix, p is the pressure
    // correction vector, and d is the velocity field divergence vector.

    int method = 1;
    switch (method)
    {
        case 1:
            {
                int nRhs = 1;       // number of rhs vectors in the

                double * p_temp;    // pointer to pressure correction array
                double * v_div_p;   // pointer to div(v) array

                p_temp =    &((*p_corr)(0));
                v_div_p =   &((*v_div)(0));

                MKL_INT opt = MKL_DSS_REFINEMENT_ON;    // intel sovler options

                total_div = 0.0;

                // velocityDivergence(grid, *v_div, v, total_div, -1.0);
                velocityDivergence(grid, *v_div, *v_temp, total_div, -1.0);
                // velocityDivergence(grid, *v_div, v, total_div, -globals->h.Value());

                iter = 0;
                do
                {


                    // Solve A*p = d to find the pressure correction
                    dss_solve_real(*solverHandle, opt, v_div_p, nRhs, p_temp);

                    for (i=0; i<node_num; i++)
                    {
                        p_temp[i] *= 1.0/dt;            // remember, p_temp points to p_corr
                        p(i) += (*p_corr)(i)*1000.0;    // multiply by density
                        (*v_div)(i) = 0.0;              // reset divergence array
                    }

                    // Apply pressure correction gradient to force div(v) = 0
                    pressureGradient(grid, v, *p_corr,  node_num, dt);
                    // pressureGradient(grid, dvdt, *p_corr,  node_num, 1.0);

                    // avgAtPoles(grid, v);

                    for (i=0; i<node_num; i++)
                    {
                        (*v_temp)(i,0) = v(i,0) * hVar[i];
                        (*v_temp)(i,1) = v(i,1) * hVar[i];
                    }
                    // Re-evalute div(v)
                    std::cout<<iter<<'\t'<<total_div<<std::endl;
                    total_div = 0.0;

                    // velocityDivergence(grid, *v_div, v, total_div, -1.0);
                    velocityDivergence(grid, *v_div, *v_temp, total_div, -1.0);
                    // velocityDivergence(grid, *v_div, v, total_div, -globals->h.Value());


                    iter++;
                    std::cout<<iter<<'\t'<<total_div<<std::endl;
                }
                while ((iter < max_iter)  && (total_div > epsilon));
            }
            break;

        case 2:
            {
                total_div = 0.0;
                velocityDivergence(grid, *v_div, v, total_div, -1.0);

                iter = 0;

                do
                {
                    // FIND AND STORE THE SPHERICAL HARMONIC COEFFICIENTS OF THE
                    // PRESSURE FIELD
                    // interpolateGG2LL(globals,
                    //                  grid,
                    //                  *ll_v_div,
                    //                  *v_div,
                    //                  grid->ll_map_coords,
                    //                  grid->V_inv,
                    //                  grid->cell_ID);

                    interpolateGG2LLConservative(globals,
                                             grid,
                                             *ll_v_div,
                                             *v_div);

                    getSHCoeffsLL(*ll_v_div, *div_lm, N_ll, l_max);

                    for (l=1; l<l_max+1; l++)
                    {
                        ll_num = (double)(l*(l+1));
                        factor = -rr / (dt*ll_num);
                        for (m=0; m<=l; m++)
                        {
                            (*p_lm)(l, m, 0) = factor * (*div_lm)(l, m, 0);  // C_lm
                            (*p_lm)(l, m, 1) = factor * (*div_lm)(l, m, 1);  // S_lm
                        }
                    }

                    double u_corr, v_corr;
                    u_corr = 0.0;
                    v_corr = 0.0;

                    // RECONSTRUCT PRESSURE CORRECTION FIELD USING SH COEFFICIENTS
                    for (i=0; i<node_num; i++)
                    {
                        u_factor     = 1.0/(r*(*trigLat)(i,0));
                        (*p_corr)(i) = 0.0;
                        u_corr       = 0.0;
                        v_corr       = 0.0;
                        for (l=1; l<l_max+1; l++)
                        {
                            for (m=0; m<=l; m++)
                            {
                                (*p_corr)(i) += (*Pbar_lm)(i, l, m) * ( (*p_lm)(l, m, 0) * (*trigMLon)(i, m, 0) +
                                                                (*p_lm)(l, m, 1) * (*trigMLon)(i, m, 1));

                                // u_corr += u_factor * (*Pbar_lm)(i, l, m)
                                //             * (-(*p_lm)(l, m, 0) * (double)m * (*trigMLon)(i, m, 1)
                                //                + (*p_lm)(l, m, 1) * (double)m * (*trigMLon)(i, m, 0));
                                //
                                // v_corr += v_factor * (*Pbar_lm_deriv)(i, l, m)
                                //            * ((*p_lm)(l, m, 0) * (*trigMLon)(i, m, 0)
                                //               + (*p_lm)(l, m, 1) * (*trigMLon)(i, m, 1));
                            }
                        }

                        // v(i,0) -= dt*u_corr;
                        // v(i,1) -= dt*v_corr;

                        p(i) += (*p_corr)(i)*1000.0;
                    }

                    // Apply pressure correction gradient to force div(v) = 0
                    pressureGradient(grid, v, *p_corr,  node_num, dt);
                    // pressureGradient(grid, dvdt, *p_corr,  node_num, 1.0);

                    // avgAtPoles(grid, v);

                    total_div = 0.0;
                    velocityDivergence(grid, *v_div, v, total_div, -1.0);

                    // std::cout<<total_div<<std::endl;

                    iter++;

                }
                while ((iter < max_iter)  && (total_div > epsilon));

            }
            break;
    }

    delete v_div;
    delete v_temp;
    delete ll_v_div;
    delete p_corr;

    delete div_lm;
    delete p_lm;

    delete hVar;

    if (total_div > epsilon) {
        // std::cout<<"PRESSURE FIELD DID NOT CONVERGE AFTER ITER="<<iter<<"! "<<std::endl;
        return 1;
    }
    else {
        // std::cout<<"PRESSURE FIELD CONVERGED AT ITER="<<iter<<"! "<<std::endl;
        return 0;

    }

};
