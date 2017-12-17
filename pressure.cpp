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


/**
*   @purpose    Function will iteratively compute the pressure and pressure
*               gradient force that exists to force the velocity solution to
*               conserve mass. The pressure is found by solving poission's
*               equation using the SIMPLE algorithm. The pressure gradient
*               is then applied to the velocity solution to update it. Three
*               options are available for the pressure solver: 2. simple
*               numerical algorithm, 1. spherical harmonic solver on the
*               geodesic grid, and 0. spherical harmonic solver on an
*               interpolated lat-lon grid.
*
*   @params
*
*   globals     pointer to the Globals object in ODIS
*   mesh        pointer to the Mesh object in ODIS
*   p           reference to the initialised pressure field on the geodesic grid
*   v           reference to the initialised velocity field on the geodesic grid
*
*   @returns
*
*   p           pressure array is updated based on mass conservation
*   v           velocity array is also updated such that the flow is divergence-free
*   int         if convergence is achieved, returns 0
*
*   @author     Hamish Hay
*
*   @history    01/07/2017  - created
*               15/09/2017  - spherical harmonic lat-lon pressure solver added
*               03/11/2017  - grad P applied directly in harmonics
*
*/
int updatePressure(Globals * globals,
                   Mesh * grid,
                   Array1D<double> & p,
                   Array2D<double> & v)
{
    int node_num, i, j, l, m, iter, max_iter, l_max, N_ll;
    double h, r, rr, u_factor, v_factor, factor, epsilon, total_div_new, total_div_old, relax_f, dt, p_min;
    double residual;
    double div_diff;

    int method = 0;

    double ll_num;
    double u_corr, v_corr;

    Array1D<double> * v_div;
    Array2D<double> * v_temp;
    Array1D<double> * p_corr;
    Array1D<double> * p_dummy;
    Array1D<double> * p_factor;
    Array2D<double> * coords;

    Array3D<double> * p_lm;         // SH coeffs of the divergence field
    Array3D<double> * div_lm;       // SH coeffs of the pressure field

    Array3D<double> * Pbar_lm;
    Array3D<double> * Pbar_lm_deriv;
    Array3D<double> * trigMLon;
    Array2D<double> * trigLat;

    Array2D<double> * ll_v_div;

    node_num = globals->node_num;
    l_max = globals->l_max.Value();
    h = globals->h.Value();
    r = globals->radius.Value();
    rr = pow(r, 2.0);

    u_factor = 1.0/r;
    v_factor = 1.0/r;

    N_ll = (int)globals->dLat.Value();

    iter = 0;
    max_iter = 5;
    epsilon = 1e-20;
    relax_f = 1.0;
    dt = globals->timeStep.Value();

    v_div = new Array1D<double>(node_num);
    ll_v_div = new Array2D<double>(180/N_ll, 360/N_ll);
    p_corr = new Array1D<double>(node_num);
    p_dummy = new Array1D<double>(node_num);

    div_lm = new Array3D<double>(2.*(l_max+1), 2.*(l_max+1), 2);
    p_lm = new Array3D<double>(2.*(l_max+1), 2.*(l_max+1), 2);

    Pbar_lm = &(grid->Pbar_lm);             // 4-pi Normalised associated leg funcs
    Pbar_lm_deriv = &(grid->Pbar_lm_deriv); // 4-pi Normalised associated leg funcs derivs
    trigMLon = &(grid->trigMLon);           // cos and sin of m*longitude
    trigLat = &(grid->trigLat);

    p_factor = &(grid->pressure_factor);

    for (i=0; i<node_num; i++)
    {
        (*v_div)(i) = 0.0;
        (*p_corr)(i) = 0.0;
        (*p_dummy)(i) = 0.0;
        p(i) = 0.0;
    }

    total_div_old = 0.0;
    switch (method) {
        /* METHOD 0: Expand divergence in harmonics on a lat-lon grid to find
        *            and apply the internal pressure gradient
        */
        case 0:
            for (i=0; i<node_num; i++)
            {
                (*v_div)(i) = 0.0;
            }

            total_div_new = 0.0;
            velocityDivergence(grid, *v_div, v, total_div_new, -1000.0);

            div_diff = total_div_new;
            total_div_old = 0.0;

            iter = 0;

            do
            {
                // FIND AND STORE THE SPHERICAL HARMONIC COEFFICIENTS OF THE
                // PRESSURE FIELD
                interpolateGG2LL(globals,
                                 grid,
                                 *ll_v_div,
                                 *v_div,
                                 grid->ll_map_coords,
                                 grid->V_inv,
                                 grid->cell_ID);

                getSHCoeffsLL(*ll_v_div, *div_lm, N_ll, l_max);

                // SOLVE POISSON EQUATION FOR PRESSURE
                for (l=1; l<l_max+1; l++)
                {
                    ll_num = (double)(l*(l+1));
                    factor = -rr * relax_f * 1.0 / (ll_num);
                    for (m=0; m<=l; m++)
                    {
                        (*p_lm)(l, m, 0) = factor * (*div_lm)(l, m, 0);  // C_lm    UNITS: Pressure * time
                        (*p_lm)(l, m, 1) = factor * (*div_lm)(l, m, 1);  // S_lm
                    }
                }


                u_corr = 0.0;
                v_corr = 0.0;

                // RECONSTRUCT PRESSURE CORRECTION FIELD USING SH COEFFICIENTS
                for (i=0; i<node_num; i++)
                {

                    u_factor = 1.0/(r*(*trigLat)(i,0));
                    (*p_corr)(i) = 0.0;
                    u_corr = 0.0;
                    v_corr = 0.0;
                    for (l=1; l<l_max+1; l++)
                    {
                        for (m=0; m<=l; m++)
                        {
                            (*p_corr)(i) += (*Pbar_lm)(i, l, m) * ( (*p_lm)(l, m, 0) * (*trigMLon)(i, m, 0) +
                                                            (*p_lm)(l, m, 1) * (*trigMLon)(i, m, 1));

                            u_corr += u_factor * (*Pbar_lm)(i, l, m)
                                        * (-(*p_lm)(l, m, 0) * (double)m * (*trigMLon)(i, m, 1)
                                           + (*p_lm)(l, m, 1) * (double)m * (*trigMLon)(i, m, 0));

                            v_corr += v_factor * (*Pbar_lm_deriv)(i, l, m)
                                       * ((*p_lm)(l, m, 0) * (*trigMLon)(i, m, 0)
                                          + (*p_lm)(l, m, 1) * (*trigMLon)(i, m, 1));

                        }
                    }

                    // if (i>1)
                    // {
                    //     v(i,0) = v(i,0) - u_corr/1000.0;
                    //     v(i,1) = v(i,1) - v_corr/1000.0;
                    //
                    // }
                    // else
                    // {
                    //     v(i,0) = 0.0;
                    //     v(i,1) = 0.0;
                    // }

                    v(i,0) = v(i,0) - u_corr/1000.0;
                    v(i,1) = v(i,1) - v_corr/1000.0;

                    p(i) += (*p_corr)(i)/dt;

                    (*v_div)(i) = 0.0;

                }

                // pressureGradient(grid, v, *p_corr, 1.0/1000.0);

                div_diff = total_div_new;
                total_div_new = 0.0;
                velocityDivergence(grid, *v_div, v, total_div_new, -1000.0);

                residual = fabs( total_div_new - total_div_old);
                total_div_old = total_div_new;

                iter++;
            }
            while ((iter < max_iter));


            break;

        /* METHOD 1: Expand divergence in harmonics on the geodesic grid to find
        *            and apply the internal pressure gradient
        */
        case 1:
            for (i=0; i<node_num; i++) (*v_div)(i) = 0.0;

            total_div_new = 0.0;
            velocityDivergence(grid, *v_div, v, total_div_new, -1.0);

            div_diff = total_div_new;

            iter = 0;

            do
            {
                // FIND AND STORE THE SPHERICAL HARMONIC COEFFICIENTS OF THE
                // PRESSURE FIELD
                getSHCoeffsGG(grid->node_pos_sph, *v_div, *div_lm, node_num, l_max);

                for (l=1; l<l_max+1; l++)
                {
                    factor = -rr / (double)( l * (l + 1));
                    for (m=0; m<=l; m++)
                    {
                        (*p_lm)(l, m, 0) = factor * (*div_lm)(l, m, 0);  // C_lm
                        (*p_lm)(l, m, 1) = factor * (*div_lm)(l, m, 1);  // S_lm

                    }
                }

                // RECONSTRUCT PRESSURE CORRECTION FIELD USING SH COEFFICIENTS
                for (i=0; i<node_num; i++)
                {
                    u_factor = 1.0/(r*(*trigLat)(i,0));
                    (*p_corr)(i) = 0.0;

                    for (l=1; l<l_max+1; l++)
                    {
                        for (m=0; m<=l; m++)
                        {
                            (*p_corr)(i) += (*Pbar_lm)(i, l, m) * ( (*p_lm)(l, m, 0) * (*trigMLon)(i, m, 0) +
                                                            (*p_lm)(l, m, 1) * (*trigMLon)(i, m, 1));

                            v(i, 0) -= u_factor * (*Pbar_lm)(i, l, m)
                                        * (-(*p_lm)(l, m, 0) * (double)m * (*trigMLon)(i, m, 1)
                                           + (*p_lm)(l, m, 1) * (double)m * (*trigMLon)(i, m, 0));

                            v(i, 1) -= v_factor * (*Pbar_lm_deriv)(i, l, m)
                                       * ((*p_lm)(l, m, 0) * (*trigMLon)(i, m, 0)
                                          + (*p_lm)(l, m, 1) * (*trigMLon)(i, m, 1));
                        }
                    }

                    p(i) += 1000.0 * (*p_corr)(i)/dt;
                    (*v_div)(i) = 0.0;
                }

                // CHECK SUM OF OCEAN DIVERGENCE
                div_diff = total_div_new;
                total_div_new = 0.0;
                velocityDivergence(grid, *v_div, v, total_div_new, -1.0);
                div_diff = fabs(total_div_new - div_diff);


                residual = fabs( total_div_new - total_div_old);
                total_div_old = total_div_new;

                iter++;
            }
            while ((iter < max_iter)  && (residual > epsilon));

            break;

        /* METHOD 2: Find pressure field using an iterative and approximate
        *            finite volume solver for the pressre Poisson equation.
        *            This method is extremely slow and requires many iterations
        *            before convergence is reached.
        */
        case 2:

            // The below commented section uses a simple iterative approach to pressure
            // solving. It is very slow.

            total_div_new = 0.0;
            velocityDivergence(grid, *v_div, v, total_div_new, -1.0);


            max_iter = 10000;
            iter = 0;
            do {
                p_min = 0.0;

                for (i=0; i<node_num; i++)
                {
                    (*p_corr)(i) = (*p_factor)(i) * (*v_div)(i);
                    p(i) += -1000.0 * (*p_corr)(i);
                }

                pressureGradient(grid, v, *p_corr, -dt);

                iter++;

                total_div_new = 0.0;
                for (i=0; i<node_num; i++) (*v_div)(i) = 0.0;
                velocityDivergence(grid, *v_div, v, total_div_new, -1.0);

                residual = fabs( total_div_new - total_div_old);
                total_div_old = total_div_new;
            }
            while ((iter < max_iter)  && (residual > epsilon));

            break;
    }

    delete v_div;
    delete ll_v_div;
    delete p_corr;

    delete div_lm;
    delete p_lm;

    if (total_div_new > epsilon && residual > epsilon)
    {
        return 1;
    }
    else
    {
        return 0;
    }

};
