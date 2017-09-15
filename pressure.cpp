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



int updatePressure(Globals * globals, Mesh * grid, Array1D<double> & p, Array2D<double> & v)
{
    int node_num, i, j, l, m, iter, max_iter, l_max, N_ll;
    double h, r, rr, u_factor, v_factor, factor, epsilon, total_div_new, total_div_old, relax_f, dt, p_min;
    double residual;
    double div_diff;

    Array1D<double> * v_div;
    Array2D<double> * v_temp;
    Array1D<double> * p_corr;
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

    max_iter = 50;
    epsilon = 1e-3;//((double)node_num)/1e9;
    relax_f = 1.0;
    dt = globals->timeStep.Value();

    v_div = new Array1D<double>(node_num);
    ll_v_div = new Array2D<double>(180/N_ll, 360/N_ll);
    p_corr = new Array1D<double>(node_num);

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
        p(i) = 0.0; // PERHAPS THE PRESSURE FIELD SHOULD ONLY DEPEND ON THE CURRENT
                    // TIMESTEPS PRESSURE CORRECTION?
    }

    //------------------------- BEGIN PRESSURE SOLVER --------------------------

    // FIND THE DIVERGENCE FIELD

    int method = 1;

    switch (method) {

        case 0:
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
                    factor = -rr / (double)( l * (l + 1) );
                    for (m=0; m<=l; m++)
                    {
                        (*p_lm)(l, m, 0) = factor * (*div_lm)(l, m, 0);  // C_lm
                        (*p_lm)(l, m, 1) = factor * (*div_lm)(l, m, 1);  // S_lm
                    }
                }

                // RECONSTRUCT PRESSURE CORRECTION FIELD USING SH COEFFICIENTS
                for (i=0; i<node_num; i++)
                {
                    (*p_corr)(i) = 0.0;
                    for (l=1; l<l_max+1; l++)
                    {
                        for (m=0; m<=l; m++)
                        {
                            (*p_corr)(i) += (*Pbar_lm)(i, l, m) * ( (*p_lm)(l, m, 0) * (*trigMLon)(i, m, 0) +
                                                            (*p_lm)(l, m, 1) * (*trigMLon)(i, m, 1));

                        }
                    }

                    p(i) += (*p_corr)(i);
                    (*v_div)(i) = 0.0;
                }

                // APPLY PRESSURE GRADIENT FORCE
                pressureGradient(grid, v, *p_corr, 1.0);

                // CHECK SUM OF OCEAN DIVERGENCE
                div_diff = total_div_new;
                total_div_new = 0.0;
                velocityDivergence(grid, *v_div, v, total_div_new, -1.0);
                div_diff = fabs(total_div_new - div_diff);

                iter++;
            }
            while ((total_div_new > epsilon) && (iter < max_iter)  && (div_diff > 1e-5));

            break;

        case 1:
            for (i=0; i<node_num; i++) (*v_div)(i) = 0.0;

            total_div_new = 0.0;
            velocityDivergence(grid, *v_div, v, total_div_new, -1.0);

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

                for (l=1; l<l_max+1; l++)
                {
                    factor = -rr * 1.0 / (double)( l * (l + 1) );
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
                    // v(i,0) = 0.0;
                    // v(i,1) = 0.0;
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

                    p(i) += 1000.0*(*p_corr)(i)/dt;
                    (*v_div)(i) = 0.0;
                }

                // }
                // CHECK SUM OF OCEAN DIVERGENCE
                div_diff = total_div_new;
                total_div_new = 0.0;
                velocityDivergence(grid, *v_div, v, total_div_new, -1.0);
                // div_diff = fabs(total_div_new - div_diff);

                residual = fabs( total_div_new - total_div_old);
                total_div_old = total_div_new;

                iter++;
            }
            while ((iter < max_iter)  && (residual > 1e-6));

            break;
        case 2:
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
                    factor = -rr / (double)( l * (l + 1) );
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
                    // v(i,0) = 0.0;
                    // v(i,1) = 0.0;
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

                    p(i) += (*p_corr)(i);
                    (*v_div)(i) = 0.0;
                }

                // CHECK SUM OF OCEAN DIVERGENCE
                div_diff = total_div_new;
                total_div_new = 0.0;
                velocityDivergence(grid, *v_div, v, total_div_new, -1.0);
                div_diff = fabs(total_div_new - div_diff);


                residual = fabs( total_div_new - total_div_old);
                total_div_old = total_div_new;
                // std::cout<<"GLOBAL DIVERGENCE: "<<total_div_new<<' '<<div_diff<<std::endl;

                iter++;
            }
            while ((total_div_new > epsilon) && (iter < max_iter)  && (div_diff > 1e-5));

            break;
        case 3:

            // The below commented section uses a simple iterative approach to pressure
            // solving. It is very slow.
            //
            max_iter = 1200;
            iter = 0;
            do {
                total_div_new = 0.0;
                p_min = 0.0;

                velocityDivergence(grid, *v_div, v, total_div_new, -1.0);

                for (i=0; i<node_num; i++)
                {
                    (*p_corr)(i) = (*p_factor)(i) * (*v_div)(i);
                    p(i) += -(*p_corr)(i);
                }

                pressureGradient(grid, v, *p_corr, -dt);

                iter++;
            }
            while ((total_div_new > epsilon) && (iter < max_iter));

            total_div_new = 0.0;
            velocityDivergence(grid, *v_div, v, total_div_new, -1.0);

            std::cout<<"GLOBAL DIVERGENCE: "<<total_div_new<<std::endl;

            break;
    }

    if (total_div_new > epsilon && residual > 1e-6) return 1;
    else {
        std::cout<<"PRESSURE FIELD CONVERGED AT ITER="<<iter<<"! "<<std::endl;
        return 0;

    }

    delete v_div;
    delete ll_v_div;
    delete p_corr;

    delete div_lm;
    delete p_lm;
}
