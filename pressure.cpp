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
#include <math.h>
#include <iostream>
#include <sstream>
#include <iomanip>



int updatePressure(Globals * globals, Mesh * grid, Array1D<double> & p, Array2D<double> & v)
{
    int node_num, i, iter, max_iter, l_max;
    double h, epsilon, total_div, relax_f, dt, p_min;

    Array1D<double> * v_div;
    Array1D<double> * p_corr;
    Array1D<double> * p_factor;

    Array2D<double> * p_lm;         // SH coeffs of the divergence field
    Array2D<double> * div_lm;       // SH coeffs of the pressure field

    node_num = globals->node_num;
    h = globals->h.Value();
    max_iter = 2000;
    epsilon = ((double)node_num)/1e9;
    relax_f = 1.0;
    dt = globals->timeStep.Value();

    v_div = new Array1D<double>(node_num);
    p_corr = new Array1D<double>(node_num);

    div_lm = new Array2D<double>(l_max+1, l_max+1);
    p_lm = new Array2D<double>(l_max+1, l_max+1);

    p_factor = &(grid->pressure_factor);

    for (i=0; i<node_num; i++)
    {
        (*v_div)(i) = 0.0;
        // (*v_div)(i,1) = 0.0;
        (*p_corr)(i) = 0.0;
        p(i) = 0.0; // PERHAPS THE PRESSURE FIELD SHOULD ONLY DEPEND ON THE CURRENT
                    // TIMESTEPS PRESSURE CORRECTION?
    }

    //------------------------- BEGIN PRESSURE SOLVER --------------------------

    // FIND THE DIVERGENCE FIELD

    total_div = 0.0;
    velocityDivergence(grid, *v_div, v, total_div, -1.0);

    // FIND AND STORE THE SPHERICAL HARMONIC COEFFICIENTS OF THE PRESSURE FIELD





    // The below commented section uses a simple iterative approach to pressure
    // solving. It is very slow.

    // iter = 0;
    // do {
    //     total_div = 0.0;
    //     p_min = 0.0;
    //
    //     velocityDivergence(grid, *v_div, v, total_div, -1.0);
    //     // velocityDivergence(grid, p, v, total_div, -1.0);
    //
    //     for (i=0; i<node_num; i++)
    //     {
    //         (*p_corr)(i) = (*p_factor)(i) * (*v_div)(i);
    //         p(i) += -(*p_corr)(i);
    //         // p(i) = std::max(0.0, p(i));
    //
    //         // p_min = std::min(p_min, p(i));
    //     }
    //
    //     pressureGradient(grid, v, *p_corr, -dt);
    //
    //     iter++;
    // }
    // while ((total_div > epsilon) && (iter < max_iter));
    //
    // // p_min = fabs(p_min);
    // // for (i=0; i<node_num; i++)
    // // {
    // //     // p_min = std::min(p_min, p(i));
    // //     p(i) += p_min;
    // // }
    //
    //
    // if (iter >= max_iter) return 1;
    // else {
    //     // std::cout<<"PRESSURE FIELD CONVERGED AT ITER="<<iter<<"! "<<std::endl;
    //     return 0;
    //
    // }

    delete v_div;
    delete p_corr;
}
