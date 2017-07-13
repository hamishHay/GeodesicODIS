#include "mesh.h"
#include "globals.h"
#include "outFiles.h"
#include "solver.h"
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include "pressure.h"
#include "spatialOperators.h"
#include <math.h>
#include <iostream>
#include <sstream>
#include <iomanip>

int updatePressure(Globals * globals, Mesh * grid, Array1D<double> & p, Array2D<double> & v)
{
    int node_num, i, iter, max_iter;
    double h, epsilon, total_div, relax_f, dt;

    Array1D<double> * v_div;
    Array1D<double> * p_corr;
    Array1D<double> * p_factor;

    node_num = globals->node_num;
    h = globals->h.Value();
    max_iter = 1000;
    epsilon = 5e-4;
    relax_f = 1.0;
    dt = globals->timeStep.Value();

    v_div = new Array1D<double>(node_num);
    p_corr = new Array1D<double>(node_num);

    p_factor = &(grid->pressure_factor);

    for (i=0; i<node_num; i++)
    {
        (*v_div)(i) = 0.0;
        // (*v_div)(i,1) = 0.0;
        (*p_corr)(i) = 0.0;
    }

    //------------------------- BEGIN PRESSURE SOLVER --------------------------

    iter = 0;
    do {
        total_div = 0.0;

        velocityDivergence(grid, *v_div, v, total_div, -1.0);
        // velocityDivergence(grid, p, v, total_div, -1.0);

        for (i=0; i<node_num; i++)
        {
            (*p_corr)(i) = (*p_factor)(i) * (*v_div)(i);
            p(i) += -(*p_corr)(i);
            p(i) = std::max(0.0, p(i));
        }

        pressureGradient(grid, v, *p_corr, -dt/1000.0);


        // std::cout<<total_div<<std::endl;
        iter++;
    }
    while ((total_div > epsilon) && (iter < max_iter));

    if (iter >= max_iter) return 1;

    // CHECK VELOCITY FIELD DIVERGENCE

    // ---- CALCULATE PRESSURE CORRECTION

    // ---- UPDATE PRESSURE FIELD

    // ---- UPDATE VELOCITY FIELD

    // ---- REPEAT UNTIL DIVERGENCE IS ZERO

    delete v_div;
    delete p_corr;

    return 0;

}
