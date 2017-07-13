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
    double h, epsilon, total_div, relax_f;

    Array1D<double> * v_div;
    Array1D<double> * p_corr;

    node_num = globals->node_num;
    h = globals->h.Value();
    max_iter = 1;
    epsilon = 1e-4;
    relax_f = 1.0;

    v_div = new Array1D<double>(node_num);
    p_corr = new Array1D<double>(node_num);

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

        velocityDivergence(grid, *v_div, v, total_div);



        std::cout<<total_div<<std::endl;
        iter++;
    }
    while ((total_div > epsilon) && (iter < max_iter));

    if (iter >= max_iter) return 0;

    // CHECK VELOCITY FIELD DIVERGENCE

    // ---- CALCULATE PRESSURE CORRECTION

    // ---- UPDATE PRESSURE FIELD

    // ---- UPDATE VELOCITY FIELD

    // ---- REPEAT UNTIL DIVERGENCE IS ZERO

    delete v_div;
    delete p_corr;

    return 1;

}
