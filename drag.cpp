#include "drag.h"
#include <math.h>
#include <iostream>
#include <sstream>


void linearDrag(int & node_num, double & drag_coeff, Array2D<double> & dv_dt, Array2D<double> & v)
{
    int i;

    for (i=0; i<node_num; ++i)
    {
        dv_dt(i,0) -= v(i,0) * drag_coeff;
        dv_dt(i,1) -= v(i,1) * drag_coeff;
    }
};
