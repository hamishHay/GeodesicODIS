#include "mesh.h"
#include "globals.h"
#include "outFiles.h"
#include "solver.h"
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include "temporalOperators.h"
#include <math.h>
#include <iostream>
#include <sstream>
#include <iomanip>

int integrateAB3vector(Globals * globals,
                 Mesh * grid,
                 Array2D<double> & v_t0,        // current solution
                 Array2D<double> & v_tm1,       // solution at previous time level
                 Array2D<double> & dv_dt_t0,    // velocity time grad at t=0
                 Array2D<double> & dv_dt_tm1,   // velocity time grad at t=-1
                 Array2D<double> & dv_dt_tm2,   // velocity time grad at t=-2
                 int iter)                      // simulation interation num
 {

     double a, b, c;        // AB3 integration constants
     double dt;
     int node_num, i;

     dt = globals->timeStep.Value();
     node_num = globals->node_num;

     a = 23./12.;
     b = -16./12.;
     c = 5./12.;

     if (iter > 1) {

         for (i = 0; i<node_num; i++)
         {
             v_t0(i,0) = (a*dv_dt_t0(i,0)
                           + b*dv_dt_tm1(i,0)
                           + c*dv_dt_tm2(i,0)) * dt
                           + v_tm1(i,0);

             v_t0(i,1) = (a*dv_dt_t0(i,1)
                           + b*dv_dt_tm1(i,1)
                           + c*dv_dt_tm2(i,1)) * dt
                           + v_tm1(i,1);

             v_tm1(i, 0) = v_t0(i, 0);
             v_tm1(i, 1) = v_t0(i, 1);

             dv_dt_tm2(i, 0) = dv_dt_tm1(i, 0);
             dv_dt_tm2(i, 1) = dv_dt_tm1(i, 1);

             dv_dt_tm1(i, 0) = dv_dt_t0(i, 0);
             dv_dt_tm1(i, 1) = dv_dt_t0(i, 1);
         }
     }
     else if (iter == 0)
     {
         for (i = 0; i<node_num; i++)
         {
             v_t0(i,0) = dv_dt_t0(i,0) * dt + v_tm1(i,0);
             v_t0(i,1) = dv_dt_t0(i,1) * dt + v_tm1(i,1);

             v_tm1(i, 0) = v_t0(i, 0);
             v_tm1(i, 1) = v_t0(i, 1);

             dv_dt_tm2(i, 0) = dv_dt_t0(i, 0);
             dv_dt_tm2(i, 1) = dv_dt_t0(i, 1);
         }
     }
     else if (iter == 1)
     {
         for (i = 0; i<node_num; i++)
         {
             v_t0(i,0) = dv_dt_t0(i,0) * dt + v_tm1(i,0);
             v_t0(i,1) = dv_dt_t0(i,1) * dt + v_tm1(i,1);

             v_tm1(i, 0) = v_t0(i, 0);
             v_tm1(i, 1) = v_t0(i, 1);

             dv_dt_tm1(i, 0) = dv_dt_t0(i, 0);
             dv_dt_tm1(i, 1) = dv_dt_t0(i, 1);
         }

     }

     return 1;
 };

 int integrateAB3scalar(Globals * globals,
                  Mesh * grid,
                  Array1D<double> & s_t0,        // current solution
                  Array1D<double> & s_tm1,       // solution at previous time level
                  Array1D<double> & ds_dt_t0,    // velocity time grad at t=0
                  Array1D<double> & ds_dt_tm1,   // velocity time grad at t=-1
                  Array1D<double> & ds_dt_tm2,   // velocity time grad at t=-2
                  int iter)                      // simulation interation num
  {

      double a, b, c;        // AB3 integration constants
      double dt;
      int node_num, i;

      dt = globals->timeStep.Value();
      node_num = globals->node_num;

      a = 23./12.;
      b = -16./12.;
      c = 5./12.;

      if (iter > 1) {
          for (i = 0; i<node_num; i++)
          {
              s_t0(i) = (a*ds_dt_t0(i)
                        + b*ds_dt_tm1(i)
                        + c*ds_dt_tm2(i)) * dt
                        + s_tm1(i);

              s_tm1(i) = s_t0(i);

              ds_dt_tm2(i) = ds_dt_tm1(i);
              ds_dt_tm1(i) = ds_dt_t0(i);
          }
      }
      else if (iter == 0)
      {
          for (i = 0; i<node_num; i++)
          {
              s_t0(i) = ds_dt_t0(i) * dt + s_tm1(i);

              s_tm1(i) = s_t0(i);
              ds_dt_tm2(i) = ds_dt_t0(i);
          }
      }
      else if (iter == 1)
      {
          for (i = 0; i<node_num; i++)
          {
              s_t0(i) = ds_dt_t0(i) * dt + s_tm1(i);

              s_tm1(i) = s_t0(i);
              ds_dt_tm1(i) = ds_dt_t0(i);
          }
      }

      return 1;
  };
