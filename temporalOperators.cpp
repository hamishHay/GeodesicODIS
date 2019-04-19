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

#include <mkl.h>
#include <omp.h>





int integrateAB3scalar(Globals * globals,
                 Mesh * grid,
                 Array1D<double> & s,        // current solution
                 Array2D<double> & ds_dt,    // velocity time grad at t=0
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
          #pragma omp parallel for
          for (i = 0; i<node_num; ++i)
          {
              s(i) += (a*ds_dt(i,0)
                            + b*ds_dt(i, 1)
                            + c*ds_dt(i, 2)) * dt;

              ds_dt(i, 2) = ds_dt(i, 1);
              ds_dt(i, 1) = ds_dt(i, 0);
          }
     }
     else if (iter == 0)
     {
         for (i = 0; i<node_num; i++)
         {
             s(i) += ds_dt(i, 0) * dt;
             ds_dt(i, 2) = ds_dt(i, 0);
             // ds_dt(i, 2) = ds_dt(i, 0);
         }
     }
     else if (iter == 1)
     {
         for (i = 0; i<node_num; i++)
         {
             s(i) += ds_dt(i, 0) * dt;
             ds_dt(i, 1) = ds_dt(i, 0);
         }

     }

     return 1;
 };


int integrateAB3vector(Globals * globals,
                 Mesh * grid,
                 Array2D<double> & v,        // current solution
                 Array3D<double> & dv_dt,    // velocity time grad at t=0
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

       // cblas_dgemm (CblasRowMajor,
       //                CblasNoTrans,
       //                CblasNoTrans,
       //                2*node_num,
       //                1,
       //                3,
       //                1.0,
       //                &(dv_dt(0,0,0)),
       //                3,
       //                coeffs,
       //                1,
       //                1.0,
       //                &(v(0,0)),
       //                1);

       // cblas_dgemv (CblasRowMajor,
       //              CblasNoTrans,
       //              2*node_num,
       //              3,
       //              1.0,
       //              &(dv_dt(0,0,0)),
       //              3,
       //              coeffs,
       //              1,
       //              1.0,
       //              &(v(0,0)),
       //              1);

          #pragma omp parallel for
          for (i = 0; i<node_num; ++i)
          {
              v(i,0) += (a*dv_dt(i,0, 0)
                            + b*dv_dt(i,0, 1)
                            + c*dv_dt(i,0, 2)) * dt;


              v(i,1) += (a*dv_dt(i, 1, 0)
                            + b*dv_dt(i, 1, 1)
                            + c*dv_dt(i, 1, 2)) * dt;

              dv_dt(i, 0, 2) = dv_dt(i, 0, 1);
              dv_dt(i, 0, 1) = dv_dt(i, 0, 0);

              dv_dt(i, 1, 2) = dv_dt(i, 1, 1);
              dv_dt(i, 1, 1) = dv_dt(i, 1, 0);

          }
     }
     else if (iter == 0)
     {
         for (i = 0; i<node_num; i++)
         {
             v(i,0) = dv_dt(i, 0, 0) * dt + v(i,0);
             v(i,1) = dv_dt(i, 1, 0) * dt + v(i,1);

             dv_dt(i, 0, 2) = dv_dt(i, 0, 0);
             dv_dt(i, 1, 2) = dv_dt(i, 1, 0);
         }
     }
     else if (iter == 1)
     {
         for (i = 0; i<node_num; i++)
         {
             v(i,0) = dv_dt(i,0, 0) * dt + v(i,0);
             v(i,1) = dv_dt(i,1, 0) * dt + v(i,1);

             dv_dt(i, 0, 1) = dv_dt(i, 0, 0);
             dv_dt(i, 1, 1) = dv_dt(i, 1, 0);
         }

     }

     return 1;
 };
