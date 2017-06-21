#include "mesh.h"
#include "globals.h"
#include "outFiles.h"
#include "solver.h"
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include "tidalPotentials.h"
#include "spatialOperators.h"
#include "drag.h"
#include <math.h>
#include <iostream>
#include <sstream>
#include <iomanip>


// Function to implement the time loop to solve mass and momentum equations over
// the geodesic grid. The time integration method is the explicit Euler method.
int eulerIntegrator(Globals * globals, Mesh * grid)
{
    std::ostringstream outstring;
    int i,j,iter,out_count;
    double ** pp, * p_t0, * p_tm1, * p_dt_t0;
    double * sinLat;

    Array2D<double> * vel_t0;      // velocity solution for current timestep
    Array2D<double> * vel_tm1;     // velocity solution at previous timestep (t minus 1)
    Array2D<double> * dvel_dt_t0;  // velocity time derivative at current timestep

    Array1D<double> * press_t0;      // displacement solution for current timestep
    Array1D<double> * press_tm1;     // displacement solution at previous timestep (t minus 1)
    Array1D<double> * dpress_dt_t0;  // displacement time derivative at current timestep

    double end_time, current_time, dt, out_frac, orbit_period, orbit_out_frac, out_time;
    double r, omega, e, obliq, drag_coeff;

    int node_num;

    OutFiles * Output;

    end_time = globals->endTime.Value();
    dt = globals->timeStep.Value();
    out_frac = globals->outputTime.Value();
    orbit_period = globals->period.Value();
    orbit_out_frac = orbit_period * out_frac;

    r = globals->radius.Value();
    omega = globals->angVel.Value();
    e = globals->e.Value();
    drag_coeff = globals->alpha.Value();

    // std::cout<<"Omega: "<<omega<<std::endl;
    // std::cout<<"End time: "<<end_time<<std::endl;

    node_num = globals->node_num;

    Output = globals->Output;
    pp = new double *[globals->out_tags.size()-1];

    std::cout<<globals->out_tags[0]<<std::endl;

    sinLat = &(grid->trigLat(0,1));

    outstring << "Defining arrays for Euler time integration..." << std::endl;

    dvel_dt_t0 = new Array2D<double>(node_num, 2);
    vel_t0 = new Array2D<double>(node_num, 2);
    vel_tm1 = new Array2D<double>(node_num, 2);

    press_t0 = new Array1D<double>(node_num);      // pressure solution for current timestep
    press_tm1 = new Array1D<double>(node_num);     // pressure solution at previous timestep (t minus 1)
    dpress_dt_t0 = new Array1D<double>(node_num);

    for (i=0; i<node_num; i++)
    {
        (*press_tm1)(i) = 0.0;//sinLat[i*2];
    }

    iter = 0;
    out_count = 1;
    current_time = 0.0;
    out_time = 0.0;
    while (current_time <= end_time)
    {

        switch (globals->tide_type)
        {
            case ECC:
                deg2Ecc(grid, *dvel_dt_t0, current_time, r, omega, e);
                break;
            case ECC_RAD:
                deg2EccRad(grid, *dvel_dt_t0, current_time, r, omega, e);
                break;
        }

        current_time += dt;
        out_time += dt;

        linearDrag(node_num, drag_coeff, *dvel_dt_t0, *vel_tm1);

        // Calculate pressure gradient
        pressureGradient(grid, *dvel_dt_t0, *press_tm1);

        // Explicit time integration
        for (i = 0; i<node_num; i++)
        {
            // Coriolis force
            // (*dvel_dt_t0)(i,0) += 2.0 * omega * sinLat[i*2] * (*vel_tm1)(i,1);
            // (*dvel_dt_t0)(i,1) -= 2.0 * omega * sinLat[i*2] * (*vel_tm1)(i,0);


            (*vel_t0)(i,0) = (*dvel_dt_t0)(i,0) * dt + (*vel_tm1)(i,0);
            (*vel_tm1)(i,0) = (*vel_t0)(i,0);

            (*vel_t0)(i,1) = (*dvel_dt_t0)(i,1) * dt + (*vel_tm1)(i,1);
            (*vel_tm1)(i,1) = (*vel_t0)(i,1);
        }


        // Calculate velocity divergence
        velocityDivergence(grid, *dpress_dt_t0, *vel_t0);

        // Update displacement solution
        for (i = 0; i<node_num; i++)
        {
            (*press_t0)(i) = (*dpress_dt_t0)(i) * dt + (*press_tm1)(i);
            (*press_tm1)(i) = (*press_t0)(i);
        }
// //
        // Check for output
        iter ++;

        if (out_time >= out_frac*orbit_period)
        {
            std::cout<<std::fixed << std::setprecision(8) <<"DUMPING DATA AT "<<current_time<<std::endl;

            out_time -= out_frac*orbit_period;
            pp[0] = &(*press_t0)(0);
            pp[1] = &(*vel_t0)(0,0);
            // pp[0] = &(*dvel_dt_t0)(0,0);

            Output->DumpData(globals, out_count, pp);
            out_count++;



            // Output->TerminateODIS();
        }
    }

    delete dvel_dt_t0;
    delete vel_t0;
    delete vel_tm1;

    delete dpress_dt_t0;
    delete press_t0;
    delete press_tm1;


    Output->Write(OUT_MESSAGE, &outstring);

    return 1;
};

// Function to implement the time loop to solve mass and momentum equations over
// the geodesic grid. The time integration method is the explicit Euler method.
// int adamsBashforth3(Globals * globals, Mesh * grid)
// {
//     std::ostringstream outstring;
//     int i,j,iter,out_count;
//     double ** pp, * p_t0, * p_tm1, * p_dt_t0;
//     double * sinLat;
//
//     Array2D<double> * vel_t0;      // velocity solution for current timestep
//     Array2D<double> * vel_tm1;     // velocity solution at previous timestep (t minus 1)
//     Array2D<double> * dvel_dt_t0;  // velocity time derivative at current timestep
//
//     Array1D<double> * press_t0;      // displacement solution for current timestep
//     Array1D<double> * press_tm1;     // displacement solution at previous timestep (t minus 1)
//     Array1D<double> * dpress_dt_t0;  // displacement time derivative at current timestep
//
//     double end_time, current_time, dt, out_frac, orbit_period, orbit_out_frac, out_time;
//     double r, omega, e, obliq, drag_coeff;
//
//     int node_num;
//
//     OutFiles * Output;
//
//     end_time = globals->endTime.Value();
//     dt = globals->timeStep.Value();
//     out_frac = globals->outputTime.Value();
//     orbit_period = globals->period.Value();
//     orbit_out_frac = orbit_period * out_frac;
//
//     r = globals->radius.Value();
//     omega = globals->angVel.Value();
//     e = globals->e.Value();
//     drag_coeff = globals->alpha.Value();
//
//     // std::cout<<"Omega: "<<omega<<std::endl;
//     // std::cout<<"End time: "<<end_time<<std::endl;
//
//     node_num = globals->node_num;
//
//     Output = globals->Output;
//     pp = new double *[globals->out_tags.size()-1];
//
//     std::cout<<globals->out_tags[0]<<std::endl;
//
//     sinLat = &(grid->trigLat(0,1));
//
//     outstring << "Defining arrays for Euler time integration..." << std::endl;
//
//     dvel_dt_t0 = new Array2D<double>(node_num, 2);
//     vel_t0 = new Array2D<double>(node_num, 2);
//     vel_tm1 = new Array2D<double>(node_num, 2);
//
//     press_t0 = new Array1D<double>(node_num);      // pressure solution for current timestep
//     press_tm1 = new Array1D<double>(node_num);     // pressure solution at previous timestep (t minus 1)
//     dpress_dt_t0 = new Array1D<double>(node_num);
//
//     for (i=0; i<node_num; i++)
//     {
//         (*press_tm1)(i) = 0.0;//sinLat[i*2];
//     }
//
//     iter = 0;
//     out_count = 1;
//     current_time = 0.0;
//     out_time = 0.0;
//     while (current_time <= end_time)
//     {
//
//         switch (globals->tide_type)
//         {
//             case ECC:
//                 deg2Ecc(grid, *dvel_dt_t0, current_time, r, omega, e);
//                 break;
//             case ECC_RAD:
//                 deg2EccRad(grid, *dvel_dt_t0, current_time, r, omega, e);
//                 break;
//         }
//
//         current_time += dt;
//         out_time += dt;
//
//         linearDrag(node_num, drag_coeff, *dvel_dt_t0, *vel_tm1);
//
//         // Calculate pressure gradient
//         pressureGradient(grid, *dvel_dt_t0, *press_tm1);
//
//         // Explicit time integration
//         for (i = 0; i<node_num; i++)
//         {
//             // Coriolis force
//             // (*dvel_dt_t0)(i,0) += 2.0 * omega * sinLat[i*2] * (*vel_tm1)(i,1);
//             // (*dvel_dt_t0)(i,1) -= 2.0 * omega * sinLat[i*2] * (*vel_tm1)(i,0);
//
//
//             (*vel_t0)(i,0) = (*dvel_dt_t0)(i,0) * dt + (*vel_tm1)(i,0);
//             (*vel_tm1)(i,0) = (*vel_t0)(i,0);
//
//             (*vel_t0)(i,1) = (*dvel_dt_t0)(i,1) * dt + (*vel_tm1)(i,1);
//             (*vel_tm1)(i,1) = (*vel_t0)(i,1);
//         }
//
//
//         // Calculate velocity divergence
//         velocityDivergence(grid, *dpress_dt_t0, *vel_t0);
//
//         // Update displacement solution
//         for (i = 0; i<node_num; i++)
//         {
//             (*press_t0)(i) = (*dpress_dt_t0)(i) * dt + (*press_tm1)(i);
//             (*press_tm1)(i) = (*press_t0)(i);
//         }
// // //
//         // Check for output
//         iter ++;
//
//         if (out_time >= out_frac*orbit_period)
//         {
//             std::cout<<std::fixed << std::setprecision(8) <<"DUMPING DATA AT "<<current_time<<std::endl;
//
//             out_time -= out_frac*orbit_period;
//             pp[0] = &(*press_t0)(0);
//             pp[1] = &(*vel_t0)(0,0);
//             // pp[0] = &(*dvel_dt_t0)(0,0);
//
//             Output->DumpData(globals, out_count, pp);
//             out_count++;
//
//
//
//             // Output->TerminateODIS();
//         }
//     }
//
//     delete dvel_dt_t0;
//     delete vel_t0;
//     delete vel_tm1;
//
//     delete dpress_dt_t0;
//     delete press_t0;
//     delete press_tm1;
//
//
//     Output->Write(OUT_MESSAGE, &outstring);
//
//     return 1;
// };
