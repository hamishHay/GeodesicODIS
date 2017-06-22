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
#include "coriolis.h"
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
    double r, omega, e, obliq, drag_coeff, h;

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
    h = globals->h.Value();

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
        // UPDATE TIDAL POTENTIAL / EXTERNAL FORCING
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

        // APPLY DRAG TO FLUID
        switch (globals->fric_type)
        {
          case LINEAR:
            linearDrag(node_num, drag_coeff, *dvel_dt_t0, *vel_tm1);
            break;
          case QUADRATIC:
            quadraticDrag(node_num, drag_coeff, h, *dvel_dt_t0, *vel_tm1);
            break;
        }

        // APPLY CORIOLIS ACCELERATION
        coriolisForce(grid, *dvel_dt_t0, *vel_tm1);


        // APPLY PRESSURE GRADIENT
        pressureGradient(grid, *dvel_dt_t0, *press_tm1);

        // INTEGRATE VELOCITIES FORWARD IN TIME
        for (i = 0; i<node_num; i++)
        {
            (*vel_t0)(i,0) = (*dvel_dt_t0)(i,0) * dt + (*vel_tm1)(i,0);
            (*vel_tm1)(i,0) = (*vel_t0)(i,0);

            (*vel_t0)(i,1) = (*dvel_dt_t0)(i,1) * dt + (*vel_tm1)(i,1);
            (*vel_tm1)(i,1) = (*vel_t0)(i,1);
        }


        // CALCULATE VELOCITY FIELD DIVERGENCE
        velocityDivergence(grid, *dpress_dt_t0, *vel_t0);

        // INTEGRATE PRESSURE FORWARD IN TIME
        for (i = 0; i<node_num; i++)
        {
            (*press_t0)(i) = (*dpress_dt_t0)(i) * dt + (*press_tm1)(i);
            (*press_tm1)(i) = (*press_t0)(i);
        }

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
int ab3Integrator(Globals * globals, Mesh * grid)
{
    std::ostringstream outstring;
    int i,j,iter,out_count;
    double ** pp, * p_t0, * p_tm1, * p_dt_t0;
    double * sinLat;

    Array2D<double> * vel_t0;      // velocity solution for current timestep
    Array2D<double> * vel_tm1;     // velocity solution at previous timestep (t minus 1)
    Array2D<double> * dvel_dt_t0;  // velocity time derivative at current timestep
    Array2D<double> * dvel_dt_tm1;  // velocity time derivative at current timestep
    Array2D<double> * dvel_dt_tm2;  // velocity time derivative at current timestep

    Array1D<double> * press_t0;      // displacement solution for current timestep
    Array1D<double> * press_tm1;     // displacement solution at previous timestep (t minus 1)
    Array1D<double> * dpress_dt_t0;  // displacement time derivative at current timestep
    Array1D<double> * dpress_dt_tm1;  // displacement time derivative at current timestep
    Array1D<double> * dpress_dt_tm2;  // displacement time derivative at current timestep

    Array1D<double> * energy_diss;
    Array1D<double> * cv_mass;


    double end_time, current_time, dt, out_frac, orbit_period, orbit_out_frac, out_time;
    double r, omega, e, obliq, drag_coeff, h;

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
    obliq = globals->theta.Value();
    h = globals->h.Value();

    node_num = globals->node_num;

    Output = globals->Output;
    pp = new double *[globals->out_tags.size()-1];

    // std::cout<<globals->out_tags[1]<<std::endl;

    // sinLat = &(grid->trigLat(0,1));

    outstring << "Defining arrays for Euler time integration..." << std::endl;

    dvel_dt_t0 = new Array2D<double>(node_num, 2);
    dvel_dt_tm1 = new Array2D<double>(node_num, 2);
    dvel_dt_tm2 = new Array2D<double>(node_num, 2);
    vel_t0 = new Array2D<double>(node_num, 2);
    vel_tm1 = new Array2D<double>(node_num, 2);

    dpress_dt_t0 = new Array1D<double>(node_num);
    dpress_dt_tm1 = new Array1D<double>(node_num);
    dpress_dt_tm2 = new Array1D<double>(node_num);
    press_t0 = new Array1D<double>(node_num);      // pressure solution for current timestep
    press_tm1 = new Array1D<double>(node_num);     // pressure solution at previous timestep (t minus 1)

    cv_mass = new Array1D<double>(node_num);

    energy_diss = new Array1D<double>(node_num);

    for (i=0; i<node_num; i++)
    {
        (*press_tm1)(i) = 0.0;
        (*cv_mass)(i) = 1000.0 * h * grid->control_volume_surf_area_map(i);

        // std::cout<<(*cv_mass)(i)<<std::endl;
    }

    double a, b, c;

    // AB3 constants for the time integration
    a = 23./12.;
    b = -16./12.;
    c = 5./12.;

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
            case OBLIQ:
                deg2Obliq(grid, *dvel_dt_t0, current_time, r, omega, obliq);
                break;

        }

        current_time += dt;
        out_time += dt;

        switch (globals->fric_type)
        {
          case LINEAR:
            linearDrag(node_num, drag_coeff, *dvel_dt_t0, *vel_tm1);
            break;
          case QUADRATIC:
            quadraticDrag(node_num, drag_coeff, h, *dvel_dt_t0, *vel_tm1);
            break;
        }


        coriolisForce(grid, *dvel_dt_t0, *vel_tm1);

        // Calculate pressure gradient
        pressureGradient(grid, *dvel_dt_t0, *press_tm1);

        // Explicit time integration
        if (iter >= 2)
        {
            for (i = 0; i<node_num; i++)
            {
                (*vel_t0)(i,0) = (a*(*dvel_dt_t0)(i,0)
                                  + b*(*dvel_dt_tm1)(i,0)
                                  + c*(*dvel_dt_tm2)(i,0)) * dt
                                 + (*vel_tm1)(i,0);

                (*vel_t0)(i,1) = (a*(*dvel_dt_t0)(i,1)
                                  + b*(*dvel_dt_tm1)(i,1)
                                  + c*(*dvel_dt_tm2)(i,1)) * dt
                                 + (*vel_tm1)(i,1);

                (*vel_tm1)(i, 0) = (*vel_t0)(i, 0);
                (*vel_tm1)(i, 1) = (*vel_t0)(i, 1);

                (*dvel_dt_tm2)(i, 0) = (*dvel_dt_tm1)(i, 0);
                (*dvel_dt_tm2)(i, 1) = (*dvel_dt_tm1)(i, 1);

                (*dvel_dt_tm1)(i, 0) = (*dvel_dt_t0)(i, 0);
                (*dvel_dt_tm1)(i, 1) = (*dvel_dt_t0)(i, 1);

                (*energy_diss)(i) = drag_coeff * (*cv_mass)(i)
                                    * ((*vel_t0)(i,0)*(*vel_t0)(i,0) + (*vel_t0)(i,1)*(*vel_t0)(i,1));
            }
        }
        else if (iter == 1)
        {
            for (i = 0; i<node_num; i++)
            {
                (*vel_t0)(i,0) = (*dvel_dt_t0)(i,0) * dt + (*vel_tm1)(i,0);
                (*vel_tm1)(i,0) = (*vel_t0)(i,0);
                (*dvel_dt_tm1)(i,0) = (*dvel_dt_t0)(i,0);

                (*vel_t0)(i,1) = (*dvel_dt_t0)(i,1) * dt + (*vel_tm1)(i,1);
                (*vel_tm1)(i,1) = (*vel_t0)(i,1);
                (*dvel_dt_tm1)(i,1) = (*dvel_dt_t0)(i,1);

                (*energy_diss)(i) = drag_coeff * (*cv_mass)(i) * ((*vel_t0)(i,0)*(*vel_t0)(i,0) + (*vel_t0)(i,1)*(*vel_t0)(i,1));
            }
        }
        else
        {
            for (i = 0; i<node_num; i++)
            {
                (*vel_t0)(i,0) = (*dvel_dt_t0)(i,0) * dt + (*vel_tm1)(i,0);
                (*vel_tm1)(i,0) = (*vel_t0)(i,0);
                (*dvel_dt_tm2)(i,0) = (*dvel_dt_t0)(i,0);

                (*vel_t0)(i,1) = (*dvel_dt_t0)(i,1) * dt + (*vel_tm1)(i,1);
                (*vel_tm1)(i,1) = (*vel_t0)(i,1);
                (*dvel_dt_tm2)(i,1) = (*dvel_dt_t0)(i,1);

                (*energy_diss)(i) = drag_coeff * (*cv_mass)(i) * ((*vel_t0)(i,0)*(*vel_t0)(i,0) + (*vel_t0)(i,1)*(*vel_t0)(i,1));

            }
        }



        // Calculate velocity divergence
        velocityDivergence(grid, *dpress_dt_t0, *vel_t0);


        if (iter >= 2)
        {
            for (i = 0; i<node_num; i++)
            {
                (*press_t0)(i) = (a*(*dpress_dt_t0)(i)
                                 + b*(*dpress_dt_tm1)(i)
                                 + c*(*dpress_dt_tm2)(i)) * dt
                                  + (*press_tm1)(i);
                (*press_tm1)(i) = (*press_t0)(i);

                (*dpress_dt_tm2)(i) = (*dpress_dt_tm1)(i);
                (*dpress_dt_tm1)(i) = (*dpress_dt_t0)(i);

            }
        }
        else if (iter == 1)
        {
            for (i = 0; i<node_num; i++)
            {
                (*press_t0)(i) = (*dpress_dt_t0)(i) * dt + (*press_tm1)(i);
                (*press_tm1)(i) = (*press_t0)(i);
                (*dpress_dt_tm1)(i) = (*dpress_dt_t0)(i);
            }
        }
        else
        {
            for (i = 0; i<node_num; i++)
            {
                (*press_t0)(i) = (*dpress_dt_t0)(i) * dt + (*press_tm1)(i);
                (*press_tm1)(i) = (*press_t0)(i);
                (*dpress_dt_tm2)(i) = (*dpress_dt_t0)(i);
            }
        }

        // Check for output
        iter ++;

        if (out_time >= out_frac*orbit_period)
        {
            std::cout<<std::fixed << std::setprecision(8) <<"DUMPING DATA AT "<<current_time/orbit_period<<std::endl;

            out_time -= out_frac*orbit_period;
            pp[2] = &(*energy_diss)(0);
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
