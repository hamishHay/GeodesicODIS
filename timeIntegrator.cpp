#include "mesh.h"
#include "globals.h"
#include "outFiles.h"
#include "solver.h"
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include "tidalPotentials.h"
#include "spatialOperators.h"
#include "temporalOperators.h"
#include "energy.h"
#include "pressure.h"
#include "drag.h"
#include "coriolis.h"
#include <math.h>
#include <iostream>
#include <sstream>
#include <iomanip>

int updateVelocity(Globals * globals, Mesh * grid, Array2D<double> & dvdt, Array2D<double> & v_tm1, Array1D<double> & p_tm1, double current_time)
{
    double r, omega, e, obliq, h, drag_coeff, visc;
    int node_num;

    r = globals->radius.Value();
    omega = globals->angVel.Value();
    e = globals->e.Value();
    drag_coeff = globals->alpha.Value();
    obliq = globals->theta.Value();
    h = globals->h.Value();
    node_num = globals->node_num;
    visc = 1e3;

    switch (globals->tide_type)
    {
        case ECC:
            deg2Ecc(grid, dvdt, current_time, r, omega, e);
            break;
        case ECC_RAD:
            deg2EccRad(grid, dvdt, current_time, r, omega, e);
            break;
        case OBLIQ:
            deg2Obliq(grid, dvdt, current_time, r, omega, obliq);
            break;
    }


    switch (globals->fric_type)
    {
        case LINEAR:
            linearDrag(node_num, drag_coeff, dvdt, v_tm1);
            break;
        case QUADRATIC:
            quadraticDrag(node_num, drag_coeff, h, dvdt, v_tm1);
            break;
    }


    coriolisForce(grid, dvdt, v_tm1);

    pressureGradient(grid, dvdt, p_tm1);

    velocityDiffusion(grid, dvdt, v_tm1, visc);

};

int updateDisplacement(Globals * globals, Mesh * grid, Array1D<double> & deta_dt, Array2D<double> & v_t0)
{
    double sum = -1.0;
    velocityDivergence(grid, deta_dt, v_t0, sum);
};

// Function to implement the time loop to solve mass and momentum equations over
// the geodesic grid. The time integration method is the explicit Euler method.
int eulerExplicit(Globals * globals, Mesh * grid)
{
    std::ostringstream outstring;
    int i,j,iter,out_count;
    double ** pp;
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
        (*press_tm1)(i) = 0.0;
        (*vel_t0)(i,0) = 0.0;
        (*vel_t0)(i,1) = 0.0;
        (*press_t0)(i) = 0.0;
        // (*cv_mass)(i) = grid->control_volume_mass(i);
    }

    iter = 0;
    out_count = 1;
    current_time = 0.0;
    out_time = 0.0;
    while (current_time <= end_time)
    {
        updateVelocity(globals, grid, *dvel_dt_t0, *vel_tm1, *press_tm1, current_time);

        // INTEGRATE VELOCITIES FORWARD IN TIME
        for (i = 0; i<node_num; i++)
        {
            (*vel_t0)(i,0) = (*dvel_dt_t0)(i,0) * dt + (*vel_tm1)(i,0);
            (*vel_tm1)(i,0) = (*vel_t0)(i,0);

            (*vel_t0)(i,1) = (*dvel_dt_t0)(i,1) * dt + (*vel_tm1)(i,1);
            (*vel_tm1)(i,1) = (*vel_t0)(i,1);
        }

        updateDisplacement(globals, grid, *dpress_dt_t0, *vel_t0);

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
int ab3Explicit(Globals * globals, Mesh * grid)
{
    //--------------------------------------------------------------------------
    //------------------------- DECLARE VARIABLES ------------------------------
    //--------------------------------------------------------------------------

    std::ostringstream outstring;
    int i,j,k,iter,out_count;
    double ** pp;

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

    //--------------------------------------------------------------------------
    //------------------------- ASSIGN VARIABLES -------------------------------
    //--------------------------------------------------------------------------

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
    pp = new double *[globals->out_tags.size()];

    // std::cout<<globals->out_tags.size()-1<<std::endl;
    // std::cout<<globals->out_tags[0]<<std::endl;
    // std::cout<<globals->out_tags[1]<<std::endl;
    // std::cout<<globals->out_tags[2]<<std::endl;
    // std::cout<<globals->out_tags[3]<<std::endl;

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

    double * total_diss;
    double e_diss;
    total_diss = new double[1];

    for (i=0; i<node_num; i++)
    {
        (*press_tm1)(i) = 0.0;
        (*vel_t0)(i,0) = 0.0;
        (*vel_t0)(i,1) = 0.0;
        (*press_t0)(i) = 0.0;
        (*cv_mass)(i) = grid->control_volume_mass(i);
    }

    //--------------------------------------------------------------------------
    //------------------------- BEGIN TIME LOOP --------------------------------
    //--------------------------------------------------------------------------

    iter = 0;
    out_count = 1;
    current_time = 0.0;
    out_time = 0.0;

    while (current_time <= end_time)
    {
        current_time += dt;
        out_time += dt;

        switch (globals->surface_type)
        {
            case FREE:
                // SOLVE THE MOMENTUM EQUATION
                updateVelocity(globals, grid, *dvel_dt_t0, *vel_tm1, *press_tm1, current_time);

                // MARCH VELOCITY FORWARD IN TIME
                integrateAB3vector(globals, grid, *vel_t0, *vel_tm1, *dvel_dt_t0, *dvel_dt_tm1, *dvel_dt_tm2, iter);

                // UPDATE ENERGY DISSIPATION
                updateEnergy(globals, e_diss, *energy_diss, *vel_t0, *cv_mass);
                total_diss[0] = e_diss;

                // SOLVE THE CONTINUITY EQUATION
                updateDisplacement(globals, grid, *dpress_dt_t0, *vel_t0);

                // MARCH DISPLACEMENT FORWARD IN TIME
                integrateAB3scalar(globals, grid, *press_t0, *press_tm1, *dpress_dt_t0, *dpress_dt_tm1, *dpress_dt_tm2, iter);

                break;

            case INF_LID:
                // SOLVE THE MOMENTUM EQUATION
                updateVelocity(globals, grid, *dvel_dt_t0, *vel_tm1, *press_tm1, current_time);

                // MARCH FORWARD VELOCITY SOLUTION
                integrateAB3vector(globals, grid, *vel_t0, *vel_tm1, *dvel_dt_t0, *dvel_dt_tm1, *dvel_dt_tm2, iter);

                // FORCE CONTINUITY EQUATION
                updatePressure(globals, grid, *press_t0, *vel_t0);

                    // CHECK VELOCITY FIELD DIVERGENCE

                    // ---- CALCULATE PRESSURE CORRECTION

                    // ---- UPDATE PRESSURE FIELD

                    // ---- UPDATE VELOCITY FIELD

                    // ---- REPEAT UNTIL DIVERGENCE IS ZERO

                // UPDATE ENERGY DISSIPATION

                break;

            case LID:
                outstring << "LID SOLVER COMING SOON!" << std::endl;
                Output->Write(ERR_MESSAGE, &outstring);
                Output->TerminateODIS();
                break;

        }

        if (iter == 1000) Output->TerminateODIS();

        // Check for output
        iter ++;

        if (out_time >= out_frac*orbit_period)
        {
            std::cout<<std::fixed << std::setprecision(8) <<"DUMPING DATA AT "<<current_time/orbit_period;
            std::cout<<" TOTAL DISS: "<<*total_diss<<" GW"<<std::endl;

            out_time -= out_frac*orbit_period;
            pp[3] = &(*energy_diss)(0);
            pp[1] = &(*press_t0)(0);
            // pp[1] = &(grid->land_mask(0));
            pp[2] = &(*vel_t0)(0,0);
            pp[0] = &total_diss[0];

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
