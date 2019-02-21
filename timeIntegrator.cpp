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
#include "initialConditions.h"
#include <math.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <signal.h>

int flag = 0;

void CatchExit(int sig) {
    printf("%s\n", "Caught Terminate Signal...");
    flag = 1;
}

int updateVelocity(Globals * globals, Mesh * grid, Array2D<double> & dvdt, Array2D<double> & v_tm1, Array1D<double> & p_tm1, double current_time)
{
    double r, omega, e, obliq, h, drag_coeff, visc, g, dt;
    int node_num;

    r =       globals->radius.Value();
    g =       globals->g.Value();
    dt =      globals->timeStep.Value();
    omega =   globals->angVel.Value();
    e =       globals->e.Value();
    drag_coeff = globals->alpha.Value();
    obliq =   globals->theta.Value();
    h =       globals->h.Value();
    node_num = globals->node_num;

    visc = 1e4;

    for (int i=0; i <node_num; i++)
    {
        dvdt(i, 0) = 0.0;
        dvdt(i, 1) = 0.0;
    }

    switch (globals->tide_type)
    {
        case ECC:
            deg2Ecc(grid, dvdt, current_time, r, omega, e);
            break;
        case ECC_LIB:
            deg2EccLib(grid, dvdt, current_time, r, omega, e);
            break;
        case ECC_WEST:
            deg2EccWest(grid, dvdt, current_time, r, omega, e);
            break;
        case ECC_EAST:
            deg2EccEast(grid, dvdt, current_time, r, omega, e);
            break;
        case ECC_RAD:
            deg2EccRad(grid, dvdt, current_time, r, omega, e);
            break;
        case OBLIQ:
            deg2Obliq(grid, dvdt, current_time, r, omega, obliq);
            break;
        case OBLIQ_WEST:
            deg2ObliqWest(grid, dvdt, current_time, r, omega, obliq);
            break;
        case OBLIQ_EAST:
            deg2ObliqEast(grid, dvdt, current_time, r, omega, obliq);
            break;
        case FULL:
            deg2Full(grid, dvdt, current_time, r, omega, obliq, e);
            // deg2Ecc(grid, dvdt, current_time, r, omega, e);
            // deg2Obliq(grid, dvdt, current_time, r, omega, obliq);
            break;
        case PLANET:
            deg2Planet(grid, dvdt, current_time, r);
            break;
        case PLANET_OBL:
            deg2PlanetObl(grid, dvdt, current_time, r, obliq);
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

    switch (globals->surface_type)
    {
        case FREE:
            pressureGradient(grid, dvdt, p_tm1, node_num, g);
            // pressureGradient(grid, dvdt, p_tm1, g);
            break;

        case FREE_LOADING:
            pressureGradient(grid, dvdt, p_tm1, node_num, g);
            pressureGradientSH(globals, grid, dvdt, p_tm1, g);
            break;

        case LID_LOVE:
            pressureGradient(grid, dvdt, p_tm1, node_num, g);
            pressureGradientSH(globals, grid, dvdt, p_tm1, g);
            break;

        case LID_MEMBR:
            pressureGradient(grid, dvdt, p_tm1, node_num, g);
            pressureGradientSH(globals, grid, dvdt, p_tm1, g);
            break;

        case LID_INF:
            // pressureGradient(grid, dvdt, p_tm1, 1.0/1000.0);
            // pressureGradientSH(globals, grid, dvdt, p_tm1, -1.0/1000.0);
            break;

    }

    velocityDiffusion(grid, dvdt, v_tm1, visc);

    return 1;

};

int updateDisplacement(Globals * globals, Mesh * grid, Array1D<double> & deta_dt, Array2D<double> & v_t0)
{
    double sum = -1.0;

    velocityDivergence(grid, deta_dt, v_t0, sum, globals->h.Value());

    return 1;
};

// Function to implement the time loop to solve mass and momentum equations over
// the geodesic grid. The time integration method is the explicit Euler method.
int eulerExplicit(Globals * globals, Mesh * grid)
{
    std::ostringstream outstring;
    int i,j,iter,out_count;
    double ** pp;

    Array2D<double> * vel_t0;      // velocity solution for current timestep
    Array2D<double> * vel_tm1;     // velocity solution at previous timestep (t minus 1)
    Array2D<double> * dvel_dt_t0;  // velocity time derivative at current timestep

    Array1D<double> * press_t0;      // displacement solution for current timestep
    Array1D<double> * press_tm1;     // displacement solution at previous timestep (t minus 1)
    Array1D<double> * dpress_dt_t0;  // displacement time derivative at current timestep

    double end_time, current_time, dt, out_frac, orbit_period, out_time;
    // double r, omega, e, obliq, drag_coeff, h;

    int node_num;

    OutFiles * Output;

    end_time = globals->endTime.Value();
    dt = globals->timeStep.Value();
    out_frac = globals->outputTime.Value();
    orbit_period = globals->period.Value();

    node_num = globals->node_num;

    Output = globals->Output;
    pp = new double *[globals->out_tags.size()-1];

    std::cout<<globals->out_tags[0]<<std::endl;

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

        current_time += dt;
        out_time += dt;

        // Check for output
        iter ++;

        if (out_time >= out_frac*orbit_period)
        {
            std::cout<<std::fixed << std::setprecision(8) <<"DUMPING DATA AT "<<current_time/orbit_period<<std::endl;
            // std::cout<<" AVG DISS: "<<*total_diss<<" W"<<std::endl;

            out_time -= out_frac*orbit_period;
            // pp[3] = &(*energy_diss)(0);
            pp[1] = &(*press_t0)(0);
            // pp[3] = &(*press_dummy)(0);
            // pp[1] = &(grid->land_mask(0));
            pp[2] = &(*vel_t0)(0,0);
            // pp[2] = &(*vel_dummy)(0,0);
            // pp[0] = &total_diss[0];

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

    // Array2D<double> * vel_dummy;

    Array1D<double> * press_t0;      // displacement solution for current timestep
    Array1D<double> * press_tm1;     // displacement solution at previous timestep (t minus 1)
    Array1D<double> * press_tm2;
    Array1D<double> * dpress_dt_t0;  // displacement time derivative at current timestep
    Array1D<double> * dpress_dt_tm1;  // displacement time derivative at current timestep
    Array1D<double> * dpress_dt_tm2;  // displacement time derivative at current timestep

    // Array1D<double> * press_dummy;

    Array1D<double> * energy_diss;
    Array1D<double> * cv_mass;

    // Array2D<int> * friend_list;
    // friend_list = &(grid->node_friends);

    double end_time, current_time, dt, out_frac, orbit_period, out_time;
    // double r, omega, e, obliq, drag_coeff, h;
    double lon, lat;
    double r = globals->radius.Value();
    double omega = globals->radius.Value();
    double obl = globals->theta.Value();

    int node_num;

    OutFiles * Output;

    signal(SIGINT, CatchExit);

    //--------------------------------------------------------------------------
    //------------------------- ASSIGN VARIABLES -------------------------------
    //--------------------------------------------------------------------------

    end_time = globals->endTime.Value();
    dt = globals->timeStep.Value();
    out_frac = globals->outputTime.Value();
    orbit_period = globals->period.Value();

    node_num = globals->node_num;

    Output = globals->Output;
    pp = new double *[globals->out_tags.size()];

    std::vector<std::string> * tags;
    tags = &globals->out_tags;

    double * total_diss;
    double e_diss;
    total_diss = new double[1];



    outstring << "Defining arrays for Adams-Bashforth time integration..." << std::endl;

    dvel_dt_t0 = new Array2D<double>(node_num, 2);
    dvel_dt_tm1 = new Array2D<double>(node_num, 2);
    dvel_dt_tm2 = new Array2D<double>(node_num, 2);
    vel_t0 = new Array2D<double>(node_num, 2);
    vel_tm1 = new Array2D<double>(node_num, 2);

    // vel_dummy = new Array2D<double>(node_num, 2);
    // press_dummy = new Array1D<double>(node_num);

    dpress_dt_t0 = new Array1D<double>(node_num);
    dpress_dt_tm1 = new Array1D<double>(node_num);
    dpress_dt_tm2 = new Array1D<double>(node_num);
    press_t0 = new Array1D<double>(node_num);      // pressure solution for current timestep
    press_tm1 = new Array1D<double>(node_num);     // pressure solution at previous timestep (t minus 1)
    press_tm2 = new Array1D<double>(node_num);
    cv_mass = new Array1D<double>(node_num);

    energy_diss = new Array1D<double>(node_num);

    for (i=0; i<globals->out_tags.size(); i++)
    {
        if ((*tags)[i] == "velocity output")             pp[i] = &(*vel_t0)(0,0);
        else if ((*tags)[i] == "displacement output")    pp[i] = &(*press_t0)(0);
        else if ((*tags)[i] == "dissipation output")     pp[i] = &(*energy_diss)(0);
        else if ((*tags)[i] == "dissipation avg output") pp[i] = &total_diss[0];
        else if ((*tags)[i] == "kinetic avg output")     pp[i] = &current_time;
    }

    for (i=0; i<node_num; i++)
    {
        (*press_tm1)(i) = 0.0;
        (*vel_t0)(i,0) = 0.0;
        (*vel_t0)(i,1) = 0.0;
        (*vel_tm1)(i,0) = 0.0;//cosLat[i*2] * 1e-6;
        (*vel_tm1)(i,1) = 0.0;
        (*press_t0)(i) = 0.0;
        (*press_tm1)(i) = 0.0;
        (*press_tm2)(i) = 0.0;
        (*cv_mass)(i) = grid->control_volume_mass(i);
    }


    if (globals->init.Value())
    {
       std::cout<<"Using initial conditions"<<std::endl;
       loadInitialConditions(globals, grid,
                             *vel_tm1,   *dvel_dt_t0,   *dvel_dt_tm1,   *dvel_dt_tm2,
                             *press_tm1, *dpress_dt_t0, *dpress_dt_tm1, *dpress_dt_tm2);

        // smoothingSH(globals, grid, *press_tm1);
        // smoothingSH(globals, grid, *dpress_dt_t0);
        // smoothingSH(globals, grid, *dpress_dt_tm1);
        // smoothingSH(globals, grid, *dpress_dt_tm2);
        // smoothingSHVector(globals, grid, *vel_tm1);
        // smoothingSHVector(globals, grid, *dvel_dt_t0);
        // smoothingSHVector(globals, grid, *dvel_dt_tm1);
        // smoothingSHVector(globals, grid, *dvel_dt_tm2);

       iter = 3;
    }
    else iter = 0;


    //--------------------------------------------------------------------------
    //------------------------- BEGIN TIME LOOP --------------------------------
    //--------------------------------------------------------------------------

    out_count = 1;
    current_time = 0.0;
    out_time = 0.0;

    int err;
    double sum_dummy = 1.0;
    double avg_val;
    int f;
    double g_old = globals->g.Value();
    bool lid_switch = false;
    // double cd_old = globals->alpha.Value();
    // double cd_max = 1e-1;
    // globals->alpha.SetValue(cd_max);

    // if we want to use the ice shell, start with free surface solution that
    // uses the approximate ice shell solution (b_2 in the pressure gradient
    // term).
    // if ((globals->surface_type == LID_LOVE) && (!globals->init.Value())) {
    //     lid_switch = true; // switch is only set to true if we want to turn on
    //                        // the ice shell later.
    //     globals->surface_type = FREE;
    //
    //     // modify gravity by the shell factor beta_2, which modifies the
    //     // restoring force on the ocean and assumes the ocean response is
    //     // governed by degree-2.
    //     globals->g.SetValue(-g_old*(globals->shell_factor_beta[2]-1.0));
    // }


    // std::cout<<g_old<<' '<<globals->g.Value()<<' '<<-(globals->shell_factor_beta[2]-1.0)<<std::endl;

    while (current_time <= end_time)
    {
        current_time += dt;
        out_time += dt;

        // if (globals->alpha.Value() > cd_old + 0.01*cd_old) globals->alpha.SetValue(cd_max*exp(-current_time/(4.*orbit_period)) + cd_old);
        // else globals->alpha.SetValue(cd_old);
        // std::cout<<globals->alpha.Value()<<std::endl;

        // check if we want to turn on the ice shell
        // if ((lid_switch) && (globals->surface_type == FREE) && (current_time >= 0.01*end_time))
        // {
        //     globals->surface_type = LID_LOVE;
        //     globals->g.SetValue(g_old);
        //     lid_switch = false;
        //     // std::cout<<g_old<<' '<<globals->g.Value()<<' '<<-(globals->shell_factor_beta[2]-1.0)<<std::endl;
        // }


        if (globals->surface_type == FREE ||
            globals->surface_type == FREE_LOADING ||
            globals->surface_type == LID_LOVE ||
            globals->surface_type == LID_MEMBR)
        {
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

            // if (out_time >= out_frac*orbit_period) smoothingSH(globals, grid, *press_tm1);
            // if (out_time >= out_frac*orbit_period) smoothingSHVector(globals, grid, *vel_tm1);
        }

        else if (globals->surface_type == LID_INF)
        {

            // TEST VERSION!!!!!!!!!!!!!!!! STILL UNDER TESTING. DON'T TRUST!

            // SOLVE THE MOMENTUM EQUATION
            updateVelocity(globals, grid, *dvel_dt_t0, *vel_tm1, *press_t0, current_time);


            // INTEGRATE VELOCITIES FORWARD IN TIME
            // for (i = 0; i<node_num; i++)
            // {
            //    (*vel_t0)(i,0) = (*dvel_dt_t0)(i,0) * dt + (*vel_tm1)(i,0);
            //    (*vel_tm1)(i,0) = (*vel_t0)(i,0);
            //
            //    (*vel_t0)(i,1) = (*dvel_dt_t0)(i,1) * dt + (*vel_tm1)(i,1);
            //    (*vel_tm1)(i,1) = (*vel_t0)(i,1);
            // }
            // MARCH FORWARD VELOCITY SOLUTION
            // integrateAB3vector(globals, grid, *vel_t0, *vel_tm1, *dvel_dt_t0, *dvel_dt_tm1, *dvel_dt_tm2, iter);

            // FORCE CONTINUITY EQUATION
            err = updatePressure(globals, grid, *press_t0, *vel_t0, *dvel_dt_t0);

            // pressureGradientSH(globals, grid, *dvel_dt_t0, *press_t0, -1.0/1000.0);

            // INTEGRATE VELOCITIES FORWARD IN TIME
            for (i = 0; i<node_num; i++)
            {
                (*vel_t0)(i,0) = (*dvel_dt_t0)(i,0) * dt + (*vel_tm1)(i,0);
                (*vel_tm1)(i,0) = (*vel_t0)(i,0);

                (*vel_t0)(i,1) = (*dvel_dt_t0)(i,1) * dt + (*vel_tm1)(i,1);
                (*vel_tm1)(i,1) = (*vel_t0)(i,1);

                // (*G_terms_tm1)(i,0) = (*G_terms_t0)(i,0);
                // (*G_terms_tm1)(i,1) = (*G_terms_t0)(i,1);
            }

            updateEnergy(globals, e_diss, *energy_diss, *vel_t0, *cv_mass);

            total_diss[0] = e_diss;


        }

        // Check for output
        iter ++;

        if (out_time >= out_frac*orbit_period)
        {
            outstring << std::fixed << std::setprecision(8) <<"DUMPING DATA AT "<<current_time/orbit_period;
            outstring << " AVG DISS: "<<*total_diss*4*pi*r*r/1e9<<" GW"<<out_count;

            Output->Write(OUT_MESSAGE, &outstring);

            out_time -= out_frac*orbit_period;

            Output->DumpData(globals, out_count, pp);
            out_count++;

        }

        if (flag) {
            outstring << "Terminate signal caught..."<<std::endl;
            Output->Write(OUT_MESSAGE, &outstring);


            break;
            // return 0;

        }
    }

    // Save initial conditions!
    // TODO - move the below code to a function in initialConditions.cpp


    FILE * outFile;
    std::string initial_condition_file;
    initial_condition_file = Output->path + SEP + "init_new.txt";
    remove(&initial_condition_file[0]);
    outFile = fopen(&initial_condition_file[0], "w");

    for (i=0; i<node_num; i++)
    {                   //   u    dudt0  dudt1  dudt2    v    dvdt0  dvdt1  dvdt2    p    dpdt0  dpdt1  dpdt2
        fprintf (outFile, "%1.6E, %1.6E, %1.6E, %1.6E, %1.6E, %1.6E, %1.6E, %1.6E, %1.6E, %1.6E, %1.6E, %1.6E\n",
                 (*vel_t0)(i,0), (*dvel_dt_t0)(i, 0), (*dvel_dt_tm1)(i, 0), (*dvel_dt_tm2)(i, 0),
                 (*vel_t0)(i,1), (*dvel_dt_t0)(i, 1), (*dvel_dt_tm1)(i, 1), (*dvel_dt_tm2)(i, 1),
                 (*press_t0)(i), (*dpress_dt_t0)(i),  (*dpress_dt_tm1)(i),  (*dpress_dt_tm2)(i) );
    }

    fclose(outFile);

    Output->Write(OUT_MESSAGE, &outstring);

    if (flag) return 1;     // return with error
    else return 0;          // return without error
};
