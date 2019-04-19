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
#include "initialConditions.h"
#include <math.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <signal.h>


#include <mkl.h>
#include <omp.h>

#include <sys/time.h>

int flag = 0;

void CatchExit(int sig) {
    printf("%s\n", "Caught Terminate Signal...");
    flag = 1;
}


typedef unsigned long long timestamp_t;
double total_time=0.0;
int time_iter=0;

static timestamp_t get_timestamp ()
{
  struct timeval now;
  gettimeofday (&now, NULL);
  return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
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

    visc = 5e3;

    Array1D<double> forcing_potential(node_num);

    forcing(globals, grid, forcing_potential, globals->tide_type, current_time, e, obliq);

    // TODO - Move these pressure grad calls to a function that makes sense
    switch (globals->surface_type)
    {
        case FREE_LOADING:
            pressureGradientSH(globals, grid, dvdt, p_tm1, forcing_potential, g);
            break;

        case LID_LOVE:
            pressureGradientSH(globals, grid, dvdt, p_tm1, forcing_potential, g);
            break;

        case LID_MEMBR:
            pressureGradientSH(globals, grid, dvdt, p_tm1, forcing_potential, g);
            break;

        case LID_INF:
            // pressureGradient(grid, dvdt, p_tm1, 1.0/1000.0);
            // pressureGradientSH(globals, grid, dvdt, p_tm1, -1.0/1000.0);
            break;
    }

    sparse_operation_t operation = SPARSE_OPERATION_NON_TRANSPOSE;
    matrix_descr descript;
    descript.type = SPARSE_MATRIX_TYPE_GENERAL;

    double alpham = 1.0;
    double betam = 0.0;

    Array2D<double> soln_old(node_num, 3);

    #pragma omp parallel for
    for (int i=0; i<node_num; ++i)
    {
        soln_old(i, 0) = v_tm1(i, 0);
        soln_old(i, 1) = v_tm1(i, 1);
        soln_old(i, 2) = p_tm1(i) - forcing_potential(i)/g;
    }

    // timestamp_t t0 = get_timestamp();

    // Perform A*d = dvdt where A is operatorMomentum, and d = [u, v, eta-U/g]
    mkl_sparse_d_mv(operation, 1.0, *(grid->operatorMomentum), descript, &(soln_old(0,0)), betam, &(dvdt(0,0)));

    // timestamp_t t1 = get_timestamp();
    //
    // total_time += (t1 - t0) / 1000000.0L;
    //
    // time_iter += 1;

    // std::cout<<"TIME: "<<total_time/((double)time_iter)<<std::endl;

    if (globals->fric_type == QUADRATIC) quadraticDrag(node_num, drag_coeff, h, dvdt, v_tm1);

    return 1;

};

int updateDisplacement(Globals * globals, Mesh * grid, Array1D<double> & deta_dt, Array2D<double> & v_t0)
{
    double sum = -1.0;

    velocityDivergence(grid, deta_dt, v_t0, sum, globals->h.Value());

    return 1;
};

// Function to implement the time loop to solve mass and momentum equations over
// the geodesic grid using the Adams-Bashforth 3rd order forward method.
int ab3Explicit(Globals * globals, Mesh * grid)
{
    //--------------------------------------------------------------------------
    //------------------------- DECLARE VARIABLES ------------------------------
    //--------------------------------------------------------------------------

    omp_set_num_threads( globals->core_num.Value() );
    mkl_set_num_threads( globals->core_num.Value() );

    std::ostringstream outstring;
    int i, j, k, iter, out_count;
    double ** pp;

    Array2D<double> * v_t0;        // velocity solution for current timestep
    Array2D<double> * dv_dt_t0;    // velocity time derivative at current timestep

    Array1D<double> * p_t0;      // displacement solution for current timestep
    Array1D<double> * dp_dt_t0;  // displacement time derivative at current timestep

    Array3D<double> * dv_dt;
    Array2D<double> * dp_dt;

    Array1D<double> * energy_diss;
    Array1D<double> * cv_mass;

    double end_time, current_time, dt, out_frac, orbit_period, out_time;

    double r = globals->radius.Value();

    int node_num;
    node_num = globals->node_num;

    OutFiles * Output;

    signal(SIGINT, CatchExit);

    //--------------------------------------------------------------------------
    //------------------------- ASSIGN VARIABLES -------------------------------
    //--------------------------------------------------------------------------

    end_time     = globals->endTime.Value();
    dt           = globals->timeStep.Value();
    out_frac     = globals->outputTime.Value();
    orbit_period = globals->period.Value();

    Output = globals->Output;
    pp = new double *[globals->out_tags.size()];

    std::vector<std::string> * tags;
    tags = &globals->out_tags;

    double * total_diss;
    total_diss = new double[1];
    double e_diss;

    outstring << "Defining arrays for Adams-Bashforth time integration..." << std::endl;

    v_t0     = new Array2D<double>(node_num, 2);
    dv_dt_t0 = new Array2D<double>(node_num, 2);

    p_t0     = new Array1D<double>(node_num);
    dp_dt_t0 = new Array1D<double>(node_num);

    dv_dt    = new Array3D<double>(node_num, 2, 3);   //[node][vel component][dvdt]
    dp_dt    = new Array2D<double>(node_num, 3);      //[node][dpdt]

    cv_mass = new Array1D<double>(node_num);
    energy_diss = new Array1D<double>(node_num);

    for (i=0; i<globals->out_tags.size(); i++)
    {
        if      ((*tags)[i] == "velocity output")        pp[i] = &(*v_t0)(0,0);
        else if ((*tags)[i] == "displacement output")    pp[i] = &(*p_t0)(0);
        else if ((*tags)[i] == "dissipation output")     pp[i] = &(*energy_diss)(0);
        else if ((*tags)[i] == "dissipation avg output") pp[i] = &total_diss[0];
        else if ((*tags)[i] == "kinetic avg output")     pp[i] = &current_time;
    }

    for (i=0; i<node_num; i++)
    {
        (*v_t0)(i,0) = 0.0;
        (*v_t0)(i,1) = 0.0;
        (*p_t0)(i) = 0.0;

        (*cv_mass)(i) = grid->control_volume_mass(i);
    }

    if (globals->init.Value())
    {
       loadInitialConditions(globals, grid, *v_t0, *dv_dt, *p_t0, *dp_dt);

       iter = 3;
    }
    else iter = 0;

    //--------------------------------------------------------------------------
    //------------------------- BEGIN TIME LOOP --------------------------------
    //--------------------------------------------------------------------------

    out_count = 1;
    current_time = 0.0;
    out_time = 0.0;

    while (current_time <= end_time)
    {
        current_time += dt;
        out_time += dt;

        if (globals->surface_type == FREE ||
            globals->surface_type == FREE_LOADING ||
            globals->surface_type == LID_LOVE ||
            globals->surface_type == LID_MEMBR)
        {
            // SOLVE THE MOMENTUM EQUATION
            updateVelocity(globals, grid, *dv_dt_t0, *v_t0, *p_t0, current_time);

            #pragma omp parallel for
            for (i=0; i<node_num; ++i) {
              (*dv_dt)(i, 0, 0) = (*dv_dt_t0)(i, 0);
              (*dv_dt)(i, 1, 0) = (*dv_dt_t0)(i, 1);
            }

            // MARCH VELOCITY FORWARD IN TIME
            integrateAB3vector(globals, grid, *v_t0, *dv_dt, iter);

            // UPDATE ENERGY DISSIPATION
            updateEnergy(globals, e_diss, *energy_diss, *v_t0, *cv_mass);
            total_diss[0] = e_diss;

            // SOLVE THE CONTINUITY EQUATION
            updateDisplacement(globals, grid, *dp_dt_t0, *v_t0);

            #pragma omp parallel for
            for (i=0; i<node_num; ++i) {
              (*dp_dt)(i, 0) = (*dp_dt_t0)(i);
            }

            integrateAB3scalar(globals, grid, *p_t0, *dp_dt, iter);
        }

        // else if (globals->surface_type == LID_INF)
        // {
        //
        //     // TEST VERSION!!!!!!!!!!!!!!!! STILL UNDER TESTING. DON'T TRUST!
        //
        //     // SOLVE THE MOMENTUM EQUATION
        //     updateVelocity(globals, grid, *dv_dt_t0, *v_tm1, *p_t0, current_time);
        //
        //
        //     // INTEGRATE VELOCITIES FORWARD IN TIME
        //     // for (i = 0; i<node_num; i++)
        //     // {
        //     //    (*v_t0)(i,0) = (*dv_dt_t0)(i,0) * dt + (*v_tm1)(i,0);
        //     //    (*v_tm1)(i,0) = (*v_t0)(i,0);
        //     //
        //     //    (*v_t0)(i,1) = (*dv_dt_t0)(i,1) * dt + (*v_tm1)(i,1);
        //     //    (*v_tm1)(i,1) = (*v_t0)(i,1);
        //     // }
        //     // MARCH FORWARD VELOCITY SOLUTION
        //     // integrateAB3vector(globals, grid, *v_t0, *v_tm1, *dv_dt_t0, *dv_dt_tm1, *dv_dt_tm2, iter);
        //
        //     // FORCE CONTINUITY EQUATION
        //     err = updatePressure(globals, grid, *p_t0, *v_t0, *dv_dt_t0);
        //
        //     // pressureGradientSH(globals, grid, *dv_dt_t0, *p_t0, -1.0/1000.0);
        //
        //     // INTEGRATE VELOCITIES FORWARD IN TIME
        //     for (i = 0; i<node_num; i++)
        //     {
        //         (*v_t0)(i,0) = (*dv_dt_t0)(i,0) * dt + (*v_tm1)(i,0);
        //         (*v_tm1)(i,0) = (*v_t0)(i,0);
        //
        //         (*v_t0)(i,1) = (*dv_dt_t0)(i,1) * dt + (*v_tm1)(i,1);
        //         (*v_tm1)(i,1) = (*v_t0)(i,1);
        //
        //     }
        //
        //     updateEnergy(globals, e_diss, *energy_diss, *v_t0, *cv_mass);
        //
        //     total_diss[0] = e_diss;
        //
        //
        // }

        iter ++;

        // if (iter >= 3)  globals->Output->TerminateODIS();

        // Check for output
        if (out_time >= out_frac*orbit_period)
        {
            outstring << std::fixed <<"DUMPING DATA AT "<<current_time/orbit_period;
            outstring << " AVG DISS: "<<std::scientific<<*total_diss*4*pi*r*r/1e9<<" GW"<<out_count;

            Output->Write(OUT_MESSAGE, &outstring);

            out_time -= out_frac*orbit_period;

            Output->DumpData(globals, out_count, pp);
            out_count++;
        }

        if (flag) {
            outstring << "Terminate signal caught..."<<std::endl;
            Output->Write(OUT_MESSAGE, &outstring);

            break;
        }
    }

    // Save initial conditions!
    writeInitialConditions(globals, grid, *v_t0, *dv_dt, *p_t0, *dp_dt);

    Output->Write(OUT_MESSAGE, &outstring);

    if (flag) return 1;     // return with error
    else return 0;          // return without error
};
