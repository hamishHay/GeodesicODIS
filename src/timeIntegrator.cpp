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
#include "interpolation.h"
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

int updateVelocity(Globals * globals, Mesh * grid, Array1D<double> & dvdt, Array1D<double> & v_tm1, Array1D<double> & p_tm1, double current_time)
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
    int face_num = globals->face_num;

    // std::cout<<omega<<' '<<r<<' '<<g<<' '<<h<<' '<<globals->period.Value()<<std::endl;

    visc = 5e3;

    Array1D<double> forcing_potential(node_num);

    forcing(globals, grid, forcing_potential, globals->tide_type, current_time, e, obliq);

    // Array1D<double> self_gravity(node_num);

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
    //
    // Array2D<double> soln_old(node_num, 3);
    Array1D<double> soln_new(node_num);

    for (int i=0; i<face_num; ++i)
    {
      dvdt(i) = 0.0;
      // dvdt(i, 1) = 0.0;
    }

    // #pragma omp parallel for
    for (int i=0; i<node_num; ++i)
    {
        // soln_old(i, 0) = v_tm1(i, 0);
        // soln_old(i, 1) = v_tm1(i, 1);
        // soln_old(i, 2) = p_tm1(i);// - forcing_potential(i)/g;
        soln_new(i) = p_tm1(i) - forcing_potential(i)/g;
    }

    mkl_sparse_d_mv(operation, 1.0, *(grid->operatorGradient), descript, &(soln_new(0)), betam, &(dvdt(0)));
    // pressureGradientN(grid, dvdt, soln_new, node_num, g);

    // mkl_sparse_d_mv(operation, 1.0, *(grid->operatorLinearDrag), descript, &(v_tm1(0)), 1.0, &(dvdt(0)));
    // // linearDrag(node_num, drag_coeff, dvdt, v_tm1);
    //
    // mkl_sparse_d_mv(operation, 1.0, *(grid->operatorCoriolis), descript, &(v_tm1(0)), 1.0, &(dvdt(0)));
    // timestamp_t t0 = get_timestamp();

    // Perform A*d = dvdt where A is operatorMomentum, and d = [u, v, eta-U/g]
    mkl_sparse_d_mv(operation, 1.0, *(grid->operatorMomentum), descript, &(v_tm1(0)), 1.0, &(dvdt(0)));


    // timestamp_t t1 = get_timestamp();
    //
    // total_time += (t1 - t0) / 1000000.0L;
    //
    // time_iter += 1;

    // std::cout<<"TIME: "<<total_time/((double)time_iter)<<std::endl;

    // if (globals->fric_type == QUADRATIC) quadraticDrag(node_num, drag_coeff, h, dvdt, v_tm1);

    return 1;

};

int updateDisplacement(Globals * globals, Mesh * grid, Array1D<double> & deta_dt, Array1D<double> & v_t0, Array1D<double> & eta)
{
    double sum = 0.0;
    int node_num = globals->node_num;
    int face_num = globals->face_num;
    double h = globals->h.Value();

    sparse_operation_t operation = SPARSE_OPERATION_NON_TRANSPOSE;
    matrix_descr descript;
    descript.type = SPARSE_MATRIX_TYPE_GENERAL;

    double alpham = 1.0;
    double betam = 0.0;

    mkl_sparse_d_mv(operation, globals->h.Value(), *(grid->operatorDivergenceX), descript, &(v_t0(0)), 0.0, &(deta_dt(0)));
    // Array1D<double> hv(globals->face_num);

    // for (int i=0; i<face_num; ++i)
    // {
    //     int node_in = grid->face_nodes(i, 0);
    //     int node_out = grid->face_nodes(i, 1);
    //     double inner = eta(node_in);
    //     double outer = eta(node_out);
    //     double diff_thickness = outer - inner;
    //     double avg_thickness = 0.5 * ( (inner + h) + (outer + h) );

    //     // check if flow is in or out of face
    //     // if outward, use inner h 
    //     // if inward, use outer h

    //     // if (v_t0(i) > 0.0) thickness = inner+h;
    //     // else thickness = outer+h;

    //     hv(i) = v_t0(i) * avg_thickness - 0.5*fabs(v_t0(i))*diff_thickness;///thickness*v_t0(i);
    //     // hv(i,1) = (eta(i)+globals->h.Value())*v_t0(i);
    // }

    // mkl_sparse_d_mv(operation, 1.0, *(grid->operatorDivergenceX), descript, &(hv(0)), 0.0, &(deta_dt(0)));

    // velocityDivergenceN(grid, deta_dt, v_t0, sum, globals->h.Value());
    // velocityDivergence(grid, deta_dt, v_t0, sum, globals->h.Value());

    // for (int i=0; i<node_num; ++i)
    // {
    //     if (grid->cell_is_boundary(i) == 2) deta_dt(i) = 0.0;
    // }

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

    // Array1D<double> * v_t0;        // velocity solution for current timestep
    // Array1D<double> * dv_dt_t0;    // velocity time derivative at current timestep


    // Array1D<double> * p_t0;      // displacement solution for current timestep
    // Array1D<double> * dp_dt_t0;  // displacement time derivative at current timestep

    // Array2D<double> * dv_dt;
    // Array2D<double> * dp_dt;

    // Array1D<double> * energy_diss;
    // Array1D<double> * cv_mass;

    double end_time, current_time, dt, out_frac, orbit_period, out_time;

    double r = globals->radius.Value();

    int node_num;
    node_num = globals->node_num;
    int face_num = globals->face_num;

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

    Array1D<double> v_t0    (face_num);     // velocity array
    Array1D<double> dv_dt_t0(face_num);     // velocity rate of change array

    Array2D<double> v_avg(face_num, 2);     // velocity components array

    Array1D<double> p_t0    (node_num);     // pressure array
    Array1D<double> dp_dt_t0(node_num);     // pressure rate of change array

    Array2D<double> dv_dt(face_num, 3);     // [face][dvdt level]
    Array2D<double> dp_dt(node_num, 3);     // [cv][dpdt level]

    Array1D<double> cv_mass     (node_num); // control volume mass array
    Array1D<double> energy_diss (face_num); // energy dissipation array

    for (i=0; i<globals->out_tags.size(); i++)
    {
        if      ((*tags)[i] == "velocity output")        pp[i] = &v_avg(0,0);
        else if ((*tags)[i] == "displacement output")    pp[i] = &p_t0(0);
        else if ((*tags)[i] == "dissipation output")     pp[i] = &energy_diss(0);
        else if ((*tags)[i] == "dissipation avg output") pp[i] = &total_diss[0];
        else if ((*tags)[i] == "kinetic avg output")     pp[i] = &current_time;
    }

    for (i=0; i<face_num; i++) {
        v_t0(i) = 0.0;
        // (*v_t0)(i) = 10.0*grid->face_normal_vec_map(i,0);
    }

    for (i=0; i<node_num; i++)
    {
        v_avg(i,0) = 0.0;
        v_avg(i,1) = 0.0;
        double lat=0*pi/180;
        double lon=150*pi/180;
        double llat=30*pi/180.0;
        double dist, lat1, lon1;
        double llon = 0.0;
        double rr=1.0;
        dist = 0.0;
        double dist2 = 0.0;
        lat1 = grid->node_pos_sph(i,0);
        lon1 = grid->node_pos_sph(i,1);

        distanceBetweenSph(dist, lat, lat1, lon, lon1, rr);
        distanceBetweenSph(dist2, llat, lat1, lon, lon1, rr);

        // p_t0(i) = 400./(1.+dist*dist+3*dist2);

        p_t0(i) = 0.0;

        cv_mass(i) = grid->control_volume_mass(i);
    }



    if (globals->init.Value())
    {
       loadInitialConditions(globals, grid, v_t0, dv_dt, p_t0, dp_dt);

       iter = 3;
    }
    else iter = 0;

    interpolateVelocity(globals, grid, v_avg, v_t0);

    //--------------------------------------------------------------------------
    //------------------------- BEGIN TIME LOOP --------------------------------
    //--------------------------------------------------------------------------

    out_count = 1;
    current_time = 0.0;
    out_time = 0.0;

    double h = globals->h.Value();
    double alpha = globals->alpha.Value();

    while (current_time <= end_time)
    {
        if (globals->surface_type == FREE ||
            globals->surface_type == FREE_LOADING ||
            globals->surface_type == LID_LOVE ||
            globals->surface_type == LID_MEMBR)
        {
            // SOLVE THE MOMENTUM EQUATION
            updateVelocity(globals, grid, dv_dt_t0, v_t0, p_t0, current_time);

            if (globals->fric_type == QUADRATIC)
                quadraticDrag(face_num, alpha, h, dv_dt_t0, v_avg, v_t0);

            #pragma omp parallel for
            for (i=0; i<face_num; ++i) dv_dt(i, 0) = dv_dt_t0(i);

            // MARCH VELOCITY FORWARD IN TIME
            integrateAB3scalar(globals, grid, v_t0, dv_dt, iter, grid->face_num);


            for (i=0; i<face_num; ++i)
            {
                int inner, outer;
                inner = grid->face_nodes(i, 0);
                outer = grid->face_nodes(i, 1);
                if ((grid->cell_is_boundary(inner) == 1 && grid->cell_is_boundary(outer) == 0) ||
                    (grid->cell_is_boundary(inner) == 0 && grid->cell_is_boundary(outer) == 1))
                {
                    v_t0(i) = 0.0;
                }
                else if (grid->cell_is_boundary(inner) && grid->cell_is_boundary(outer))
                {
                    v_t0(i) = 0.0;
                }
            }

            // UPDATE ENERGY DISSIPATION
            interpolateVelocity(globals, grid, v_avg, v_t0);
            updateEnergy(globals, e_diss, energy_diss, v_avg, grid->face_area);
            total_diss[0] = e_diss;

            // SOLVE THE CONTINUITY EQUATION
            updateDisplacement(globals, grid, dp_dt_t0, v_t0, p_t0);

            #pragma omp parallel for
            for (i=0; i<node_num; ++i) dp_dt(i, 0) = dp_dt_t0(i);

            // MARCH PRESSURE FORWARD IN TIME
            integrateAB3scalar(globals, grid, p_t0, dp_dt, iter, grid->node_num);

            for (i=0; i<node_num; ++i)
            {
                if (grid->cell_is_boundary(i) == 1)
                {
                    p_t0(i) = 0.0;
                }
            }

        }

        current_time += dt;
        out_time += dt;

        iter ++;

        // Check for output
        if (out_time >= out_frac*orbit_period)
        {
            outstring << std::fixed <<"DUMPING DATA AT "<<current_time/orbit_period;
            outstring << " AVG DISS: "<<std::scientific<<*total_diss*4*pi*r*r/1e9<<" GW"<<out_count;

            Output->Write(OUT_MESSAGE, &outstring);

            out_time -= out_frac*orbit_period;

            Output->DumpData(globals, out_count, pp);
            out_count++;
            // Output->TerminateODIS();
        }


        if (flag) {
            outstring << "Terminate signal caught..."<<std::endl;
            Output->Write(OUT_MESSAGE, &outstring);

            break;
        }
    }

    // Save initial conditions!
    writeInitialConditions(globals, grid, v_t0, dv_dt, p_t0, dp_dt);

    Output->Write(OUT_MESSAGE, &outstring);

    if (flag) return 1;     // return with error
    else return 0;          // return without error
};


// Function to implement the time loop to solve mass and momentum equations over
// the geodesic grid using the Runge-Kutta 4th order forward method.
int rk4Explicit(Globals * globals, Mesh * grid)
{
    //--------------------------------------------------------------------------
    //------------------------- DECLARE VARIABLES ------------------------------
    //--------------------------------------------------------------------------

    omp_set_num_threads( globals->core_num.Value() );
    mkl_set_num_threads( globals->core_num.Value() );

    std::ostringstream outstring;
    int i, j, k, iter, out_count;
    double ** pp;

    double end_time, current_time, dt, out_frac, orbit_period, out_time;

    double r = globals->radius.Value();

    int node_num;
    int face_num;
    node_num = globals->node_num;
    face_num = globals->face_num;

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

    outstring << "Defining arrays for Runge-Kutta time integration..." << std::endl;

    Array1D<double> v_t0    (face_num);
    Array2D<double> v_avg   (face_num, 2);
    Array1D<double> v_t0_old(face_num);
    Array1D<double> dv_dt_t0(face_num);

    Array1D<double> p_t0    (node_num);
    Array1D<double> p_t0_old(node_num);
    Array1D<double> dp_dt_t0(node_num);

    Array2D<double> dv_dt(face_num, 4);   //[node][vel component][dvdt]
    Array2D<double> dp_dt(node_num, 4);      //[node][dpdt]

    Array1D<double> cv_mass    (node_num);
    Array1D<double> energy_diss(face_num);

    for (i=0; i<globals->out_tags.size(); i++)
    {
        if      ((*tags)[i] == "velocity output")        pp[i] = &v_avg(0,0);
        else if ((*tags)[i] == "displacement output")    pp[i] = &p_t0(0);
        else if ((*tags)[i] == "dissipation output")     pp[i] = &energy_diss(0);
        else if ((*tags)[i] == "dissipation avg output") pp[i] = &total_diss[0];
        else if ((*tags)[i] == "kinetic avg output")     pp[i] = &current_time;
    }

    for (i=0; i<face_num; i++)
        v_t0(i) = 0.00;

    for (i=0; i<node_num; i++)
    {
        double lat=0*pi/180;
        double lon=150*pi/180;
        double llat=30*pi/180.0;
        double dist, lat1, lon1;
        double rr=1.0;
        dist = 0.0;
        double dist2 = 0.0;
        lat1 = grid->node_pos_sph(i,0);
        lon1 = grid->node_pos_sph(i,1);

        distanceBetweenSph(dist, lat, lat1, lon, lon1, rr);
        // distanceBetweenSph(dist2, llat, lat1, lon, lon1, rr);

        p_t0(i) = 0.0;//-200./(1.+2*dist+3*dist2);
        // (*v_t0)(i,1) = 0.0;


        cv_mass(i) = grid->control_volume_mass(i);
    }

    if (globals->init.Value())
    {
       loadInitialConditions(globals, grid, v_t0, dv_dt, p_t0, dp_dt);
    
       iter = 0;
    }
    else iter = 0;

    //--------------------------------------------------------------------------
    //------------------------- BEGIN TIME LOOP --------------------------------
    //--------------------------------------------------------------------------

    out_count = 1;
    current_time = 0.0;
    out_time = 0.0;


    double alpha = globals->alpha.Value();
    double h = globals->h.Value();

    while (current_time <= end_time)
    {
        if (globals->surface_type == FREE ||
            globals->surface_type == FREE_LOADING ||
            globals->surface_type == LID_LOVE ||
            globals->surface_type == LID_MEMBR)
        {
            for (i=0; i<node_num; ++i)
              p_t0_old(i) = p_t0(i);

            for (i=0; i<face_num; ++i)
              v_t0_old(i) = v_t0(i);


            for (int j=0; j<4; j++)
            {
              double x = 1.0;
              if (j>0 && j<3) x = 0.5;

              // SOLVE THE MOMENTUM EQUATION
              updateVelocity(globals, grid, dv_dt_t0, v_t0, p_t0, current_time + x*dt);

              if (globals->fric_type == QUADRATIC)
                  quadraticDrag(face_num, alpha, h, dv_dt_t0, v_avg, v_t0);

              // MARCH VELOCITY FORWARD IN TIME
              #pragma omp parallel for
              for (i=0; i<face_num; ++i) {
                dv_dt(i, j) = dv_dt_t0(i);
                v_t0(i) = v_t0_old(i) + dv_dt_t0(i)*dt*x;
              }

              // SOLVE THE CONTINUITY EQUATION
              updateDisplacement(globals, grid, dp_dt_t0, v_t0, p_t0);

              // MARCH PRESSURE FORWARD IN TIME
              #pragma omp parallel for
              for (i=0; i<node_num; ++i) {
                dp_dt(i, j) = dp_dt_t0(i);
                p_t0(i)     = p_t0_old(i) + dp_dt_t0(i)*dt*x;
              }
            }

            #pragma omp parallel for
            for (i=0; i<node_num; ++i)
              p_t0(i) = p_t0_old(i) + dt/6.*(dp_dt(i,0)+2.*dp_dt(i,1)+2.*dp_dt(i,2)+dp_dt(i,3));

            #pragma omp parallel for
            for (i=0; i<face_num; ++i)
              v_t0(i) = v_t0_old(i) + dt/6.*(dv_dt(i,0)+2.*dv_dt(i,1)+2.*dv_dt(i,2)+dv_dt(i,3));


            // UPDATE ENERGY DISSIPATION
            interpolateVelocity(globals, grid, v_avg, v_t0);
            updateEnergy(globals, e_diss, energy_diss, v_avg, grid->face_area);
            total_diss[0] = e_diss;
        }

        current_time += dt;
        out_time += dt;

        iter ++;

        // if (iter >= 3)  globals->Output->TerminateODIS();

        // Check for output
        if (out_time >= out_frac*orbit_period)
        // if (true)
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
    writeInitialConditions(globals, grid, v_t0, dv_dt, p_t0, dp_dt);

    Output->Write(OUT_MESSAGE, &outstring);

    if (flag) return 1;     // return with error
    else return 0;          // return without error
};
