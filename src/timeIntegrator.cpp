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
#include "momAdvection.h"
#include "gridConstants.h"
#include "updateMomentum.h"
#include "updateEta.h"
#include <math.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <signal.h>


//#include <mkl.h>
//#include <omp.h>

#include <sys/time.h>

int flag = 0;

void CatchExit(int sig) {
    printf("%s\n", "Caught Terminate Signal...");
    flag = 1;
}

int IMPLICIT = 0;

double GAMMA = 0.5;
double BETA = 0.5;

typedef unsigned long long timestamp_t;
double total_time=0.0;
int time_iter=0;

static timestamp_t get_timestamp ()
{
  struct timeval now;
  gettimeofday (&now, NULL);
  return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

// Function to implement the time loop to solve mass and momentum equations over
// the geodesic grid using the Adams-Bashforth 3rd order forward method.
int ab3Explicit(Globals * globals, Mesh * grid)
{
    //--------------------------------------------------------------------------
    //------------------------- DECLARE VARIABLES ------------------------------
    //--------------------------------------------------------------------------

    // omp_set_num_threads( globals->core_num.Value() );
    // mkl_set_num_threads( globals->core_num.Value() );

    std::ostringstream outstring;
    int i, j, k, iter, out_count;
    double ** pp;

    Array1D<double> v_t0    (grid->face_num);     // velocity array
    Eigen::Map<Eigen::VectorXd> v_t0_eig(&v_t0(0), grid->face_num_ng); // Misses out the ghost nodes

    Array1D<double> v_guess    (FACE_NUM);     // velocity array
    Eigen::Map<Eigen::VectorXd> v_guess_eig(&v_guess(0), FACE_NUM);

    Array1D<double> drag_term    (grid->face_num_ng);     // velocity array
    Eigen::Map<Eigen::VectorXd> drag_term_eig(&drag_term(0), grid->face_num_ng);
    Array1D<double> dv_dt_t0(grid->face_num_ng);     // velocity rate of change array

    Array2D<double> v_xyz   (grid->node_num_ng, 3);  // Cartesian velocity components 

    Array2D<double> v_avg(grid->face_num_ng, 2);     // velocity components array

    Array1D<double> p_t0    (grid->node_num);     // pressure array
    Eigen::Map<Eigen::VectorXd> p_t0_eig(&p_t0(0), grid->node_num);
    Array1D<double> p_guess    (NODE_NUM);     // pressure array
    Eigen::Map<Eigen::VectorXd> p_guess_eig(&p_guess(0), NODE_NUM);
    Array1D<double> dp_dt_t0(grid->node_num_ng);     // pressure rate of change array

    Array2D<double> dv_dt(grid->face_num_ng, 3);     // [face][dvdt level]
    Array2D<double> dp_dt(grid->node_num_ng, 3);     // [cv][dpdt level]

    Array1D<double> cv_mass     (NODE_NUM); // control volume mass array
    Array1D<double> energy_diss (FACE_NUM); // energy dissipation array

    Array1D<double> dummy_vel(FACE_NUM);

    Array1D<double> ekin(NODE_NUM);
    Array1D<double> h_total(NODE_NUM);

    Array1D<double> forcing_potential(grid->node_num);
    Eigen::Map<Eigen::VectorXd> forcing_potential_eig(&forcing_potential(0), grid->node_num);
    // Array1D<double> 

#ifdef TESTS

    Array1D<double> test_v(FACE_NUM);
    Array2D<double> test_dvdt(FACE_NUM, 3);
    Array1D<double> test_p(NODE_NUM);
    Array2D<double> test_dpdt(NODE_NUM, 3);

#endif

    double end_time, current_time, dt, orbit_period;

    double r = globals->radius.Value();

    OutFiles * Output;

    signal(SIGINT, CatchExit);

    //--------------------------------------------------------------------------
    //------------------------- ASSIGN VARIABLES -------------------------------
    //--------------------------------------------------------------------------

    end_time     = globals->endTime.Value();
    dt           = globals->timeStep.Value();
    orbit_period = globals->period.Value();

    Output = globals->Output;
    pp = new double *[globals->out_tags.size()];

    std::vector<std::string> * tags;
    tags = &globals->out_tags;

    double * total_diss;
    total_diss = new double[1];
    double e_diss;

    outstring << "Defining arrays for Adams-Bashforth time integration..." << std::endl;


    for (i=0; i<globals->out_tags.size(); i++)
    {
        if      ((*tags)[i] == "velocity output")        pp[i] = &v_avg(0,0);
        else if ((*tags)[i] == "velocity cartesian output")   pp[i] = &v_xyz(0,0);
        else if ((*tags)[i] == "displacement output")    pp[i] = &p_t0(0);
        else if ((*tags)[i] == "dissipation output")     pp[i] = &energy_diss(0);
        else if ((*tags)[i] == "dissipation avg output") pp[i] = &total_diss[0];
        else if ((*tags)[i] == "kinetic avg output")     pp[i] = &current_time;
        else if ((*tags)[i] == "dummy1 output")          pp[i] = &cv_mass(0);
    }

    Eigen::Map<Eigen::VectorXd> div(&cv_mass(0), NODE_NUM);
    Eigen::Map<Eigen::VectorXd> velocity(&dummy_vel(0), FACE_NUM);


    // for (int i=0; i<FACE_NUM; ++i)
    //     {
    //         int node_in = grid->face_nodes(i, 0);
    //         int node_out = grid->face_nodes(i, 1);

    //         dummy_vel(i) = v_t0(i) * ((p_t0(node_in) + p_t0(node_out))*0.5);    // * avg_thickness - 0.5*fabs(v_t0(i))*diff_thickness;///thickness*v_t0(i);
    //     }   

    // div = grid->operatorDivergence*velocity;

    getInitialConditions(globals, grid, v_t0, dv_dt, p_t0, dp_dt);

    iter = 0;

    interpolateVelocity(globals, grid, v_avg, v_t0);
    interpolateVelocityCartRBF(globals, grid, v_xyz, v_t0);

    

    // updateEnergy(globals, e_diss, energy_diss, v_avg, grid->face_area);
    total_diss[0] = e_diss;

    

    //--------------------------------------------------------------------------
    //------------------------- BEGIN TIME LOOP --------------------------------
    //--------------------------------------------------------------------------

    out_count = 1;
    current_time = dt*iter;
    int out_freq = globals->totalIter.Value() / globals->outputTime.Value();

    outstring << std::fixed << "DUMPING DATA AT " << current_time / orbit_period;
    outstring << " AVG DISS: " << std::scientific << *total_diss * 4 * pi * r * r / 1e9 << " GW" << out_count;

    Output->Write(OUT_MESSAGE, &outstring);

    Output->DumpData(globals, out_count, pp);
    out_count++;

    double h = globals->h.Value();
    double alpha = globals->alpha.Value();
    double gr = 1.0/globals->g.Value();
    double g = globals->g.Value();

    // for (int i=0; i < NODE_NUM; i++) h_total(i) = h + p_t0(i);

    while (iter < globals->totalIter.Value()*globals->endTime.Value())
    {
        if (globals->surface_type == FREE ||
            globals->surface_type == FREE_LOADING ||
            globals->surface_type == LID_LOVE ||
            globals->surface_type == LID_MEMBR)
        {
            std::cout<<"HERE"<<std::endl;
            // SOLVE THE MOMENTUM EQUATION
            updateMomentum(globals, grid, dv_dt_t0, v_t0, p_t0, h_total, ekin, GAMMA, IMPLICIT);

            std::cout<<"updated momentum!"<<std::endl;

            for (i=0; i<grid->face_num_ng; ++i) dv_dt(i, 0) = dv_dt_t0(i);

            // Apply drag and forcing
            forcing(globals, grid, forcing_potential, globals->tide_type, current_time+dt, globals->e.Value(), globals->theta.Value());
            drag_term_eig = grid->operatorLinearDrag * v_t0_eig + grid->operatorGradient * forcing_potential_eig; 
            

            std::cout<<"updated drag!"<<std::endl;
            // // Implicit approach - does not work with advection yet!!
            // if (IMPLICIT) {
            //      v_guess_eig = v_t0_eig;// - dt*g*(1.0 - GAMMA)*grid->operatorGradient*p_t0_eig;
            //     integrateAB3scalar(globals, grid, v_guess, dv_dt, iter, FACE_NUM);

            //     // Now remove drag, as it is outside of the AB3 loop
            //     v_guess_eig += dt*drag_term_eig;  
            //     p_guess_eig = p_t0_eig - -h*dt*grid->operatorDivergence*(BETA*v_guess_eig+(1.0-BETA)*v_t0_eig);

            //     // p_t0_eig = grid->cg.solveWithGuess(p_guess_eig, p_t0_eig);
            //     p_t0_eig = grid->cg.solve(p_guess_eig);

            //     v_t0_eig = v_guess_eig - GAMMA*dt*g*grid->operatorGradient*p_t0_eig;

            // }
            // else {


            integrateAB3scalar(globals, grid, v_t0, dv_dt, iter, grid->face_num_ng);

            std::cout<<"updated v!"<<std::endl;

            // Now remove drag, as it is outside of the AB3 loop
            v_t0_eig += dt*drag_term_eig;

            std::cout<<"added drag!"<<std::endl;

#if defined(TEST_SW1) || defined(TEST_GAUSS_HILLS) 
            // Overwrite velocity field                 dummy    dummy     dummy 
            getInitialConditions(globals, grid, v_t0, test_dvdt, test_p, test_dpdt, current_time);
#endif

            updateEta(globals, grid, dp_dt_t0, v_t0, p_t0, h_total);

            std::cout<<"updated continuity!"<<std::endl;

            for (i=0; i<grid->node_num_ng; ++i) dp_dt(i, 0) = dp_dt_t0(i);

            integrateAB3scalar(globals, grid, p_t0, dp_dt, iter, grid->node_num_ng);

            // }


           

            // UPDATE ENERGY DISSIPATION
            interpolateVelocity(globals, grid, v_avg, v_t0);
            // updateEnergy(globals, e_diss, energy_diss, v_avg, grid->face_area);
            total_diss[0] = e_diss;

            // UPDATE COLUMN THICKNESS
            if (globals->advection.Value())
            {
                for (int i=0; i < NODE_NUM; i++) h_total(i) = h + p_t0(i);
            }
            

        }

        

        iter ++;
        current_time = dt*iter;

        // Check for output
        if (true)//(iter%out_freq == 0)
        {
            // for (int i=0; i<FACE_NUM; ++i)
            // {
            //     int node_in = grid->face_nodes(i, 0);
            //     int node_out = grid->face_nodes(i, 1);

            //     dummy_vel(i) = v_t0(i) * ((p_t0(node_in) + p_t0(node_out))*0.5);    // * avg_thickness - 0.5*fabs(v_t0(i))*diff_thickness;///thickness*v_t0(i);
            // }   

            interpolateVelocityCartRBF(globals, grid, v_xyz, v_t0);

            // std::cout<<v_xyz(15, 0)<<std::endl;

            // div = grid->operatorDivergence*velocity;

            outstring << std::fixed <<"DUMPING DATA AT "<<current_time/orbit_period;
            outstring << " AVG DISS: "<<std::scientific<<*total_diss*4*pi*r*r/1e9<<" GW"<<out_count;

            Output->Write(OUT_MESSAGE, &outstring);

            Output->DumpData(globals, out_count, pp);
            out_count++;

            Output->TerminateODIS();
            
        }


        if (flag) {
            outstring << "Terminate signal caught..."<<std::endl;
            Output->Write(OUT_MESSAGE, &outstring);

            break;
        }
    }

    // Save initial conditions!
    // writeInitialConditions(globals, grid, v_t0, dv_dt, p_t0, dp_dt);

    Output->Write(OUT_MESSAGE, &outstring);

    if (flag) return 1;     // return with error
    else return 0;          // return without error
};


// Function to implement the time loop to solve mass and momentum equations over
// the geodesic grid using the Runge-Kutta 4th order forward method.
// int rk4Explicit(Globals * globals, Mesh * grid)
// {
//     //--------------------------------------------------------------------------
//     //------------------------- DECLARE VARIABLES ------------------------------
//     //--------------------------------------------------------------------------

//     // omp_set_num_threads( globals->core_num.Value() );
//     // mkl_set_num_threads( globals->core_num.Value() );

//     std::ostringstream outstring;
//     int i, j, k, iter, out_count;
//     double ** pp;

//     Array1D<double> v_t0    (FACE_NUM);     // velocity array
//     Eigen::Map<Eigen::VectorXd> v_t0_eig(&v_t0(0), FACE_NUM);

//     Array1D<double> v_guess    (FACE_NUM);     // velocity array
//     Eigen::Map<Eigen::VectorXd> v_guess_eig(&v_guess(0), FACE_NUM);

//     Array1D<double> drag_term    (FACE_NUM);     // velocity array
//     Eigen::Map<Eigen::VectorXd> drag_term_eig(&drag_term(0), FACE_NUM);
//     Array1D<double> dv_dt_t0(FACE_NUM);     // velocity rate of change array

//     Array2D<double> v_xyz   (NODE_NUM, 3);  // Cartesian velocity components 

//     Array2D<double> v_avg(FACE_NUM, 2);     // velocity components array

//     Array1D<double> p_t0    (NODE_NUM);     // pressure array
//     Eigen::Map<Eigen::VectorXd> p_t0_eig(&p_t0(0), NODE_NUM);
//     Array1D<double> p_guess    (NODE_NUM);     // pressure array
//     Eigen::Map<Eigen::VectorXd> p_guess_eig(&p_guess(0), NODE_NUM);
//     Array1D<double> dp_dt_t0(NODE_NUM);     // pressure rate of change array

//     Array2D<double> dv_dt(FACE_NUM, 3);     // [face][dvdt level]
//     Array2D<double> dp_dt(NODE_NUM, 3);     // [cv][dpdt level]

//     Array1D<double> cv_mass     (NODE_NUM); // control volume mass array
//     Array1D<double> energy_diss (FACE_NUM); // energy dissipation array

//     Array1D<double> dummy_vel(FACE_NUM);


//     Array1D<double> v_t0_old(FACE_NUM);
//         // Array1D<double> dv_dt_t0(face_num);

//         // Array1D<double> p_t0    (node_num);
//     Array1D<double> p_t0_old(NODE_NUM);
//     // Array1D<double> dp_dt_t0(node_num);

//     Array2D<double> dv_dt4(FACE_NUM, 4);   //[node][vel component][dvdt]
//     Array2D<double> dp_dt4(NODE_NUM, 4); 


//     double end_time, current_time, dt, out_frac, orbit_period, out_time;

//     double r = globals->radius.Value();

//     int node_num;
//     int face_num;
//     node_num = globals->node_num;
//     face_num = globals->face_num;

//     OutFiles * Output;

//     signal(SIGINT, CatchExit);

//     //--------------------------------------------------------------------------
//     //------------------------- ASSIGN VARIABLES -------------------------------
//     //--------------------------------------------------------------------------

//     end_time     = globals->endTime.Value();
//     dt           = globals->timeStep.Value();
//     out_frac     = globals->outputTime.Value();
//     orbit_period = globals->period.Value();

//     Output = globals->Output;
//     pp = new double *[globals->out_tags.size()];

//     std::vector<std::string> * tags;
//     tags = &globals->out_tags;

//     double * total_diss;
//     total_diss = new double[1];
//     double e_diss;

//     outstring << "Defining arrays for Runge-Kutta time integration..." << std::endl;

//     // Array1D<double> v_t0    (face_num);
//     // Array2D<double> v_avg   (face_num, 2);
//     //[node][dpdt]

//     // Array1D<double> cv_mass    (node_num);
//     // Array1D<double> energy_diss(face_num);

//     for (i=0; i<globals->out_tags.size(); i++)
//     {
//         if      ((*tags)[i] == "velocity output")        pp[i] = &v_avg(0,0);
//         else if ((*tags)[i] == "velocity cartesian output")   pp[i] = &v_xyz(0,0);
//         else if ((*tags)[i] == "displacement output")    pp[i] = &p_t0(0);
//         else if ((*tags)[i] == "dissipation output")     pp[i] = &energy_diss(0);
//         else if ((*tags)[i] == "dissipation avg output") pp[i] = &total_diss[0];
//         else if ((*tags)[i] == "kinetic avg output")     pp[i] = &current_time;
//         else if ((*tags)[i] == "dummy1 output")          pp[i] = &cv_mass(0);
//     }

//     for (i=0; i<face_num; i++)
//         v_t0(i) = 0.00;

//     // for (i=0; i<node_num; i++)
//     // {
//     //     double lat=0*pi/180;
//     //     double lon=150*pi/180;
//     //     double llat=30*pi/180.0;
//     //     double dist, lat1, lon1;
//     //     double rr=1.0;
//     //     dist = 0.0;
//     //     double dist2 = 0.0;
//     //     lat1 = grid->node_pos_sph(i,0);
//     //     lon1 = grid->node_pos_sph(i,1);

//     //     distanceBetweenSph(dist, lat, lat1, lon, lon1, rr);
//     //     // distanceBetweenSph(dist2, llat, lat1, lon, lon1, rr);

//     //     // p_t0(i) = 0.0;//-200./(1.+2*dist+3*dist2);
//     //     // (*v_t0)(i,1) = 0.0;


//     //     cv_mass(i) = grid->control_volume_mass(i);
//     // }

   
//     getInitialConditions(globals, grid, v_t0, dv_dt, p_t0, dp_dt);


//     interpolateVelocity(globals, grid, v_avg, v_t0);
//     interpolateVelocityCartRBF(globals, grid, v_xyz, v_t0);

//     //--------------------------------------------------------------------------
//     //------------------------- BEGIN TIME LOOP --------------------------------
//     //--------------------------------------------------------------------------

//     out_count = 1;
//     current_time = 0.0;
//     out_time = 0.0;

//     out_count = 1;
//     current_time = dt*iter;
//     int out_freq = globals->totalIter.Value() / globals->outputTime.Value();

//     outstring << std::fixed << "DUMPING DATA AT " << current_time / orbit_period;
//     outstring << " AVG DISS: " << std::scientific << *total_diss * 4 * pi * r * r / 1e9 << " GW" << out_count;

//     Output->Write(OUT_MESSAGE, &outstring);

//     Output->DumpData(globals, out_count, pp);
//     out_count++;

//     double alpha = globals->alpha.Value();
//     double h = globals->h.Value();
    

//     while (iter < globals->totalIter.Value()*globals->endTime.Value())
//     {
//         if (globals->surface_type == FREE ||
//             globals->surface_type == FREE_LOADING ||
//             globals->surface_type == LID_LOVE ||
//             globals->surface_type == LID_MEMBR)
//         {
//             for (i=0; i<NODE_NUM; ++i)
//               p_t0_old(i) = p_t0(i);

//             for (i=0; i<FACE_NUM; ++i)
//               v_t0_old(i) = v_t0(i);


//             for (int j=0; j<4; j++)
//             {
//               double x = 1.0;
//               if (j>0 && j<3) x = 0.5;

//               // SOLVE THE MOMENTUM EQUATION
//               updateVelocity(globals, grid, dv_dt_t0, v_t0, p_t0, current_time + x*dt);

//               if (globals->fric_type == QUADRATIC)
//                   quadraticDrag(face_num, alpha, h, dv_dt_t0, v_avg, v_t0);

//               // MARCH VELOCITY FORWARD IN TIME
//               #pragma omp parallel for
//               for (i=0; i<FACE_NUM; ++i) {
//                 dv_dt4(i, j) = dv_dt_t0(i);
//                 v_t0(i) = v_t0_old(i) + dv_dt_t0(i)*dt*x;
//               }

//               // SOLVE THE CONTINUITY EQUATION
//               updateDisplacement(globals, grid, dp_dt_t0, v_t0, p_t0);

//               // MARCH PRESSURE FORWARD IN TIME
//               #pragma omp parallel for
//               for (i=0; i<NODE_NUM; ++i) {
//                 dp_dt4(i, j) = dp_dt_t0(i);
//                 p_t0(i)     = p_t0_old(i) + dp_dt_t0(i)*dt*x;
//               }
//             }

//             #pragma omp parallel for
//             for (i=0; i<NODE_NUM; ++i)
//               p_t0(i) = p_t0_old(i) + dt/6.*(dp_dt4(i,0)+2.*dp_dt4(i,1)+2.*dp_dt4(i,2)+dp_dt4(i,3));

//             #pragma omp parallel for
//             for (i=0; i<FACE_NUM; ++i)
//               v_t0(i) = v_t0_old(i) + dt/6.*(dv_dt4(i,0)+2.*dv_dt4(i,1)+2.*dv_dt4(i,2)+dv_dt4(i,3));


//             // UPDATE ENERGY DISSIPATION
//             interpolateVelocity(globals, grid, v_avg, v_t0);
//             updateEnergy(globals, e_diss, energy_diss, v_avg, grid->face_area);
//             total_diss[0] = e_diss;
//         }

//         iter ++;
//         current_time = dt*iter;

//         // if (iter >= 3)  globals->Output->TerminateODIS();

//         // Check for output
//         if (iter%out_freq == 0)
//         // if (true)
//         {

//             interpolateVelocityCartRBF(globals, grid, v_xyz, v_t0);
//             outstring << std::fixed <<"DUMPING DATA AT "<<current_time/orbit_period;
//             outstring << " AVG DISS: "<<std::scientific<<*total_diss*4*pi*r*r/1e9<<" GW"<<out_count;

//             Output->Write(OUT_MESSAGE, &outstring);

//             out_time -= out_frac*orbit_period;

//             Output->DumpData(globals, out_count, pp);
//             out_count++;
//         }

//         if (flag) {
//             outstring << "Terminate signal caught..."<<std::endl;
//             Output->Write(OUT_MESSAGE, &outstring);

//             break;
//         }
//     }

//     // Save initial conditions!
//     writeInitialConditions(globals, grid, v_t0, dv_dt, p_t0, dp_dt);

//     Output->Write(OUT_MESSAGE, &outstring);

//     if (flag) return 1;     // return with error
//     else return 0;          // return without error
// };


// int updateVelocity(Globals * globals, Mesh * grid, Array1D<double> & dvdt, Array1D<double> & v_tm1, Array1D<double> & p_tm1, double current_time)
// {
//     double r, omega, e, obliq, h, drag_coeff, visc, g, gr, dt;
//     int node_num;

//     r =       globals->radius.Value();
//     g =       globals->g.Value();
//     gr = 1.0/g;
//     dt =      globals->timeStep.Value();
//     omega =   globals->angVel.Value();
//     e =       globals->e.Value();
//     drag_coeff = globals->alpha.Value();
//     obliq =   globals->theta.Value();
//     h =       globals->h.Value();
//     node_num = globals->node_num;
//     int face_num = globals->face_num;


//     visc = 5e3;

//     Array1D<double> soln_new(NODE_NUM);
//     Array1D<double> total_thickness(NODE_NUM);

//     Array1D<double> forcing_potential(NODE_NUM);
//     Eigen::Map<Eigen::VectorXd> forcing_potential_eig(&forcing_potential(0), NODE_NUM);

//     Array1D<double> advection_term(FACE_NUM);
//     Eigen::Map<Eigen::VectorXd> adv_eig(&advection_term(0), FACE_NUM);
//     Array1D<double> Ekin(NODE_NUM);

    
//     forcing(globals, grid, forcing_potential, globals->tide_type, current_time, e, obliq);
    

//     // Array1D<double> self_gravity(node_num);
// /*
//     // TODO - Move these pressure grad calls to a function that makes sense
//     switch (globals->surface_type)
//     {
//         case FREE_LOADING:
//             pressureGradientSH(globals, grid, dvdt, p_tm1, forcing_potential, g);
//             break;

//         case LID_LOVE:
//             pressureGradientSH(globals, grid, dvdt, p_tm1, forcing_potential, g);
//             break;

//         case LID_MEMBR:
//             pressureGradientSH(globals, grid, dvdt, p_tm1, forcing_potential, g);
//             break;

//         case LID_INF:
//             // pressureGradient(grid, dvdt, p_tm1, 1.0/1000.0);
//             // pressureGradientSH(globals, grid, dvdt, p_tm1, -1.0/1000.0);
//             break;
//     }
// */

//     for (int i=0; i<NODE_NUM; ++i)
//     {
//         soln_new(i) = (1-GAMMA*(double)IMPLICIT)*p_tm1(i);// - forcing_potential(i)*gr;
//         total_thickness(i) = h + p_tm1(i);

//     }

//     Eigen::Map<Eigen::VectorXd> solnEigen(&soln_new(0), NODE_NUM);
//     Eigen::Map<Eigen::VectorXd> dvdtEigen(&dvdt(0), FACE_NUM);
//     Eigen::Map<Eigen::VectorXd> velEigen(&v_tm1(0), FACE_NUM);

//     // Perform sparse matrix * vector operation
//     // dvdtEigen.setZero();
//     dvdtEigen = -g*grid->operatorGradient * solnEigen;// + grid->operatorLinearDrag * velEigen;

//     // dvdtEigen += grid->operatorLinearDrag * velEigen;

//     if (globals->advection.Value()) 
//     {
//         calculateMomentumAdvection(globals, grid, v_tm1, total_thickness, advection_term, Ekin);
        
//         dvdtEigen -= adv_eig;
//     }
//     else {
//         dvdtEigen += grid->operatorCoriolis * velEigen;
//     }
    
//     return 1;

// };


// int updateDisplacement(Globals * globals, Mesh * grid, Array1D<double> & deta_dt, Array1D<double> & v_t0, Array1D<double> & eta, Array1D<double> &h_total)
// {
//     double sum = 0.0;
//     int node_num = globals->node_num;
//     int face_num = globals->face_num;
//     double h = globals->h.Value();
//     int advection = globals->advection.Value();

//     double alpham = 1.0;
//     double betam = 0.0;

    

//     Eigen::Map<Eigen::VectorXd> detadtEigen(&deta_dt(0), NODE_NUM);

    

//     if (globals->advection.Value()) 
//     {
//         Array1D<double> hv(FACE_NUM);
//         Eigen::Map<Eigen::VectorXd> hvEigen(&hv(0), FACE_NUM);


//         // Get thickness flux normal to each face
//         interpolateLSQFlux(globals, grid, hv, v_t0, h_total);

//         // Get divergence of thickness flux to find dEta/dt
//         detadtEigen = grid->operatorDivergence * hvEigen;
//     }
//     else 
//     {
//         Eigen::Map<Eigen::VectorXd> hvEigen(&v_t0(0), FACE_NUM);

//         // Get divergence of velocity field to find dEta/dt
//         detadtEigen = h * grid->operatorDivergence * hvEigen;   
//     }
     

//     return 1;
// };