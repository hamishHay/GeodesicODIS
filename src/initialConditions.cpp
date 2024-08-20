#include "mesh.h"
#include "globals.h"
#include "gridConstants.h"
#include "array3d.h"
#include "array2d.h"
#include "array1d.h"
#include "analyticalLTE.h"
#include "mathRoutines.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <array>

#include <sys/stat.h>



int loadInitialConditions(Globals * globals, Mesh * mesh,
                          Array1D<double> & v, Array2D<double> & dvdt,
                          Array1D<double> & p, Array2D<double> & dpdt)
{
    int i, j;
    std::ostringstream outstring;

    std::string line, val;                 // strings for column and individual number
    std::string file_str_vel;              // string with path to mesh file.
    std::string file_str_pres;

    file_str_vel = globals->path + SEP + "InitialConditions/vel_init.txt";
    file_str_pres = globals->path + SEP + "InitialConditions/pres_init.txt";

    // in stream for input.in file
    std::ifstream velFile(file_str_vel, std::ifstream::in);
    std::ifstream presFile(file_str_pres, std::ifstream::in);

    if (velFile.is_open())
    {
        outstring << std::endl << "Found initial conditions file: " + file_str_vel << std::endl;
        globals->Output->Write(OUT_MESSAGE, &outstring);

        // std::getline(velFile, line);                               // READ HEADER
        i = 0;
        while (std::getline(velFile, line))
        {
            if (i >= globals->face_num) break;

            std::istringstream line_ss(line);

            std::getline(line_ss >> std::ws,val,' ');
            v(i) = std::stod(val);

            std::getline(line_ss >> std::ws,val,' ');
            dvdt(i,0) = std::stod(val);

            std::getline(line_ss >> std::ws,val,' ');
            dvdt(i,1) = std::stod(val);

            std::getline(line_ss >> std::ws,val,' ');
            dvdt(i,2) = std::stod(val);

            if (globals->solver_type == RK4) {
                std::getline(line_ss >> std::ws,val,' ');
                dvdt(i,3) = std::stod(val);
            }

            i++;
        }

        velFile.close();
    }
    else
    {

        v(i) = 0.0;
        dvdt(i,0) = 0.0;
        dvdt(i,1) = 0.0;
        dvdt(i,2) = 0.0;

        outstring << "WARNING: NO INITIAL CONDITION FILE FOUND " + file_str_vel << std::endl;
        globals->Output->Write(ERR_MESSAGE, &outstring);
    }

    if (presFile.is_open())
    {
        outstring << std::endl << "Found initial conditions file: " + file_str_pres << std::endl;
        globals->Output->Write(OUT_MESSAGE, &outstring);

        // std::getline(presFile, line);                               // READ HEADER
        i = 0;
        while (std::getline(presFile, line))
        {
            if (i >= globals->node_num) break;

            std::istringstream line_ss(line);

            std::getline(line_ss >> std::ws,val,' ');
            p(i) = std::stod(val);

            std::getline(line_ss >> std::ws,val,' ');
            dpdt(i,0) = std::stod(val);

            std::getline(line_ss >> std::ws,val,' ');
            dpdt(i,1) = std::stod(val);

            std::getline(line_ss >> std::ws,val,' ');
            dpdt(i,2) = std::stod(val);

            if (globals->solver_type == RK4) {
                std::getline(line_ss >> std::ws,val,' ');
                dpdt(i,3) = std::stod(val);
            }

            i++;
        }

        presFile.close();
    }
    else
    {

        v(i) = 0.0;
        dvdt(i,0) = 0.0;
        dvdt(i,1) = 0.0;
        dvdt(i,2) = 0.0;

        p(i) = 0.0;
        dpdt(i,0) = 0.0;
        dpdt(i,1) = 0.0;
        dpdt(i,2) = 0.0;

        if (globals->solver_type == RK4) {
            dvdt(i, 3) = 0.0;
            dpdt(i, 3) = 0.0;
        }


        outstring << "WARNING: NO INITIAL CONDITION FILE FOUND " + file_str_pres << std::endl;
        globals->Output->Write(ERR_MESSAGE, &outstring);
    }

    return 1;

};

int analyticalInitialConditions(Globals * globals, Mesh * mesh,
                          Array1D<double> & v, Array2D<double> & dvdt,
                          Array1D<double> & p, Array2D<double> & dpdt)
{
    int i, j;

    double t = 0.0*2*pi/globals->angVel.Value();
    double lat, lon, nx, ny;
    double dt = globals->timeStep.Value();

    std::array<double, 6> initArr;

    for (i=0; i<FACE_NUM; i++) {
        lat = mesh->face_centre_pos_sph(i, 0);
        lon = mesh->face_centre_pos_sph(i, 1);

        nx = mesh->face_normal_vec_map(i,0);
        ny = mesh->face_normal_vec_map(i,1);

        initArr = analyticalLTE(globals, mesh, globals->tide_type, lat, lon, t);

        v(i) = initArr[0]*nx + initArr[1]*ny;
        dvdt(i,0) = initArr[2]*nx + initArr[3]*ny;

        initArr = analyticalLTE(globals, mesh, globals->tide_type, lat, lon, t-dt);
        dvdt(i,1) = initArr[2]*nx + initArr[3]*ny;
        
        initArr = analyticalLTE(globals, mesh, globals->tide_type, lat, lon, t-2*dt);
        dvdt(i,2) = initArr[2]*nx + initArr[3]*ny;

        if (globals->solver_type == RK4) {
            initArr = analyticalLTE(globals, mesh, globals->tide_type, lat, lon, t-3*dt); 
            dvdt(i,3) = initArr[2]*nx + initArr[3]*ny;
        }

    }

    for (i=0; i<NODE_NUM; i++) {
        lat = mesh->node_pos_sph(i, 0);
        lon = mesh->node_pos_sph(i, 1);

        initArr = analyticalLTE(globals, mesh, globals->tide_type, lat, lon, t);
        p(i) = initArr[4];
        dpdt(i,0) = initArr[5];

        initArr = analyticalLTE(globals, mesh, globals->tide_type, lat, lon, t-dt);
        dpdt(i,1) = initArr[5];

        initArr = analyticalLTE(globals, mesh, globals->tide_type, lat, lon, t-2*dt);
        dpdt(i,2) = initArr[5];

        if (globals->solver_type == RK4) {
            initArr = analyticalLTE(globals, mesh, globals->tide_type, lat, lon, t-3*dt);
            dpdt(i,3) = initArr[5];
        }

    }

    return 1;

};


int writeInitialConditions(Globals * globals, Mesh * mesh,
                          Array1D<double> & v, Array2D<double> & dvdt,
                          Array1D<double> & p, Array2D<double> & dpdt)
{
    mkdir(&(globals->path +  SEP + "InitialConditions" + SEP)[0], S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP);

    FILE * outFile;
    std::string initial_condition_file;
    
    initial_condition_file = globals->path + SEP + "InitialConditions" + SEP + "vel_init.txt";
    remove(&initial_condition_file[0]);

    outFile = fopen(&initial_condition_file[0], "w");

    switch (globals->solver_type)
    {
    case AB3:
        for (int i=0; i<FACE_NUM; i++)
        {                   //   u    dudt0  dudt1  dudt2    p    dpdt0  dpdt1  dpdt2
            fprintf (outFile, "%1.6E, %1.6E, %1.6E, %1.6E\n",
                    v(i),   dvdt(i, 0),    dvdt(i, 1),    dvdt(i, 2));
        }
        break;
    case RK4:
        for (int i=0; i<FACE_NUM; i++)
        {                   //   u    dudt0  dudt1  dudt2    p    dpdt0  dpdt1  dpdt2
            fprintf (outFile, "%1.6E, %1.6E, %1.6E, %1.6E. %1.6E\n",
                    v(i),   dvdt(i, 0),    dvdt(i, 1),    dvdt(i, 2), dvdt(i, 3));
        }
        break;
    
    default:
        break;
    }
    

    fclose(outFile);

    initial_condition_file = globals->path + SEP + "InitialConditions" + SEP + "pres_init.txt";
    remove(&initial_condition_file[0]);

    outFile = fopen(&initial_condition_file[0], "w");

    switch (globals->solver_type)
    {
    case AB3:
        for (int i=0; i<NODE_NUM; i++)
        {                   //   u    dudt0  dudt1  dudt2    p    dpdt0  dpdt1  dpdt2
            fprintf (outFile, "%1.6E, %1.6E, %1.6E, %1.6E\n",
                    p(i),   dpdt(i, 0),    dpdt(i, 1),    dpdt(i, 2));
        }
        break;
    case RK4:
        for (int i=0; i<NODE_NUM; i++)
        {                   //   u    dudt0  dudt1  dudt2    p    dpdt0  dpdt1  dpdt2
            fprintf (outFile, "%1.6E, %1.6E, %1.6E, %1.6E. %1.6E\n",
                    p(i),   dpdt(i, 0),    dpdt(i, 1),    dpdt(i, 2), dpdt(i, 3));
        }
        break;
    
    default:
        break;
    }

    fclose(outFile);
    return 1;

};


int getInitialConditions(Globals * globals, Mesh * mesh,
                          Array1D<double> & v, Array2D<double> & dvdt,
                          Array1D<double> & p, Array2D<double> & dpdt, double time=0.0)
{
    
#ifndef TESTS
    // Set all arrays to zero
    if (globals->initial_condition == INIT_NONE) {

        for (int i=0; i<NODE_NUM; i++)
        {
            p(i) = 0.0;
            dpdt(i,0) = 0.0;
            dpdt(i,1) = 0.0;
            dpdt(i,2) = 0.0;
        }

        for (int i=0; i<FACE_NUM; i++)
        {
            v(i) = 0.0;
            dvdt(i,0) = 0.0;
            dvdt(i,1) = 0.0;
            dvdt(i,2) = 0.0;
        }
    }

    // Load arrays from file
    else if (globals->initial_condition == INIT_LOAD) {
        loadInitialConditions(globals, mesh, v, dvdt, p, dpdt);
    }

    // Use analytical formula for initial conditions
    else if (globals->initial_condition == INIT_ANALYTICAL) {
        analyticalInitialConditions(globals, mesh, v, dvdt, p, dpdt);
    }

#elif defined(TEST_SW1)

    // Get SW test initial conditions
    double a, ut, vt, u0, lat, lon, nx, ny;
    double r = globals->radius.Value();
    double h = globals->h.Value();

        
    // Cosine bell will advect around once every 12 orbital periods
    u0 = r*globals->angVel.Value()/12.0;
    
    a = 0.05;

    // Assign velocity field
    for (int i=0; i<FACE_NUM; i++) 
    {
        v(i) = 0.0;

        lat = mesh->face_centre_pos_sph(i, 0);
        lon = mesh->face_centre_pos_sph(i, 1);

        nx = mesh->face_normal_vec_map(i, 0);
        ny = mesh->face_normal_vec_map(i, 1);

        ut = cos(lat)*cos(a) + sin(lat)*cos(lon)*sin(a);
        ut *= u0;

        vt = -u0*sin(lon)*sin(a);

        v(i) = ut*nx + vt*ny;
    }

    double dist = 0.0;
    double lat_bell = 0*pi/180;
    double lon_bell = 180.0*pi/180;
    // Assign height field
    for (int i=0; i<NODE_NUM; i++) 
    {
        lat = mesh->node_pos_sph(i,0);
        lon = mesh->node_pos_sph(i,1);

        distanceBetweenSph(dist, lat_bell, lat, lon_bell, lon, r);

        if (dist< r/3.0) {
            p(i) = h*0.5 * ( 1 + cos(pi*dist/(r/3.0) ) );
        }
        else {
            p(i) = 0.0;
        }
    }

#elif defined(TEST_SW2)
        // Get SW test initial conditions
    double a, ut, vt, u0, lat, lon, nx, ny;
    double r = globals->radius.Value();
    double h = globals->h.Value();
    double g = globals->g.Value();
    double Omega = globals->angVel.Value();

        
    // Cosine bell will advect around once every 12 orbital periods
    u0 = r*globals->angVel.Value()/12.0;
    
    a = pi/4.0;

    // Assign velocity field
    for (int i=0; i<FACE_NUM; i++) 
    {
        v(i) = 0.0;

        lat = mesh->face_centre_pos_sph(i, 0);
        lon = mesh->face_centre_pos_sph(i, 1);

        nx = mesh->face_normal_vec_map(i, 0);
        ny = mesh->face_normal_vec_map(i, 1);

        ut = cos(lat)*cos(a) + sin(lat)*cos(lon)*sin(a);
        ut *= u0;

        vt = -u0*sin(lon)*sin(a);

        v(i) = ut*nx + vt*ny;
    }


    // Assign height field
    for (int i=0; i<NODE_NUM; i++) 
    {
        lat = mesh->node_pos_sph(i,0);
        lon = mesh->node_pos_sph(i,1);

        p(i) = g*h - (r*Omega*u0 + u0*u0*0.5)*pow(-cos(lon)*cos(lat)*sin(a) + sin(lat)*cos(a), 2.0);
        p(i) = p(i)/g - h;
    }

#elif defined(TEST_SW6)
        // Get SW test initial conditions
    double a, ut, vt, u0, lat, lon, nx, ny;
    double r = globals->radius.Value();
    double h = globals->h.Value();
    double g = globals->g.Value();
    double Omega = globals->angVel.Value();
    double omega = 0.10762479 * Omega;
    double R = 4.0;
    double K = omega;

        
    // Cosine bell will advect around once every 12 orbital periods
    u0 = r*globals->angVel.Value()/12.0;
    
    a = 0.05;

    // Assign velocity field
    double clat, slat;
    for (int i=0; i<FACE_NUM; i++) 
    {
        v(i) = 0.0;

        lat = mesh->face_centre_pos_sph(i, 0);
        lon = mesh->face_centre_pos_sph(i, 1);

        nx = mesh->face_normal_vec_map(i, 0);
        ny = mesh->face_normal_vec_map(i, 1);

        clat = cos(lat);
        slat = sin(lat);

        ut = r*(omega*cos(lat) + K*pow(clat, R-1) * (R*slat*slat - clat*clat)*cos(R*lon));

        vt = -r*K*R*pow(clat,R-1) * slat * sin(R*lon);

        v(i) = ut*nx + vt*ny;
    }


    // Assign height field
    double A, B, C;
    double rr = r*r;
    for (int i=0; i<NODE_NUM; i++) 
    {
        lat = mesh->node_pos_sph(i,0);
        lon = mesh->node_pos_sph(i,1);

        clat = cos(lat);


        A = 0.5*omega*(2*Omega + omega)*clat*clat 
            + 0.25*K*K * pow(clat, 2*R) 
                * ( (R+1) *clat*clat + (2*R*R - R - 2)
                    -2*R*R*1.0/(clat*clat));

        B = 2*(Omega + omega)*K / ((R+1)*(R+2)) * pow(clat, R) 
            * ( (R*R + 2*R + 2) - (R+1)*(R+1)*clat*clat);

        C = 0.25 * K*K * pow(clat, 2*R) * ( (R+1)*clat*clat - (R+2)); 


        p(i) = g*h + rr*A + rr*B*cos(R*lon) + rr*C*cos(2*R*lon);
        p(i) = p(i)/g - h;
    }



#elif defined(TEST_GAUSS_HILLS)
        // Get SW test initial conditions
    double k, ut, vt, u0, lat, lon, nx, ny;
    double r = globals->radius.Value();
    double h = globals->h.Value();
    double g = globals->g.Value();
    double Omega = globals->angVel.Value();
    double T =  12 * globals->period.Value();
        
    // Cosine bell will advect around once every 12 orbital periods
    u0 = r/T;
    
    k = 6.5;

    // Assign velocity field
    for (int i=0; i<FACE_NUM; i++) 
    {
        lat = mesh->face_centre_pos_sph(i, 0);
        lon = mesh->face_centre_pos_sph(i, 1) - 2*pi*time/T;

        nx = mesh->face_normal_vec_map(i, 0);
        ny = mesh->face_normal_vec_map(i, 1);

        ut = u0*k*sin(lon)*sin(lon)*sin(2*lat)*cos(pi*time/T) + r*2*pi*cos(lat)/T;

        vt = u0*k*sin(2*lon)*cos(lat)*cos(pi*time/T);

        v(i) = ut*nx + vt*ny;
    }

    if ( !(time > 0.01) ) 
    {
        // Assign height field
        double lat0_1 = 0.0;
        double lon0_1 = 5*pi/6.0;
        double lat0_2 = 0.0;
        double lon0_2 = 7*pi/6.0;
        for (int i=0; i<NODE_NUM; i++) 
        {
            lat = mesh->node_pos_sph(i,0);
            lon = mesh->node_pos_sph(i,1);

            p(i) = gaussian(lat, lon, lat0_1, lon0_1, r, 0.001*h) 
                + gaussian(lat, lon, lat0_2, lon0_2, r, 0.001*h);
            
        }
    }
    

#endif

    return 1;
}