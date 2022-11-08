#include "mesh.h"
#include "globals.h"
#include "array3d.h"
#include "array2d.h"
#include "array1d.h"
#include "analyticalLTE.h"

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
    int face_num = globals->face_num;
    int node_num = globals->node_num;

    double t = 0.0*2*pi/globals->angVel.Value();
    double lat, lon, nx, ny;
    double dt = globals->timeStep.Value();

    std::array<double, 6> initArr;

    for (i=0; i<face_num; i++) {
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

    for (i=0; i<node_num; i++) {
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

    int face_num = globals->face_num;
    int node_num = globals->node_num;

    FILE * outFile;
    std::string initial_condition_file;
    
    initial_condition_file = globals->path + SEP + "InitialConditions" + SEP + "vel_init.txt";
    remove(&initial_condition_file[0]);

    outFile = fopen(&initial_condition_file[0], "w");

    switch (globals->solver_type)
    {
    case AB3:
        for (int i=0; i<face_num; i++)
        {                   //   u    dudt0  dudt1  dudt2    p    dpdt0  dpdt1  dpdt2
            fprintf (outFile, "%1.6E, %1.6E, %1.6E, %1.6E\n",
                    v(i),   dvdt(i, 0),    dvdt(i, 1),    dvdt(i, 2));
        }
        break;
    case RK4:
        for (int i=0; i<face_num; i++)
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
        for (int i=0; i<node_num; i++)
        {                   //   u    dudt0  dudt1  dudt2    p    dpdt0  dpdt1  dpdt2
            fprintf (outFile, "%1.6E, %1.6E, %1.6E, %1.6E\n",
                    p(i),   dpdt(i, 0),    dpdt(i, 1),    dpdt(i, 2));
        }
        break;
    case RK4:
        for (int i=0; i<node_num; i++)
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


