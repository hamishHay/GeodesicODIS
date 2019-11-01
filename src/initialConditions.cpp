#include "mesh.h"
#include "globals.h"
#include "array3d.h"
#include "array2d.h"
#include "array1d.h"

#include <iostream>
#include <sstream>
#include <iomanip>

int loadInitialConditions(Globals * globals, Mesh * mesh,
                          Array1D<double> & v, Array2D<double> & dvdt,
                          Array1D<double> & p, Array2D<double> & dpdt)
{
    int i, j;
    std::ostringstream outstring;

    std::string line, val;                 // strings for column and individual number
    std::string file_str_vel;                  // string with path to mesh file.
    std::string file_str_pres;

    // file_str_vel = globals->path + SEP + "initial_conditions.txt";
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


        outstring << "WARNING: NO INITIAL CONDITION FILE FOUND " + file_str_pres << std::endl;
        globals->Output->Write(ERR_MESSAGE, &outstring);
    }

    return 1;

};


int writeInitialConditions(Globals * globals, Mesh * mesh,
                          Array1D<double> & v, Array2D<double> & dvdt,
                          Array1D<double> & p, Array2D<double> & dpdt)
{
    int face_num = globals->face_num;
    int node_num = globals->node_num;

    FILE * outFile;
    std::string initial_condition_file;
    
    initial_condition_file = globals->path + SEP + "InitialConditions/vel_init.txt";
    remove(&initial_condition_file[0]);

    outFile = fopen(&initial_condition_file[0], "w");

    for (int i=0; i<face_num; i++)
    {                   //   u    dudt0  dudt1  dudt2    p    dpdt0  dpdt1  dpdt2
        // fprintf (outFile, "%1.6E, %1.6E, %1.6E, %1.6E, %1.6E, %1.6E, %1.6E, %1.6E\n",
        //          v(i),   dvdt(i, 0),    dvdt(i, 1),    dvdt(i, 2),
        //          p(i),   dpdt(i, 0),    dpdt(i, 1),    dpdt(i, 2) );
        fprintf (outFile, "%1.6E, %1.6E, %1.6E, %1.6E\n",
                 v(i),   dvdt(i, 0),    dvdt(i, 1),    dvdt(i, 2));
    }

    fclose(outFile);

    initial_condition_file = globals->path + SEP + "InitialConditions/pres_init.txt";
    remove(&initial_condition_file[0]);

    outFile = fopen(&initial_condition_file[0], "w");

    for (int i=0; i<node_num; i++)
    {                   //   u    dudt0  dudt1  dudt2    p    dpdt0  dpdt1  dpdt2
        fprintf (outFile, "%1.6E, %1.6E, %1.6E, %1.6E\n",
                 p(i),   dpdt(i, 0),    dpdt(i, 1),    dpdt(i, 2) );
    }

    fclose(outFile);
    return 1;

};
