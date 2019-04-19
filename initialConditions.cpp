#include "mesh.h"
#include "globals.h"
#include "array3d.h"
#include "array2d.h"
#include "array1d.h"

#include <iostream>
#include <sstream>
#include <iomanip>

int loadInitialConditions(Globals * globals, Mesh * mesh,
                          Array2D<double> & v, Array3D<double> & dvdt,
                          Array1D<double> & p, Array2D<double> & dpdt)
{
    int i, j;
    std::ostringstream outstring;

    std::string line, val;                 // strings for column and individual number
    std::string file_str;                  // string with path to mesh file.

    // file_str = globals->path + SEP + "initial_conditions.txt";
    file_str = globals->path + SEP + "init_new.txt";

    // in stream for input.in file
    std::ifstream gridFile(file_str, std::ifstream::in);

    if (gridFile.is_open())
    {
        outstring << std::endl << "Found initial conditions file: " + file_str << std::endl;
        globals->Output->Write(OUT_MESSAGE, &outstring);

        // std::getline(gridFile, line);                               // READ HEADER
        i = 0;
        while (std::getline(gridFile, line))
        {
            if (i >= globals->node_num) break;

            std::istringstream line_ss(line);

            std::getline(line_ss >> std::ws,val,' ');
            v(i,0) = std::stod(val);

            std::getline(line_ss >> std::ws,val,' ');
            dvdt(i,0,0) = std::stod(val);

            std::getline(line_ss >> std::ws,val,' ');
            dvdt(i,0,1) = std::stod(val);

            std::getline(line_ss >> std::ws,val,' ');
            dvdt(i,0,2) = std::stod(val);




            std::getline(line_ss >> std::ws,val,' ');
            v(i,1) = std::stod(val);

            std::getline(line_ss >> std::ws,val,' ');
            dvdt(i,1,0) = std::stod(val);

            std::getline(line_ss >> std::ws,val,' ');
            dvdt(i,1,1) = std::stod(val);

            std::getline(line_ss >> std::ws,val,' ');
            dvdt(i,1,2) = std::stod(val);




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

        gridFile.close();
    }
    else
    {

        v(i,0) = 0.0;
        dvdt(i,0,0) = 0.0;
        dvdt(i,0,1) = 0.0;
        dvdt(i,0,2) = 0.0;

        v(i,1) = 0.0;
        dvdt(i,1,0) = 0.0;
        dvdt(i,1,1) = 0.0;
        dvdt(i,1,2) = 0.0;

        p(i) = 0.0;
        dpdt(i,0) = 0.0;
        dpdt(i,1) = 0.0;
        dpdt(i,2) = 0.0;


        outstring << "WARNING: NO INITIAL CONDITION FILE FOUND " + file_str << std::endl;
        globals->Output->Write(ERR_MESSAGE, &outstring);
    }

    return 1;

};


int writeInitialConditions(Globals * globals, Mesh * mesh,
                          Array2D<double> & v, Array3D<double> & dvdt,
                          Array1D<double> & p, Array2D<double> & dpdt)
{
    int node_num = globals->node_num;

    FILE * outFile;
    std::string initial_condition_file;
    initial_condition_file = globals->path + SEP + "init_new.txt";
    remove(&initial_condition_file[0]);

    outFile = fopen(&initial_condition_file[0], "w");

    for (int i=0; i<node_num; i++)
    {                   //   u    dudt0  dudt1  dudt2    v    dvdt0  dvdt1  dvdt2    p    dpdt0  dpdt1  dpdt2
        fprintf (outFile, "%1.6E, %1.6E, %1.6E, %1.6E, %1.6E, %1.6E, %1.6E, %1.6E, %1.6E, %1.6E, %1.6E, %1.6E\n",
                 v(i,0), dvdt(i, 0, 0), dvdt(i, 0, 1), dvdt(i, 0, 2),
                 v(i,1), dvdt(i, 1, 0), dvdt(i, 1, 1), dvdt(i, 1, 2),
                 p(i),   dpdt(i, 0),    dpdt(i, 1),    dpdt(i, 2) );
    }

    fclose(outFile);
    return 1;

};
