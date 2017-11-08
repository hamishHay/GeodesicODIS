#include "mesh.h"
#include "globals.h"
#include "array3d.h"
#include "array2d.h"
#include "array1d.h"

int loadInitialConditions(Globals * globals, Mesh * mesh, Array2D<double> & vel, Array1D<double> & pres)
{
    int node_num;
    int i, j;
    std::ostringstream outstring;

    node_num = mesh->node_num;
    std::string line, val;                 // strings for column and individual number
    std::string file_str;                  // string with path to mesh file.

    file_str = globals->path + SEP + "init.txt";

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
            // std::cout<<line<<std::endl;
            std::istringstream line_ss(line);
            std::getline(line_ss >> std::ws,val,' ');                 // COL 0: Node ID number
            vel(i,0) = std::stof(val);
            //
            std::getline(line_ss >> std::ws,val,' ');                 // COL 1: Node Latitude
            vel(i,1) = std::stof(val);

            std::getline(line_ss >> std::ws,val,' ');                 // COL 2: Node Longitude
            pres(i) = std::stof(val);

            i++;
        }

        gridFile.close();
    }
    else
    {
        outstring << "ERROR: GRID FILE NOT FOUND AT " + file_str << std::endl;
        globals->Output->Write(ERR_MESSAGE, &outstring);
        globals->Output->TerminateODIS();
    }

    return 1;

};
