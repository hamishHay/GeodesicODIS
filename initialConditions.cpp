#include "mesh.h"
#include "globals.h"
#include "array3d.h"
#include "array2d.h"
#include "array1d.h"

int loadInitialConditions(Globals * globals, Mesh * mesh,
                          Array2D<double> & v, Array2D<double> & dvdt1, Array2D<double> & dvdt2, Array2D<double> & dvdt3,
                          Array1D<double> & p, Array1D<double> & dpdt1, Array1D<double> & dpdt2, Array1D<double> & dpdt3)
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
            // std::cout<<line<<std::endl;
            std::istringstream line_ss(line);

            std::getline(line_ss >> std::ws,val,' ');
            v(i,0) = std::stod(val);

            std::getline(line_ss >> std::ws,val,' ');
            dvdt1(i,0) = std::stod(val);

            std::getline(line_ss >> std::ws,val,' ');
            dvdt2(i,0) = std::stod(val);

            std::getline(line_ss >> std::ws,val,' ');
            dvdt3(i,0) = std::stod(val);




            std::getline(line_ss >> std::ws,val,' ');
            v(i,1) = std::stod(val);

            std::getline(line_ss >> std::ws,val,' ');
            dvdt1(i,1) = std::stod(val);

            std::getline(line_ss >> std::ws,val,' ');
            dvdt2(i,1) = std::stod(val);

            std::getline(line_ss >> std::ws,val,' ');
            dvdt3(i,1) = std::stod(val);




            std::getline(line_ss >> std::ws,val,' ');
            p(i) = std::stod(val);

            std::getline(line_ss >> std::ws,val,' ');
            dpdt1(i) = std::stod(val);

            std::getline(line_ss >> std::ws,val,' ');
            dpdt2(i) = std::stod(val);

            std::getline(line_ss >> std::ws,val,' ');
            dpdt3(i) = std::stod(val);

            i++;
        }

        gridFile.close();
    }
    else{

      v(i,0) = 0.0;
      dvdt1(i,0) = 0.0;
      dvdt2(i,0) = 0.0;
      dvdt3(i,0) = 0.0;

      v(i,1) = 0.0;
      dvdt1(i,1) = 0.0;
      dvdt2(i,1) = 0.0;
      dvdt3(i,1) = 0.0;

      p(i) = 0.0;
      dpdt1(i) = 0.0;
      dpdt2(i) = 0.0;
      dpdt3(i) = 0.0;


        outstring << "WARNING: NO INITIAL CONDITION FILE FOUND " + file_str << std::endl;
        globals->Output->Write(ERR_MESSAGE, &outstring);
    }

    return 1;

};


// int saveInitialConditions(Globals * globals, Mesh * mesh,
//                           Array2D<double> & v, Array2D<double> & dvdt1, Array2D<double> & dvdt2, Array2D<double> & dvdt3,
//                           Array1D<double> & p, Array1D<double> & dpdt1, Array1D<double> & dpdt2, Array1D<double> & dpdt3)
// {
//     int node_num;
//     int i, j;
//     std::ostringstream outstring;
//
//     node_num = mesh->node_num;
//     std::string line, val;                 // strings for column and individual number
//     std::string file_str;                  // string with path to mesh file.
//
//     // file_str = globals->path + SEP + "initial_conditions.txt";
//     file_str = globals->path + SEP + "init.txt";
//
//     // in stream for input.in file
//     std::ifstream gridFile(file_str, std::ifstream::in);
//
//     if (gridFile.is_open())
//     {
//         outstring << std::endl << "Found initial conditions file: " + file_str << std::endl;
//         globals->Output->Write(OUT_MESSAGE, &outstring);
//
//         // std::getline(gridFile, line);                               // READ HEADER
//         i = 0;
//         while (std::getline(gridFile, line))
//         {
//             // std::cout<<line<<std::endl;
//             std::istringstream line_ss(line);
//
//             std::getline(line_ss >> std::ws,val,' ');
//             v(i,0) = std::stod(val);
//
//             std::getline(line_ss >> std::ws,val,' ');
//             dvdt1(i,0) = std::stod(val);
//
//             std::getline(line_ss >> std::ws,val,' ');
//             dvdt2(i,0) = std::stod(val);
//
//             std::getline(line_ss >> std::ws,val,' ');
//             dvdt3(i,0) = std::stod(val);
//
//
//
//
//             std::getline(line_ss >> std::ws,val,' ');
//             v(i,1) = std::stod(val);
//
//             std::getline(line_ss >> std::ws,val,' ');
//             dvdt1(i,1) = std::stod(val);
//
//             std::getline(line_ss >> std::ws,val,' ');
//             dvdt2(i,1) = std::stod(val);
//
//             std::getline(line_ss >> std::ws,val,' ');
//             dvdt3(i,1) = std::stod(val);
//
//
//
//
//             std::getline(line_ss >> std::ws,val,' ');
//             p(i) = std::stod(val);
//
//             std::getline(line_ss >> std::ws,val,' ');
//             dpdt1(i) = std::stod(val);
//
//             std::getline(line_ss >> std::ws,val,' ');
//             dpdt2(i) = std::stod(val);
//
//             std::getline(line_ss >> std::ws,val,' ');
//             dpdt3(i) = std::stod(val);
//
//             i++;
//         }
//
//         gridFile.close();
//     }
//     else
//     {
//         outstring << "ERROR: GRID FILE NOT FOUND AT " + file_str << std::endl;
//         globals->Output->Write(ERR_MESSAGE, &outstring);
//         globals->Output->TerminateODIS();
//     }
//
//     return 1;
//
// };
