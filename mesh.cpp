#include "mesh.h"
#include "globals.h"
#include "math.h"
#include "array2d.h"

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cstring>

//Mesh::Mesh():Mesh(2., 2.) {}; //default constructor

Mesh::Mesh(Globals * Globals, int N):node_pos_sph(N,2), node_pos_map(N,2), node_friends(N,6) //non-default constructor
{
  int i;


	globals = Globals;                   // define reference to all constants

  // InitializeArrays();

  // Read in grid file
  ReadMeshFile();

  // Calculate mapping coordinates
  CalcMappingCoords();

};

int Mesh::CalcMappingCoords(void)
{

}

// Function to read in text file containing mesh information
// The three pieces of information read and stored are:
//      1. Node ID numbers
//      2. Node positions in latitude and longitude
//      3. Neighbouring node ID numbers
int Mesh::ReadMeshFile(void)
{
  std::string line, val;                 // strings for column and individual number
  std::string file_str;                  // string with path to mesh file.
  int i;
  int node_id;

  file_str = globals->path + SEP + "input_files" + SEP + "grid_l" + std::to_string(globals->geodesic_l.Value()) + ".txt";

  // in stream for input.in file
  std::ifstream gridFile(file_str, std::ifstream::in);

  if (gridFile.is_open())
  {
    outstring << "Found mesh file: " + file_str << std::endl;
    globals->Output.Write(OUT_MESSAGE, &outstring);

    while (std::getline(gridFile, line))
    {
      std::istringstream line_ss(line);
      std::getline(line_ss >> std::ws,val,' ');                 // COL 0: Node ID number

      node_id = std::stoi(val);

      std::getline(line_ss >> std::ws,val,' ');                 // COL 1: Node Latitude
      // node_pos_sph[node_id][0] = std::atof(val.c_str());
      node_pos_sph(node_id,0) = std::atof(val.c_str());

      std::getline(line_ss >> std::ws,val,' ');                 // COL 2: Node Longitude
      // node_pos_sph[node_id][1] = std::atof(val.c_str());
      node_pos_sph(node_id,1) = std::atof(val.c_str());

      std::getline(line_ss >> std::ws,val,'{');                 // Read up to friends open bracket
      for (i = 0; i<5; i++)
      {
        std::getline(line_ss >> std::ws,val,',');               // Read each friend ID
        node_friends(node_id,i) = std::stoi(val);
      }
      std::getline(line_ss >> std::ws,val,'}');                 // Read last friend ID
      node_friends(node_id,5) = std::stoi(val);

    }

    // test(node_id, 0) = node_pos_sph[node_id][0];
    // test(node_id, 1) = node_pos_sph[node_id][1];

    gridFile.close();
  }
  else
  {
    outstring << "ERROR: GRID FILE NOT FOUND AT " + file_str << std::endl;
    globals->Output.Write(ERR_MESSAGE, &outstring);
    globals->Output.TerminateODIS();
  }
}

// int Mesh::CalculateCellArea(void)
// {
//         cellArea = new double*[ReturnLatLen()];
//         for (int i=0; i<ReturnLatLen()-1; i++) {
//            cellArea[i] = new double[ReturnLonLen()];
//            for (int j=0; j<ReturnLonLen(); j++) {
//               cellArea[i][j] = globals->radius.Value()*globals->radius.Value();
//               cellArea[i][j] *= sin(radConv*lat[i]) - sin(radConv*lat[i+1]);
//               cellArea[i][j] *= radConv*dLon;
//            }
//         }
// };

// int Mesh::CalculateDt(void){
// 	double waveSpeed;
// 	double dt;
// 	double minDistance;
//
// 	outstring << std::endl << "Calculating maximum time step: " << std::endl;
//
// 	waveSpeed = sqrt(globals->g.Value()*globals->h.Value());
//
// 	outstring << std::endl << "\t\t Average gravitational wave speed = sqrt(gh): " << waveSpeed << " m/s."<<std::endl;
//
// 	minDistance = globals->radius.Value() * sin(dLat*radConv) * dLon*radConv;
//
// 	dt = minDistance/waveSpeed*1.5;
//
// 	outstring << "\t\t Time step calculated: " << dt <<" s."<< std::endl;
//
// 	if (dt<globals->timeStep.Value()) {
// 		outstring << std::endl << "\t\t WARNING: Calculated timestep is less than user selected dt = " << globals->timeStep.Value() <<" s." <<std::endl;
// 		outstring << "\t\t Setting dt to calculated value of " << dt <<" s." <<std::endl;
//
// 		globals->timeStep.SetValue(dt);
//
// 		globals->Output.Write(OUT_MESSAGE,&outstring);
//
// 		outstring << std::endl << "WARNING: Using maximum calculated timestep as user defined dt is too large." << std::endl;
// 		globals->Output.Write(ERR_MESSAGE,&outstring);
// 	}
// 	else {
// 		outstring << std::endl << "\t\t Calculated timestep is greater than user selected dt = " << globals->timeStep.Value() <<" s." <<std::endl;
// 		outstring << "\t\t ODIS will keep user defined timestep. You could save time if you chose a larger timestep."<<std::endl;
//
// 		globals->Output.Write(OUT_MESSAGE,&outstring);
// 	}
//
//
//
//
//
// 	return 1;
// }

//
// int Mesh::InitializeArrays(void)
// {
//   int node_num, i;
//
//   // geodesic node num expression (Lee and Macdonald, 2008)
//   node_num = 10 * pow(pow(2, globals->geodesic_l.Value() - 1), 2) + 2;
//
//   // Initialize node position array
//
//   // node_pos_sph = new double*[node_num];
//   // node_pos_sph[0] = new double[node_num * 2];
//   // for (i=1; i < node_num; i++) {
//   //   node_pos_sph[i] = &node_pos_sph[0][i*2];
//   // }
//
//   //
//   // node_pos_sph = new double*[node_num];
//   // // node_pos_sph[0] = new double[node_num * 2];
//   // for (i=0; i < node_num; i++) {
//   //   node_pos_sph[i] = new double[2];
//   // }
//
//
//   // Initialize node position array
//   // node_pos_map = new double*[node_num];
//   // node_pos_map[0] = new double[node_num * 2];
//   // for (i=1; i < node_num; i++) {
//   //   node_pos_map[i] = &node_pos_map[0][i*2];
//   // }
//
//   // Initialize node friends array. Each row in node_friends contains the ID (index)
//   // to all surrounding friends. E.g, the ID's of the neighbouring nodes to node 4
//   // would be accessed via node_friends[4][0], node_friends[4][1], etc.
//   // node_friends = new int*[node_num];
//   // node_friends[0] = new int[node_num * 6];
//   // for (i=1; i < node_num; i++) {
//   //   node_friends[i] = &node_friends[0][i*6];
//   // }
//
//   return 1;
// };
