// FILE: globals.cpp
// DESCRIPTION: main file for the globals class. Globals is an object designed
//              to store useful model constants that, if necessary, are accessible
//              anywhere in the code. This file contains the constructors and member
//              functions defined in globals.h, which populate the list of global
//              variables and read their input values from the input.in file.
//
//   Version              Date                   Programmer
//     0.1      |        13/09/2016        |        H. Hay        |
// ---- Initial version of ODIS, described in Hay and Matsuyama (2016)

#ifdef _WIN32
#include <direct.h>
#include <Windows.h>
#define getcwd _getcwd
#define SEP "\\"

#elif _WIN64
#include <direct.h>
#define getcwd _getcwd
#define SEP "\\"

#elif __linux__
#include <unistd.h>
#define SEP "/"

#else
#error "OS not supported!"
#endif

#include "globals.h"
#include "outFiles.h"
#include "array1d.h"
#include "boundaryConditions.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cstring>
#include <math.h>

// Points constructor to second, optional constructor. --> Not actually necessary?
Globals::Globals() :Globals(1) {};

Globals::Globals(int action) {
  // Constructor assings stringIDs and automatically reads from the input file.

  Output = new OutFiles;

  // set simulation path to that from the Output class.
  path = Output->path;

  // copy the string path to a character array
  strcpy(cpath, path.c_str());

  // Yes, currently action is redundant.
  if (action) SetDefault();
  else {
    SetDefault();

    // Append variable IDs to all global variables, then add them to variable list.
    // Variables require a string ID so the correct variable can be matched to that
    // being read from input.in.

    // All variables described in globals.h

    // Please keep *all* string ID's in lower case.
    radius.SetStringID("radius");
    allGlobals.push_back(&radius);

    angVel.SetStringID("angular velocity");
    allGlobals.push_back(&angVel);

    loveK2.SetStringID("k2");
    allGlobals.push_back(&loveK2);

    loveH2.SetStringID("h2");
    allGlobals.push_back(&loveH2);

    loveReduct.SetStringID("love reduction factor");
    allGlobals.push_back(&loveReduct);

    h.SetStringID("ocean thickness");
    allGlobals.push_back(&h);

    H1.SetStringID("layer 1 thickness");
    allGlobals.push_back(&H1);

    H2.SetStringID("layer 2 thickness");
    allGlobals.push_back(&H2);

    shell_thickness.SetStringID("shell thickness");
    allGlobals.push_back(&shell_thickness);

    g.SetStringID("surface gravity");
    allGlobals.push_back(&g);

    g_reduced.SetStringID("reduced gravity");
    allGlobals.push_back(&g_reduced);

    a.SetStringID("semimajor axis");
    allGlobals.push_back(&a);

    e.SetStringID("eccentricity");
    allGlobals.push_back(&e);

    den1.SetStringID("density layer 1");
    allGlobals.push_back(&den1);

    den2.SetStringID("density layer 2");
    allGlobals.push_back(&den2);

    den_ratio.SetStringID("density ratio");
    allGlobals.push_back(&den_ratio);

    radius_bottom.SetStringID("radius bottom");
    allGlobals.push_back(&radius_bottom);

    radius_top.SetStringID("radius top");
    allGlobals.push_back(&radius_top);

    radius_ratio.SetStringID("radius ratio");
    allGlobals.push_back(&radius_ratio);

    theta.SetStringID("obliquity");
    allGlobals.push_back(&theta);

    timeStep.SetStringID("time step");
    allGlobals.push_back(&timeStep);

    converge.SetStringID("converge");
    allGlobals.push_back(&converge);

    alpha.SetStringID("friction coefficient");
    allGlobals.push_back(&alpha);

	l_max.SetStringID("sh degree");
	allGlobals.push_back(&l_max);

    dLat.SetStringID("latitude spacing");     //TODO remove
    allGlobals.push_back(&dLat);

    dLon.SetStringID("longitude spacing");
    allGlobals.push_back(&dLon);              //TODO remove

    geodesic_l.SetStringID("geodesic grid level");
    allGlobals.push_back(&geodesic_l);

    period.SetStringID("orbital period");
    allGlobals.push_back(&period);

    endTime.SetStringID("simulation end time");
    allGlobals.push_back(&endTime);

    potential.SetStringID("potential");
    allGlobals.push_back(&potential);

    friction.SetStringID("friction type");
    allGlobals.push_back(&friction);

    surface.SetStringID("surface type");
    allGlobals.push_back(&surface);

    solver.SetStringID("solver type");
    allGlobals.push_back(&solver);

    init.SetStringID("initial conditions");
    allGlobals.push_back(&init);

    diss_avg.SetStringID("dissipation avg output");
    allGlobals.push_back(&diss_avg);

    kinetic_avg.SetStringID("kinetic avg output");
    allGlobals.push_back(&kinetic_avg);

    work.SetStringID("work flux output");
    allGlobals.push_back(&work);

    field_displacement_output.SetStringID("displacement output");
    allGlobals.push_back(&field_displacement_output);

    field_velocity_output.SetStringID("velocity output");
    allGlobals.push_back(&field_velocity_output);

    field_pressure_output.SetStringID("pressure output");
    allGlobals.push_back(&field_pressure_output);

    field_diss_output.SetStringID("dissipation output");
    allGlobals.push_back(&field_diss_output);

    field_kinetic_output.SetStringID("kinetic output");
    allGlobals.push_back(&field_kinetic_output);

    field_dummy1_output.SetStringID("dummy1 output");
    allGlobals.push_back(&field_dummy1_output);

    field_dummy2_output.SetStringID("dummy2 output");
    allGlobals.push_back(&field_dummy2_output);

    sh_coeff_output.SetStringID("sh coefficient output");
    allGlobals.push_back(&sh_coeff_output);

    layer_num.SetStringID("layer number");
    allGlobals.push_back(&layer_num);

    outputTime.SetStringID("output time");
    allGlobals.push_back(&outputTime);

    ReadGlobals(); //Read globals from input.in file
  }



  // Convert end time from units of orbital period to seconds.
  endTime.SetValue(endTime.Value()*period.Value());

  // geodesic node num expression (Lee and Macdonald, 2008)
  node_num = 10 * pow(pow(2, geodesic_l.Value() - 1), 2) + 2;

  Output->CreateHDF5Framework(this);

    // identify friction type and select the corresponding enum value.
    if (friction.Value() == "LINEAR") fric_type = LINEAR;
    else if (friction.Value() == "QUADRATIC") fric_type = QUADRATIC;
    else {
        outstring << "ERROR: NO DRAG MODEL FOUND!" << std::endl;
        Output->Write(ERR_MESSAGE, &outstring);
        Output->TerminateODIS();
    }

    if (solver.Value() == "EULER") solver_type = EULER;
    else if (solver.Value() == "AB3") solver_type = AB3;
    else {
        outstring << "ERROR: NO SOLVER FOUND!" << std::endl;
        Output->Write(ERR_MESSAGE, &outstring);
        Output->TerminateODIS();
    }

    if (potential.Value() == "ECC_RAD")         tide_type = ECC_RAD;
    else if (potential.Value() == "ECC_LIB")    tide_type = ECC_LIB;
    else if (potential.Value() == "ECC")        tide_type = ECC;
    else if (potential.Value() == "ECC_WEST")   tide_type = ECC_WEST;
    else if (potential.Value() == "ECC_EAST")   tide_type = ECC_EAST;
    else if (potential.Value() == "OBLIQ")      tide_type = OBLIQ;
    else if (potential.Value() == "OBLIQ_WEST") tide_type = OBLIQ_WEST;
    else if (potential.Value() == "OBLIQ_EAST") tide_type = OBLIQ_EAST;
    else if (potential.Value() == "FULL")       tide_type = FULL;
    else if (potential.Value() == "TOTAL")      tide_type = TOTAL;
    else if (potential.Value() == "ECC_W3")     tide_type = ECC_W3;
    else if (potential.Value() == "OBLIQ_W3")   tide_type = OBLIQ_W3;
    else {
        outstring << "ERROR: NO POTENTIAL FORCING FOUND!" << std::endl;
        Output->Write(ERR_MESSAGE, &outstring);
        Output->TerminateODIS();
    }

    if (surface.Value() == "FREE")              surface_type = FREE;
    else if (surface.Value() == "FREE_LOADING") surface_type = FREE_LOADING;
    else if (surface.Value() == "LID_LOVE")     surface_type = LID_LOVE;
    else if (surface.Value() == "LID_MEMBR")    surface_type = LID_MEMBR;
    else if (surface.Value() == "LID_NUM")      surface_type = LID_NUM;
    else if (surface.Value() == "LID_INF")      surface_type = LID_INF;
    else {
        outstring << "ERROR: NO OCEAN SURFACE BOUNDARY CONDITION FOUND!" << std::endl;
        Output->Write(ERR_MESSAGE, &outstring);
        Output->TerminateODIS();
    }

    // IF TWO LAYER MODEL
    if (layer_num.Value() == 2)
    {
        outstring << "Warning: two layer model selected. Switching to rigid lid boundary." << std::endl;
        Output->Write(OUT_MESSAGE, &outstring);

        surface_type = LID_INF_2L;
        surface.SetValue("LID_INF_2L");

        // Set ocean thickness based on top and bottom ocean radii
        // h.SetValue(radius_top.Value() - radius_bottom.Value());
        // radius_bottom.SetValue(radius_top.Value() - h.Value()*0.5);


        H1.SetValue(h.Value()*0.5);
        H2.SetValue(h.Value()*0.5);
        radius_top.SetValue(radius.Value() - shell_thickness.Value() - 0.25*h.Value());
        radius_bottom.SetValue(radius.Value() - shell_thickness.Value() - 0.75*h.Value());

        // H2.SetValue(h.Value() - H1.Value());    // ensure bottom layer thickness is consistent with h and H1

        if (H1.Value() > h.Value() || H2.Value() > h.Value())
        {
            outstring << "ERROR: OCEAN LAYER THICKNESS GREATER THAN TOTAL OCEAN THICKNESS!" << std::endl;
            Output->Write(ERR_MESSAGE, &outstring);
            Output->TerminateODIS();
        }


        radius_ratio.SetValue(radius_bottom.Value()/radius_top.Value());
        radius.SetValue(radius_top.Value());    // set main radius to the ocean top

        // std::cout<<radius_ratio.Value()<<std::endl;

        den_ratio.SetValue(den2.Value()/den1.Value());

        g_reduced.SetValue( (den2.Value() - den1.Value())/den2.Value() * g.Value() );

        loveReduct.SetValue(1. - den_ratio.Value() * pow(radius_ratio.Value(), 2.0));

        if (loveReduct.Value() < 0.0)
        {
            outstring << "ERROR: TWO LAYER FORCING FACTOR IS NEGATIVE (UNPHYSICAL)" << std::endl;
            Output->Write(ERR_MESSAGE, &outstring);
            Output->TerminateODIS();
        }
    }

    if (surface_type == FREE_LOADING) loading_factor = new double[l_max.Value() + 1];
    else if (surface_type == LID_LOVE ||
             surface_type == LID_MEMBR)
    {
        shell_factor_beta = new double[l_max.Value() + 1];
    }


    applySurfaceBCs(this);

    // Print out all constants to output.txt
    OutputConsts();
};

int Globals::ReadGlobals(void)
{
  // Function to read model input parameters from input.in only.
  // Avoids the need to recompile for each different model setup.

  // in stream for input.in file
  std::ifstream inputFile(path + SEP + "input.in",std::ifstream::in);

  // varID: string to contain the variable ID from input.in.
  // varVal: string to contain the actual variable value from input.in
  std::string varID, varVal;

  // variables to store any inputs of specific types
  int valInt;
  double valDouble;
  bool valBool = true;
  std::string valStr;

  // boolean to record whether any constant variables have been added to the list
  // or not.
  bool added=false;

  if (inputFile.is_open())
  {

    outstring << "Found input file." << std::endl;
    Output->Write(OUT_MESSAGE, &outstring);


    while (std::getline(inputFile >> std::ws, varID, ';')) //>> std::ws removes any leading whitespace
    {
      // varID now contains string up to the first ";" --> this is the stringID.

      // Convert string to lower case
      std::transform(varID.begin(), varID.end(), varID.begin(), ::tolower);
      std::getline(inputFile >> std::ws, varVal, ';');

      // string stream to convert string into numerical value
      std::stringstream value(varVal);

      added = false;
      unsigned int i = 0; //unsigned to remove warnings

      do
      {
        // iterate until correct stringID found
        if (varID == allGlobals[i]->StringID())
        {
          // assign global variable to that in val based on type
          if (allGlobals[i]->IsType("double"))
          {
            value >> valDouble; // converts stringstream to double
            // casting required here to correctly access elements in allGlobals
            // vector.
            ((GlobalVar<double > *) allGlobals[i])->SetValue(valDouble);
            added = true;
            allGlobals[i]->Added(added);
          }
          else if (allGlobals[i]->IsType("int"))
          {
            value >> valInt;
            ((GlobalVar<int> *) allGlobals[i])->SetValue(valInt);
            added = true;
            allGlobals[i]->Added(added);
          }
          else if (allGlobals[i]->IsType("bool"))
          {
            value >> valStr; // this should be valBool? maybe >> doesn't work with bool
            if (valStr == "false") valBool = false;
            else if (valStr == "true") valBool = true;
            ((GlobalVar<bool> *) allGlobals[i])->SetValue(valBool);
            added = true;
            allGlobals[i]->Added(added);
          }
          else if (allGlobals[i]->IsType("string")) {
            value >> valStr;
            ((GlobalVar<std::string> *) allGlobals[i])->SetValue(valStr);
            added = true;
            allGlobals[i]->Added(added);
          }
          else {
            outstring << "Unsused type read for " << allGlobals[i]->StringID() << ". " << std::endl;
            Output->Write(ERR_MESSAGE, &outstring);
            Output->TerminateODIS();
          };
        }
        i++;
      } while (i < allGlobals.size() && !added);

      // The end of each line in input.in is the purely for verbosity. This line
      // assigns the third column in input.in to varID, but it is never used.
      std::getline(inputFile, varID, ';');
    };

    // reading of input file is complete.
    inputFile.close();

    // Set obliquity to radians
    theta.SetValue(theta.Value()*pi/180.);

    // Overwrite period to match angular velocity
    period.SetValue(2.*pi/angVel.Value());


    // Add all output tags to the string array out_tags
    if (field_velocity_output.Value())      out_tags.push_back(field_velocity_output.StringID());
    if (field_displacement_output.Value())  out_tags.push_back(field_displacement_output.StringID());
    if (field_pressure_output.Value())      out_tags.push_back(field_pressure_output.StringID());
    if (field_diss_output.Value())          out_tags.push_back(field_diss_output.StringID());
    if (field_kinetic_output.Value())       out_tags.push_back(field_kinetic_output.StringID());
    if (diss_avg.Value())                   out_tags.push_back(diss_avg.StringID());
    if (kinetic_avg.Value())                out_tags.push_back(kinetic_avg.StringID());
    if (field_dummy1_output.Value())        out_tags.push_back(field_dummy1_output.StringID());
    if (field_dummy2_output.Value())        out_tags.push_back(field_dummy2_output.StringID());
    // if (work.Value())                       out_tags.push_back(work.StringID());
    // if (sh_coeff_output.Value())            out_tags.push_back(sh_coeff_output.StringID());

    // Check for any unassigned global variables
    for (unsigned int i = 0; i < allGlobals.size(); i++){
      if (!allGlobals[i]->IsAdded()) {
        outstring << "WARNING: Unassigned global constant: " << allGlobals[i]->StringID() << ". " << std::endl;
        outstring << "Assigning default value could be risky..." << std::endl;
        outstring << "Using Titan value of ";

        if (allGlobals[i]->IsType("double")) {
          outstring << ((GlobalVar<double > *) allGlobals[i])->Value() << std::endl << std::endl;
        }
        else if (allGlobals[i]->IsType("int")) {
          outstring << ((GlobalVar<int > *) allGlobals[i])->Value() << std::endl << std::endl;
        }
        else if (allGlobals[i]->IsType("bool")) {
          outstring << ((GlobalVar<bool > *) allGlobals[i])->Value() << std::endl << std::endl;
        }
        else if (allGlobals[i]->IsType("string")) {
          outstring << ((GlobalVar<std::string > *) allGlobals[i])->Value() << std::endl << std::endl;
        }

        Output->Write(ERR_MESSAGE, &outstring);
      }
    }
  }
  else
  {
    // if input.in cannot be found, terminate the program.
    outstring << "Unable to open 'input.in' file." << std::endl << std::endl;
    Output->Write(ERR_MESSAGE, &outstring);
    Output->TerminateODIS();
  }

  return 0;
};

void Globals::SetDefault(void)
{
  // Member function acts as a backup in case no variables are assigned. Perhaps
  // it would be safer to simply require a fully populated input.in file?

  //Satellite Radius
  radius.SetValue(2574.73e3);

  //k_2 love number
  loveK2.SetValue(.120);

  //h_2 love number
  loveH2.SetValue(.2);

  //Love reduction factor
  //loveReduct.SetValue(1 + loveK2.Value() - loveH2.Value());
  loveReduct.SetValue(0.920012);
  //Ocean thickness
  h.SetValue(400);

  //Surface gravity
  g.SetValue(1.35);

  //Semimajor axis
  a.SetValue(1221.87e6);

  //Eccentricity
  e.SetValue(0.0288);

  //Obliquity
  theta.SetValue(0.32*pi / 180.);

  //Coefficient of linear friction
  alpha.SetValue(2.28e-7);

  //Maximum spherical harmonic degree for expansion
  l_max.SetValue(10);

  //Lat and long spacing
  dLat.SetValue(45);
  dLon.SetValue(90);

  //Geodesic grid resolution
  geodesic_l.SetValue(1);

  //Orbital period
  period.SetValue(16. * 24 * 60 * 60);

  //Time step
  timeStep.SetValue(40.);

  //Convergence Criteria
  converge.SetValue(1e-7);

  //Angular velocity
  angVel.SetValue(4.56e-6);

  //Orbital period
  period.SetValue(2*pi/angVel.Value());

  //Simulation end time
  endTime.SetValue(80.);

  //Potential
  potential.SetValue("ECC");

  //Friction Type
  friction.SetValue("QUADRATIC");

  //Initial conditions
  init.SetValue(false);

  //Output dissipated energy?
  diss_avg.SetValue(true);

  //Output kinetic energy?
  kinetic_avg.SetValue(false);

  //Output work flux?
  work.SetValue(false);

  layer_num.SetValue(1);

  //Output time in fraction of orbital period. Default is set to output
  //every period.
  outputTime.SetValue(1);
};

void Globals::OutputConsts(void)
{
  const int varWidth = 30;
  const char separator = ' ';
  // Print all constants to output file.
  outstring << "Model parameters: " << std::endl;
  for (unsigned int i = 0; i < allGlobals.size(); i++){

    outstring <<"\t\t "<<std::left<< std::setw(varWidth) << std::setfill(separator)<<allGlobals[i]->StringID();
    if (allGlobals[i]->IsType("double")) {
      outstring << ((GlobalVar<double > *) allGlobals[i])->Value();
    }
    else if (allGlobals[i]->IsType("int")) {
      outstring << ((GlobalVar<int > *) allGlobals[i])->Value();
    }
    else if (allGlobals[i]->IsType("bool")) {
      outstring << ((GlobalVar<bool > *) allGlobals[i])->Value();
    }
    else if (allGlobals[i]->IsType("string")) {
      outstring << ((GlobalVar<std::string > *) allGlobals[i])->Value();
    }

    Output->Write(OUT_MESSAGE, &outstring);
  }
};


// int Globals::GetShellCoeffs(void)
// {
//   double hs;
//
//   hs = hShell.Value();
//
//   // Check if ice shell value is within table range
//   if ((hs < 1e3) || (hs > 100e3)) {
//     outstring << "ERROR:\t\t User selected ice shell thickness of " << hs/1e3;
//     outstring << " km is outside of calculated shell thickness range. " << std::endl;
//     outstring << "\t\t Shell behaviour is therefore undefined. "<<std::endl;
//
//     Output.Write(ERR_MESSAGE, &outstring);
//
//     Output.TerminateODIS();
//   }
//   // if ice shell thickness is not an integer value in kilometers, we must
//   // interpolate
//   else if (fmod(hs/1e3,1.0) > 1e-8) {
//     double x_l, x_r, x_interp;       // i and i+1 and interpolated beta values
//     double dhs;                      // distance from nearest int shell thickness
//     int count_col;                   // counter for counting file column position
//     int count_row;                   // counter for counting file row position
//     std::string line, val;           // strings for column and individual number
//
//     outstring << "WARNING:\t User selected shell thickness of " << hs/1e3;
//     outstring << " km is not an integer value of kilometers." << std::endl;
//     outstring << "\t\t ODIS will interpolate shell coefficients from between ";
//     outstring << int(hs)/1000 << " km and " << int(hs)/1000 + 1 << " km." << std::endl;
//
//     Output.Write(ERR_MESSAGE, &outstring);
//
//     // in stream for input.in file
//     std::ifstream betaFile(path + SEP + "SHELL_COEFFS" + SEP + "ENCELADUS" + SEP + "beta_hs1km_to_hs100km_lmax30.txt",std::ifstream::in);
//
//     dhs = fmod(hs/1e3,1.0);
//
//     // --------- h1 --------- h ------------- h2 --------
//     //             <--------->
//     //                 dhs
//
//     count_col = 1;
//     count_row = 1;
//     if (betaFile.is_open())
//     {
//
//       outstring << "Beta coefficients found." << std::endl;
//       Output.Write(OUT_MESSAGE, &outstring);
//
//       while (std::getline(betaFile, line) && count_row <= l_max.Value())
//       {
//         x_l = 0.0;
//         x_r = 0.0;
//         x_interp = 0.0;
//
//         std::istringstream line_ss(line);
//         while (std::getline(line_ss,val,'\t'))
//         {
//           if (count_col == int(hs)/1000) {          // get h1 value
//             x_l = std::atof(val.c_str());
//           }
//           else if (count_col == int(hs)/1000+1)     // get h2 value
//           {
//             x_r = std::atof(val.c_str());
//             count_col = 1;
//             break;
//           }
//           count_col++;
//         }
//
//         x_interp = x_l + dhs * (x_r - x_l)/1.0;     // linear interpolate to find h
//         beta_shell[count_row] = x_interp;
//         count_row++;
//       };
//
//       betaFile.close();
//     }
//
//     std::ifstream nuFile(path + SEP + "SHELL_COEFFS" + SEP + "ENCELADUS" + SEP + "nu_hs1km_to_hs100km_l2.txt",std::ifstream::in);
//
//     count_col = 1;
//     count_row = 1;
//     if (nuFile.is_open())
//     {
//
//       outstring << "Nu coefficients found." << std::endl;
//       Output.Write(OUT_MESSAGE, &outstring);
//
//       while (std::getline(nuFile, line) && count_row <= 2)
//       {
//         x_l = 0.0;
//         x_r = 0.0;
//         x_interp = 0.0;
//
//         std::istringstream line_ss(line);
//         while (std::getline(line_ss,val,'\t'))
//         {
//           if (count_col == int(hs)/1000) {
//             x_l = std::atof(val.c_str());
//           }
//           else if (count_col == int(hs)/1000+1)
//           {
//             x_r = std::atof(val.c_str());
//             count_col = 1;
//             break;
//           }
//           count_col++;
//         }
//
//         x_interp = x_l + dhs * (x_r - x_l)/1.0;
//         nu_shell = x_interp;
//
//         count_row++;
//       };
//
//       nuFile.close();
//     }
//     loveReduct.SetValue(nu_shell);
//
//
//   }
//   else
//   {
//     double x;                        // beta values
//     int count_col;                   // counter for counting file column position
//     int count_row;                   // counter for counting file row position
//     std::string line, val;           // strings for column and individual number
//
//     std::ifstream betaFile(path + SEP + "SHELL_COEFFS" + SEP + "ENCELADUS" + SEP + "beta_hs1km_to_hs100km_lmax30.txt",std::ifstream::in);
//
//     count_col = 1;
//     count_row = 1;
//     if (betaFile.is_open())
//     {
//
//       outstring << "Beta coefficients found." << std::endl;
//       Output.Write(OUT_MESSAGE, &outstring);
//
//       while (std::getline(betaFile, line) && count_row <= l_max.Value())
//       {
//         x = 0.0;
//         std::istringstream line_ss(line);
//         while (std::getline(line_ss,val,'\t'))
//         {
//           if (count_col == int(hs)/1000) {
//             x = std::atof(val.c_str());
//             count_col = 1;
//             break;
//
//           }
//           count_col++;
//         }
//
//         beta_shell[count_row] = x;
//         count_row++;
//       };
//
//       betaFile.close();
//     }
//
//     std::ifstream nuFile(path + SEP + "SHELL_COEFFS" + SEP + "ENCELADUS" + SEP + "nu_hs1km_to_hs100km_l2.txt",std::ifstream::in);
//
//     count_col = 1;
//     count_row = 1;
//     if (nuFile.is_open())
//     {
//
//       outstring << "Nu coefficients found." << std::endl;
//       Output.Write(OUT_MESSAGE, &outstring);
//
//       while (std::getline(nuFile, line) && count_row <= 2)
//       {
//         x = 0.0;
//
//         std::istringstream line_ss(line);
//         while (std::getline(line_ss,val,'\t'))
//         {
//           if (count_col == int(hs)/1000) {
//             x = std::atof(val.c_str());
//
//             count_col = 1;
//             break;
//
//           }
//           count_col++;
//         }
//         nu_shell = x;
//
//         count_row++;
//       };
//
//       nuFile.close();
//     }
//
//     loveReduct.SetValue(nu_shell);
//   }
//
//   beta_shell[0] = 0.0;
//   // loveReduct.SetValue(1.0);
//
// };


Globals::~Globals()
{
    delete loading_factor;
}
