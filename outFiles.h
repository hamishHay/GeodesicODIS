#ifndef OUTFILES_H
#define OUTFILES_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "H5Cpp.h"

#include "globals.h"
#include "H5Cpp.h"

enum mess_type{ OUT_MESSAGE, ERR_MESSAGE, GRID_MESSAGE };

class Globals;        // Forward declare globals class

class OutFiles {
private:
	std::string outName;
	std::string errName;
	std::ofstream output;
	std::ofstream error;

    int node_num;

	void WriteMessage(std::ostringstream * sstream);
	void WriteError(std::ostringstream * sstream);

	void ClearSStream(std::ostringstream * sstream);

public:
	std::string path;
	std::string dataPath;

  std::vector<std::string> * tags;

	OutFiles();

	void Write(mess_type message, std::ostringstream * sstream);
	void WelcomeMessage(void);

	void TerminateODIS(void);

  void CreateHDF5Framework(Globals *);

  void DumpData(Globals *, int, double **);

  //------------------------Objects for HDF5 Storage----------------------------

  float * u_1D;
  float * v_1D;
  float * eta_1D;
  float * press_1D;
  float * diss_1D;
  float * kinetic_1D;
  float * diss_avg_1D;
  float * kinetic_avg_1D;
  float * dummy1_1D;
  float * dummy2_2D;
  // float * harm_coeff_1D;

  hsize_t rank_field;
  hsize_t rank_harm;
  hsize_t rank_1D;

  hid_t file;

  hid_t data_space_eta;
  hid_t data_space_u;
  hid_t data_space_v;
  hid_t data_space_press;
  hid_t data_space_diss;
  hid_t data_space_kinetic;
  hid_t data_space_1D_avg;
  hid_t data_space_1D_kinetic_avg;
  hid_t data_space_dummy1;
  hid_t data_space_dummy2;
  // hid_t data_space_harm_coeff;

  hid_t mem_space_eta;
  hid_t mem_space_u;
  hid_t mem_space_v;
  hid_t mem_space_press;
  hid_t mem_space_diss;
  hid_t mem_space_kinetic;
  hid_t mem_space_1D_avg;
  hid_t mem_space_1D_kinetic_avg;
  hid_t mem_space_dummy1;
  hid_t mem_space_dummy2;
  // hid_t mem_space_harm_coeff;

  hid_t data_set_eta;
  hid_t data_set_u;
  hid_t data_set_v;
  hid_t data_set_press;
  hid_t data_set_diss;
  hid_t data_set_kinetic;
  hid_t data_set_1D_avg;
  hid_t data_set_1D_kinetic_avg;
  hid_t data_set_dummy1;
  hid_t data_set_dummy2;
  // hid_t data_set_harm_coeff;

  char * dataFilePath;

  hsize_t * dims_field;
  hsize_t * dims_1D_diss_avg;
  hsize_t * dims_harm_coeff;

  hsize_t * max_dims_field;
  hsize_t * max_dims_1D_diss_avg;
  hsize_t * max_dims_harm_coeff;

  hsize_t * start;
  hsize_t * count;

  hsize_t * start_1D;
  hsize_t * count_1D;

  hsize_t * start_harm;
  hsize_t * count_harm;

  hsize_t time_slices;

};

#endif
