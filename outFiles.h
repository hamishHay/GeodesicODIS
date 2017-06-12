#ifndef OUTFILES_H
#define OUTFILES_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

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

	void WriteMessage(std::ostringstream * sstream);
	void WriteError(std::ostringstream * sstream);

	void ClearSStream(std::ostringstream * sstream);

public:
	std::string path;
	std::string dataPath;

	OutFiles();

	void Write(mess_type message, std::ostringstream * sstream);
	void WelcomeMessage(void);

	void TerminateODIS(void);

  void CreateHDF5Framework(Globals *);

  //------------------------Objects for HDF5 Storage----------------------------

  float * eta_1D;
  float * u_1D;
  float * v_1D;
  float * diss_1D;
  float * diss_avg_1D;
  float * harm_coeff_1D;

  hsize_t rank_field;
  hsize_t rank_harm;
  hsize_t rank_1D;

  hid_t file;

  hid_t data_space_eta;
  hid_t data_space_u;
  hid_t data_space_v;
  hid_t data_space_diss;
  hid_t data_space_1D_avg;
  hid_t data_space_harm_coeff;

  hid_t mem_space_eta;
  hid_t mem_space_u;
  hid_t mem_space_v;
  hid_t mem_space_diss;
  hid_t mem_space_1D_avg;
  hid_t mem_space_harm_coeff;

  hid_t data_set_eta;
  hid_t data_set_u;
  hid_t data_set_v;
  hid_t data_set_diss;
  hid_t data_set_1D_avg;
  hid_t data_set_harm_coeff;

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
