#ifdef _WIN32
#include <direct.h>
#include <Windows.h>
#define getcwd _getcwd

#elif _WIN64
#include <direct.h>
#define getcwd _getcwd

#elif __linux__
#include <unistd.h>

#else
#error "OS not supported!"
#endif

#include "outFiles.h"
#include "globals.h"
#include "H5Cpp.h"
#include <cstring>
#include <sstream>
#include <vector>
#include <unistd.h>
#include <libgen.h>

#include <sys/stat.h>

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

OutFiles::OutFiles() {
	char buffer[PATH];

#ifdef _WIN64
	GetModuleFileName(NULL, buffer, sizeof(buffer)-1);
#elif __linux__
	readlink("/proc/self/exe", buffer, sizeof(buffer)-1);
#endif

  char buff[PATH];
  memset(buff, 0, sizeof(buff));

  /* Note we use sizeof(buf)-1 since we may need an extra char for NUL. */
  if (readlink("/proc/self/exe", buff, sizeof(buff)-1) < 0)
  {
    // TODO - add some clever error checks in here
   /* There was an error...  Perhaps the path does not exist
    * or the buffer is not big enough.  errno has the details. */
   perror("readlink");
  }

  // char *buff_copy;
  // buff_copy = strdup(buff);

  // dirname(buff_copy);
  std::cout<<buff<<std::endl;
  // getcwd( buff, PATH );

  path = buff; // set char array to std::string
  path = path.substr(0, path.find_last_of("\\/"));

	outName = path + SEP + "OUTPUT.txt";
	errName = path + SEP + "ERROR.txt";
	dataPath = path + SEP + "DATA" + SEP + "data.h5";

  std::cout<<path + SEP + "OUTPUT.txt"<<std::endl;

	remove(&outName[0]); //Converts std::string to char array
	output.open(&outName[0], std::ofstream::out | std::ofstream::app);
	if (!output.is_open()) {
		std::cerr << "Couldn't open '" << outName << "'." << std::endl;
		TerminateODIS();
	}

	output.close();

	WelcomeMessage();

	remove(&errName[0]); //Converts std::string to char array
	error.open(&errName[0], std::ofstream::out | std::ofstream::app);
	if (!error.is_open()) {
		std::cerr << "Couldn't open '" << outName << "'." << std::endl;
		TerminateODIS();
	}

	error << "ODIS error file. Warnings and model termination errors will be written here." << std::endl << std::endl;
	error.close();
};

void OutFiles::WriteMessage(std::ostringstream * sstream) {
	output.open(&outName[0], std::ofstream::out | std::ofstream::app);
	output << (*sstream).str() << std::endl;
	output.close();
};

void OutFiles::WriteError(std::ostringstream * sstream) {
	error.open(&errName[0], std::ofstream::out | std::ofstream::app);
	error << (*sstream).str() << std::endl;
	error.close();
}

void OutFiles::Write(mess_type message, std::ostringstream * sstream) {
	switch (message) {
	case OUT_MESSAGE:
		WriteMessage(sstream);
		break;

	case ERR_MESSAGE:
		WriteError(sstream);
		break;
	default:
        break;
    }

	ClearSStream(sstream);
};

void OutFiles::WelcomeMessage(void) {
	std::ostringstream sstream;

	sstream <<std::endl<<std::endl<< "        | Welcome to Ocean Dissipation in Icy Satellites (ODIS)! |" << std::endl;
	sstream << "        |--------------------------------------------------------|" << std::endl;
	sstream << "        |               ___________ _____ _____                  |" << std::endl;
	sstream << "        |              |  _  |  _  \\_   _/  ___|                 |" << std::endl;
	sstream << "        |              | | | | | | | | | \\ `--.                  |" << std::endl;
	sstream << "        |              | | | | | | | | |  `--. \\                 |" << std::endl;
	sstream << "        |              \\ \\_/ / |/ / _| |_/\\__/ /                 |" << std::endl;
	sstream << "        |               \\___/|___/  \\___/\\____/                  |" << std::endl;
	sstream << "        |                                                        |" << std::endl;
	sstream << "        |--------------------------------------------------------|" << std::endl;
	sstream << "        | ODIS is a numerical finite volume ocean dynamics       |" << std::endl;
	sstream << "        | code for simulating ocean tides on icy satellites and  |" << std::endl;
	sstream << "        | other bodies.                                          |" << std::endl;
	sstream << "        |                                                        |" << std::endl;
	sstream << "        |       For additional info, contact Hamish Hay at       |" << std::endl;
	sstream << "        |       hhay@lpl.arizona.edu.                            |" <<std::endl<<std::endl<<std::endl;

	printf((sstream.str()).c_str());

	Write(OUT_MESSAGE, &sstream);
};

void OutFiles::TerminateODIS(void) {
	std::ostringstream sstream;

	sstream << "TERMINATING ODIS." << std::endl;
  std::cout << "ODIS HAS FOUND AN ERROR. TERMINATING PROGRAM." << std::endl;
	WriteError(&sstream);
	std::exit(0);
};

void OutFiles::ClearSStream(std::ostringstream * sstream) {
	(*sstream).str(std::string());
	(*sstream).clear();
};


void OutFiles::CreateHDF5Framework(Globals * globals)
{

  unsigned int i, l_max;
  double end_time, orbit_period, output_time;

  tags = &globals->out_tags;
  node_num = globals->node_num;
  l_max = globals->l_max.Value();

  end_time = (double)globals->endTime.Value();
  orbit_period = globals->period.Value();
  output_time = globals->outputTime.Value();

  #if _WIN32
    // mkdir(&(globals->path + SEP + "Grid" + SEP)[0], NULL);

    mkdir(&(globals->path +  SEP + "DATA" + SEP)[0], NULL);

  #elif __linux__
    // mkdir(&(globals->path +  SEP + "Grid" + SEP)[0], S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP);

    mkdir(&(globals->path +  SEP + "DATA" + SEP)[0], S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP);

  #endif

    // eta_1D
    eta_1D = new float[node_num];
    v_1D = new float[node_num];
    u_1D = new float[node_num];
    press_1D = new float[node_num];
    diss_1D = new float[node_num];
    kinetic_1D = new float[node_num];
    diss_avg_1D = new float[1];
    kinetic_avg_1D = new float[node_num];
    // harm_coeff_1D = new float[2*(l_max+1)*(l_max+1)];

    time_slices = (end_time/orbit_period)/output_time + 1;

    //
    rank_field = 2;
    // hsize_t rank_size = 1;
    rank_1D = 1;
    // rank_harm = 4;

    char dataFile[1024];
    std::strcpy(dataFile, (this->dataPath).c_str());

    start = new hsize_t[2];
    count = new hsize_t[2];
    //
    start_1D = new hsize_t[1];
    count_1D = new hsize_t[1];
    //
    // start_harm = new hsize_t[4];
    // count_harm = new hsize_t[4];
    //
    // Create HDF5 file
    file = H5Fcreate(dataFile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    dims_field = new hsize_t[2];
    dims_field[0] = 1;
    dims_field[1] = node_num;

    max_dims_field = new hsize_t[2];
    max_dims_field[0] = time_slices;
    max_dims_field[1] = node_num;

    dims_1D_diss_avg = new hsize_t[1];
    dims_1D_diss_avg[0] = 1;

    max_dims_1D_diss_avg = new hsize_t[1];
    max_dims_1D_diss_avg[0] = time_slices;
    // max_dims_harm_coeff = new hsize_t[4];
    // max_dims_harm_coeff[0] = time_slices;
    // max_dims_harm_coeff[1] = 2;
    // max_dims_harm_coeff[2] = harm_rows;
    // max_dims_harm_coeff[3] = harm_cols;
    //
    if (globals->field_velocity_output.Value())
    {
        data_space_u = H5Screate_simple(rank_field, max_dims_field, NULL);
        data_space_v = H5Screate_simple(rank_field, max_dims_field, NULL);
        data_set_u = H5Dcreate(file, "east velocity", H5T_NATIVE_FLOAT, data_space_u, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        data_set_v = H5Dcreate(file, "north velocity", H5T_NATIVE_FLOAT, data_space_v, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        mem_space_u = H5Screate_simple(rank_field, dims_field, NULL);
        mem_space_v = H5Screate_simple(rank_field, dims_field, NULL);
    }
    if (globals->field_displacement_output.Value())
    {
        data_space_eta = H5Screate_simple(rank_field, max_dims_field, NULL); // 2D data space
        data_set_eta = H5Dcreate(file, "displacement", H5T_NATIVE_FLOAT, data_space_eta, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        mem_space_eta = H5Screate_simple(rank_field, dims_field, NULL);
    }
    if (globals->field_pressure_output.Value())
    {
        data_space_press = H5Screate_simple(rank_field, max_dims_field, NULL); // 2D data space
        data_set_press = H5Dcreate(file, "displacement", H5T_NATIVE_FLOAT, data_space_press, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        mem_space_press = H5Screate_simple(rank_field, dims_field, NULL);
    }
    if (globals->field_pressure_output.Value())
    {
        data_space_press = H5Screate_simple(rank_field, max_dims_field, NULL); // 2D data space
        data_set_press = H5Dcreate(file, "displacement", H5T_NATIVE_FLOAT, data_space_press, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        mem_space_press = H5Screate_simple(rank_field, dims_field, NULL);
    }
    if (globals->field_diss_output.Value())
    {
        data_space_diss = H5Screate_simple(rank_field, max_dims_field, NULL);
        data_set_diss = H5Dcreate(file, "dissipated energy", H5T_NATIVE_FLOAT, data_space_diss, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        mem_space_diss = H5Screate_simple(rank_field, dims_field, NULL);
    }
    if (globals->field_kinetic_output.Value())
    {
        data_space_kinetic = H5Screate_simple(rank_field, max_dims_field, NULL); // 2D data space
        data_set_kinetic = H5Dcreate(file, "displacement", H5T_NATIVE_FLOAT, data_space_kinetic, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        mem_space_kinetic = H5Screate_simple(rank_field, dims_field, NULL);
    }
    if (globals->diss_avg.Value())
    {
        data_space_1D_avg = H5Screate_simple(rank_1D, max_dims_1D_diss_avg, NULL);
        data_set_1D_avg = H5Dcreate(file, globals->diss_avg.StringID().c_str(), H5T_NATIVE_FLOAT, data_space_1D_avg, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        mem_space_1D_avg = H5Screate_simple(rank_1D, dims_1D_diss_avg, NULL);
    }
    if (globals->kinetic_avg.Value())
    {
        data_space_1D_kinetic_avg = H5Screate_simple(rank_1D, max_dims_1D_diss_avg, NULL);
        data_set_1D_kinetic_avg = H5Dcreate(file, globals->kinetic_avg.StringID().c_str(), H5T_NATIVE_FLOAT, data_space_1D_kinetic_avg, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        mem_space_1D_kinetic_avg = H5Screate_simple(rank_1D, dims_1D_diss_avg, NULL);
    }
    if (globals->field_dummy1_output.Value())
    {
        data_space_dummy1 = H5Screate_simple(rank_field, max_dims_field, NULL); // 2D data space
        data_set_dummy1 = H5Dcreate(file, "displacement", H5T_NATIVE_FLOAT, data_space_dummy1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        mem_space_dummy1 = H5Screate_simple(rank_field, dims_field, NULL);
    }
    if (globals->field_dummy2_output.Value())
    {
        data_space_dummy2 = H5Screate_simple(rank_field, max_dims_field, NULL); // 2D data space
        data_set_dummy2 = H5Dcreate(file, "displacement", H5T_NATIVE_FLOAT, data_space_dummy2, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        mem_space_dummy2 = H5Screate_simple(rank_field, dims_field, NULL);
    }
    //
    // if (consts->sh_coeff_output.Value()) {
    //     data_space_harm_coeff = H5Screate_simple(rank_harm, max_dims_harm_coeff, NULL);
    //     data_set_harm_coeff = H5Dcreate(file, "sh coefficients", H5T_NATIVE_FLOAT, data_space_harm_coeff, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     mem_space_harm_coeff = H5Screate_simple(rank_harm, dims_harm_coeff, NULL);
    // }
    //
    start[0] = 0;
    start[1] = 0;
    //
    count[0] = 1;
    count[1] = node_num;
    // count[2] = eta_cols;
    //
    start_1D[0] = 0;
    count_1D[0] = 1;
    //
    // start_harm[0] = 0;
    // start_harm[1] = 0;
    // start_harm[2] = 0;
    // start_harm[3] = 0;
    //
    // count_harm[0] = 1;
    // count_harm[1] = 2;
    // count_harm[2] = harm_rows;
    // count_harm[3] = harm_cols;
    //
    // //-----------------------Write dataset attributes---------------------------
    //
    // hsize_t dims_coord[1];
    // dims_coord[0] = 2;
    //
    // hid_t attr_space;
    // hid_t attr;
    //
    // int data_set_size[2];
    // float data_set_dx[2];
    //
    //
    // //---------------------------Displacement-----------------------------------
    // if (consts->field_displacement_output.Value()) {
    //     attr_space = H5Screate_simple(rank_size, dims_coord, NULL);
    //
    //     attr = H5Acreate(data_set_eta, "nlat; nlon", H5T_NATIVE_INT, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    //
    //     data_set_size[0] = etaLatLen;
    //     data_set_size[1] = etaLonLen;
    //
    //     H5Awrite(attr, H5T_NATIVE_INT, data_set_size);
    //
    //     attr = H5Acreate(data_set_eta, "dlat; dlon", H5T_NATIVE_FLOAT, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    //
    //     data_set_dx[0] = (float)etadLat*1.0/radConv;
    //     data_set_dx[1] = (float)etadLon*1.0/radConv;
    //
    //     H5Awrite(attr, H5T_NATIVE_FLOAT, data_set_dx);
    // }
    //
    // //---------------------------Velocity---------------------------------------
    // if (consts->field_velocity_output.Value()) {
    //     attr = H5Acreate(data_set_v, "nlat; nlon", H5T_NATIVE_INT, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    //
    //     data_set_size[0] = vLatLen;
    //     data_set_size[1] = vLonLen;
    //
    //     H5Awrite(attr, H5T_NATIVE_INT, data_set_size);
    //
    //     attr = H5Acreate(data_set_v, "dlat; dlon", H5T_NATIVE_FLOAT, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    //
    //     data_set_dx[0] = (float)vdLat*1.0/radConv;
    //     data_set_dx[1] = (float)vdLon*1.0/radConv;
    //
    //     H5Awrite(attr, H5T_NATIVE_FLOAT, data_set_dx);
    //
    //     attr = H5Acreate(data_set_u, "nlat; nlon", H5T_NATIVE_INT, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    //
    //     data_set_size[0] = uLatLen;
    //     data_set_size[1] = uLonLen;
    //
    //     H5Awrite(attr, H5T_NATIVE_INT, data_set_size);
    //
    //     attr = H5Acreate(data_set_u, "dlat; dlon", H5T_NATIVE_FLOAT, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    //
    //     data_set_dx[0] = (float)udLat*1.0/radConv;
    //     data_set_dx[1] = (float)udLon*1.0/radConv;
    //
    //     H5Awrite(attr, H5T_NATIVE_FLOAT, data_set_dx);
    // }
    //
    // //-------------------------Dissipated energy----------------------------------
    // if (consts->field_diss_output.Value()) {
    //     attr = H5Acreate(data_set_diss, "nlat; nlon", H5T_NATIVE_INT, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    //
    //     data_set_size[0] = etaLatLen;
    //     data_set_size[1] = etaLonLen;
    //
    //     H5Awrite(attr, H5T_NATIVE_INT, data_set_size);
    //
    //     attr = H5Acreate(data_set_diss, "dlat; dlon", H5T_NATIVE_FLOAT, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    //
    //     data_set_dx[0] = (float)etadLat*1.0/radConv;
    //     data_set_dx[1] = (float)etadLon*1.0/radConv;
    //
    //     H5Awrite(attr, H5T_NATIVE_FLOAT, data_set_dx);
    // }

};

void OutFiles::DumpData(Globals * globals, int time_level, double ** data)
{
    std::string out = std::to_string(time_level);
    unsigned int i, j;
    double * p;

    start[0] = time_level-1;
    start[1] = 0;

    count[0] = 1;

    start_1D[0] = time_level - 1;
    count_1D[0] = 1;

    for (j=0; j<tags->size(); j++)
    {
        //std::cout<<"HERE: "<<(*tags)[j]<<std::endl;
        p = data[j];
        if ((*tags)[j] == "velocity output")
        {
            for (i=0; i<node_num; i++) {
                u_1D[i] = (float)(*p);
                p++;

                v_1D[i] = (float)(*p);
                p++;
            }

            count[1] = node_num;

            H5Sselect_hyperslab(data_space_u, H5S_SELECT_SET, start, NULL, count, NULL);
            H5Dwrite(data_set_u, H5T_NATIVE_FLOAT, mem_space_u, data_space_u, H5P_DEFAULT, u_1D);

            count[1] = node_num;

            H5Sselect_hyperslab(data_space_v, H5S_SELECT_SET, start, NULL, count, NULL);
            H5Dwrite(data_set_v, H5T_NATIVE_FLOAT, mem_space_v, data_space_v, H5P_DEFAULT, v_1D);
        }

        else if ((*tags)[j] == "displacement output")
        {
            for (i=0; i<node_num; i++) {
                eta_1D[i] = (float)(*p);
                p++;
            }

            count[1] = node_num;

            H5Sselect_hyperslab(data_space_eta, H5S_SELECT_SET, start, NULL, count, NULL);
            H5Dwrite(data_set_eta, H5T_NATIVE_FLOAT, mem_space_eta, data_space_eta, H5P_DEFAULT, eta_1D);
        }

        else if ((*tags)[j] == "pressure output")
        {
            for (i=0; i<node_num; i++) {
                press_1D[i] = (float)(*p);
                p++;
            }

            count[1] = node_num;

            H5Sselect_hyperslab(data_space_press, H5S_SELECT_SET, start, NULL, count, NULL);
            H5Dwrite(data_set_press, H5T_NATIVE_FLOAT, mem_space_press, data_space_press, H5P_DEFAULT, press_1D);
        }

        else if ((*tags)[j] == "kinetic output")
        {
            for (i=0; i<node_num; i++) {
                kinetic_1D[i] = (float)(*p);
                p++;
            }

            count[1] = node_num;

            H5Sselect_hyperslab(data_space_kinetic, H5S_SELECT_SET, start, NULL, count, NULL);
            H5Dwrite(data_set_kinetic, H5T_NATIVE_FLOAT, mem_space_kinetic, data_space_kinetic, H5P_DEFAULT, kinetic_1D);
        }

        else if ((*tags)[j] == "dissipation output")
        {
            for (i=0; i<node_num; i++) {
                diss_1D[i] = (float)(*p);
                p++;
            }

            count[1] = node_num;

            H5Sselect_hyperslab(data_space_diss, H5S_SELECT_SET, start, NULL, count, NULL);
            H5Dwrite(data_set_diss, H5T_NATIVE_FLOAT, mem_space_diss, data_space_diss, H5P_DEFAULT, diss_1D);
        }

        else if ((*tags)[j] == "dissipation avg output")
        {
            diss_avg_1D[0] = (float)(*p);

            H5Sselect_hyperslab(data_space_1D_avg, H5S_SELECT_SET, start_1D, NULL, count_1D, NULL);
            H5Dwrite(data_set_1D_avg, H5T_NATIVE_FLOAT, mem_space_1D_avg, data_space_1D_avg, H5P_DEFAULT, diss_avg_1D);
        }

        else if ((*tags)[j] == "kinetic avg output")
        {
            kinetic_avg_1D[0] = (float)(*p);

            H5Sselect_hyperslab(data_space_1D_kinetic_avg, H5S_SELECT_SET, start_1D, NULL, count_1D, NULL);
            H5Dwrite(data_set_1D_kinetic_avg, H5T_NATIVE_FLOAT, mem_space_1D_kinetic_avg, data_space_1D_kinetic_avg, H5P_DEFAULT, kinetic_avg_1D);
        }
    }

    // // ------------------- Write global average dissipation ----------------------
    // if (consts->diss.Value()) {
    //
    //     H5Sselect_hyperslab(data_space_1D_avg, H5S_SELECT_SET, start_1D, NULL, count_1D, NULL);
    //
    //     H5Dwrite(data_set_1D_avg, H5T_NATIVE_FLOAT, mem_space_1D_avg, data_space_1D_avg, H5P_DEFAULT, diss_avg_1D);
    // }

    //
    // for (int i=0; i<vLatLen; i++) {
    //     for (int j=0; j<vLonLen; j++) {
    //         v_1D[i*vLonLen + j] = (float)vNewArray[i][j];
    //     }
    // }
    // double cd = 0.0;
    // for (int i=0; i<etaLatLen; i++) {
    //     i_h = i*2;
    //     for (int j=0; j<etaLonLen; j++) {
    //         j_h = j*2;
    //         switch (consts->fric_type) {
    //             case QUADRATIC:
    //               cd = cdArray[i_h][j_h];
    //               //std::cout<<cd;
    //               if (cd <= 0.0 || depthArray[i_h][j_h] <= 0.0) diss_1D[i*etaLonLen + j] = 0.0;
    //               else diss_1D[i*etaLonLen + j] = (float)energy->solution[i][j]*2.0*cd/(energy->cellArea[i][j]*depthArray[i_h][j_h]);
    //               break;
    //             case LINEAR:
    //               diss_1D[i*etaLonLen + j] = (float)energy->solution[i][j]*2.0*consts->alpha.Value()/energy->cellArea[i][j];
    //               break;
    //
    //         // if (!(depthArray[i_h][j_h] > 0.0)) diss_1D[i*etaLonLen + j] = 0.0;
    //         //std::cout<<energy->solution[i][j]*cd<<std::endl;
    //         }
    //     }
    // }
    //
    //
    // for (int i=0; i<l_max+1; i++)
    // {
    //     for (int j=0; j<l_max+1; j++)
    //     {
    //         // harm_coeff_1D[i*harm_cols + j] = (float)SH_cos_coeff[i][j];
    //         harm_coeff_1D[i*harm_cols + j] = (float)SH_cos_coeff[i*(l_max+1) + j];
    //     }
    // }
    //
    // for (int i=0; i<l_max+1; i++)
    // {
    //     for (int j=0; j<l_max+1; j++)
    //     {
    //         // harm_coeff_1D[(l_max+1)*(l_max+1) + i*harm_cols + j] = (float)SH_sin_coeff[i][j];
    //         harm_coeff_1D[(l_max+1)*(l_max+1) + i*harm_cols + j] = (float)SH_sin_coeff[i*(l_max+1) + j];
    //     }
    // }
    //
    //
    // diss_avg_1D[0] = energy->currentDissEAvg;
    //
    // start[0] = time_level - 1;
    // start[1] = 0;
    // // start[2] = 0;
    // //
    // count[0] = 1;
    //
    // start_1D[0] = output_num - 1;
    //
    // start_harm[0] = output_num - 1;
    //
    // // count 1D is already set to 1.
    //
    //----------------------------------------------------------------------------
    //------------------------WRITE ARRAYS TO DATA FILE---------------------------
    //----------------------------------------------------------------------------
    //
    //
    // // ----------------------- Write displacement field --------------------------
    //
    // if (consts->field_displacement_output.Value()) {
    //     count[1] = eta_rows;
    //     count[2] = eta_cols;
    //
    //     H5Sselect_hyperslab(data_space_eta, H5S_SELECT_SET, start, NULL, count, NULL);
    //
    //     H5Dwrite(data_set_eta, H5T_NATIVE_FLOAT, mem_space_eta, data_space_eta, H5P_DEFAULT, eta_1D);
    // }
    //
    //
    // ----------------------- Write velocity fields -----------------------------

    // if (globals->field_velocity_output.Value()) {
    //
    //
    //
    //     // ----------------------- Write north velocity field ------------------------
    //
    //     // count[1] = v_rows;
    //     // count[2] = v_cols;
    //     //
    //     // H5Sselect_hyperslab(data_space_v, H5S_SELECT_SET, start, NULL, count, NULL);
    //     //
    //     // H5Dwrite(data_set_v, H5T_NATIVE_FLOAT, mem_space_v, data_space_v, H5P_DEFAULT, v_1D);
    //
    //     // ----------------------- Write east velocity field -------------------------
    //
    //     count[1] = node_num;
    //     // count[2] = u_cols;
    //
    //     H5Sselect_hyperslab(data_space_u, H5S_SELECT_SET, start, NULL, count, NULL);
    //
    //     H5Dwrite(data_set_u, H5T_NATIVE_FLOAT, mem_space_u, data_space_u, H5P_DEFAULT, u_1D);
    // }
    //
    // // ----------------------- Write dissipated energy field ---------------------
    // if (consts->field_diss_output.Value()) {
    //
    //     count[1] = eta_rows;
    //     count[2] = eta_cols;
    //
    //     H5Sselect_hyperslab(data_space_diss, H5S_SELECT_SET, start, NULL, count, NULL);
    //
    //     H5Dwrite(data_set_diss, H5T_NATIVE_FLOAT, mem_space_diss, data_space_diss, H5P_DEFAULT, diss_1D);
    // }
    //
    // // ------------------- Write global average dissipation ----------------------
    // if (consts->diss.Value()) {
    //
    //     H5Sselect_hyperslab(data_space_1D_avg, H5S_SELECT_SET, start_1D, NULL, count_1D, NULL);
    //
    //     H5Dwrite(data_set_1D_avg, H5T_NATIVE_FLOAT, mem_space_1D_avg, data_space_1D_avg, H5P_DEFAULT, diss_avg_1D);
    // }
    //
    // // ----------------------- Write SH coefficients ------------------------
    // if (consts->sh_coeff_output.Value()) {
    //     H5Sselect_hyperslab(data_space_harm_coeff, H5S_SELECT_SET, start_harm, NULL, count_harm, NULL);
    //
    //     H5Dwrite(data_set_harm_coeff, H5T_NATIVE_FLOAT, mem_space_harm_coeff, data_space_harm_coeff, H5P_DEFAULT, harm_coeff_1D);
    // }

};
