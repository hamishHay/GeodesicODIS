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
#include <sstream>
#include <string>
#include <vector>
#include "H5Cpp.h"

OutFiles::OutFiles() {
	char buffer[PATH];

#ifdef _WIN64
	GetModuleFileName(NULL, buffer, sizeof(buffer)-1);
#elif __linux__
	readlink("/proc/self/exe", buffer, sizeof(buffer)-1);
#endif

	std::string::size_type pos = std::string(buffer).find_last_of(SEP);
	path = std::string(buffer).substr(0, pos);


	outName = path + SEP + "OUTPUT.txt";
	errName = path + SEP + "ERROR.txt";
	dataPath = path + SEP + "DATA" + SEP + "data.h5";

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
	sstream << "        | ODIS is a numerical finite difference ocean dynamics   |" << std::endl;
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

// int OutFiles::DumpData(double ** data_pointer, std::vector<std::string> * tags_pointer)
// {
//
// };
//
// int OutFiles::CreateHDF5FrameWork(void)
// {
//
// };
