#ifndef OUTFILES_H
#define OUTFILES_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "globals.h"

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

};

#endif
