#ifndef VERBOSE_H
#define VERBOSE_H

#include <iostream>

static void TerminateODIS(void) {
	std::cout << "Terminating ODIS. Press enter to continue." << std::endl;
	getchar();
	std::exit(0);
}

#endif