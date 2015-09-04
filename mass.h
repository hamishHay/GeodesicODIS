#ifndef MASS_H
#define MASS_H

#include "vector"
#include "field.h"
#include "globals.h"
#include "mesh.h"

class Mass : public Field {
private:
	Globals * consts;
	Field * eta;

	double totalMass=0;

	void UpdateCellMass(int i, int j);
	void UpdateTotalMass(void);

public:
	Mass(Mesh *, int, int, Globals * Consts, Field * Eta);

	void UpdateMass();
};

#endif
