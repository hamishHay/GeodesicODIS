#ifndef MASS_H
#define MASS_H

#include "vector"
#include "field.h"
#include "globals.h"
#include "mesh.h"

class Mass : public Field {
private:
	double totalMass=0;
	void UpdateCellMass(int i, int j);
public:
	Mass(Mesh *, int, int, Globals * Consts, Field * Eta);

	Globals * consts;
	Field * eta;

	void UpdateTotalMass(void);
	void UpdateMass();
};

#endif