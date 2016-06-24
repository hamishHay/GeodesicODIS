#ifndef MASS_H
#define MASS_H

#include "vector"
#include "field.h"
#include "depth.h"
#include "globals.h"
#include "mesh.h"

class Mass : public Field {
private:
	Globals * consts;
	Field * eta;
	Depth * thickness;

	double totalMass=0;

	void UpdateCellMass(int i, int j);
	void UpdateTotalMass(void);

public:
	Mass(Mesh *, int, int, Globals * Consts, Field * Eta, Depth * Thickness);

	void UpdateMass();
};

#endif
