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



	void UpdateCellMass(int i, int j);

public:
	double totalMass;
	Mass(Mesh *, int, int, Globals * Consts, Field * Eta, Depth * Thickness);

	void UpdateMass();
    void UpdateTotalMass(void);
};

#endif
