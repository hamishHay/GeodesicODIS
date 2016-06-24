#include "field.h"
#include "depth.h"
#include "energy.h"
#include "mass.h"
#include "globals.h"

#include <math.h>

Mass::Mass(Mesh * mesh, int lat, int lon, Globals * Consts, Field * Eta, Depth * Thickness) : Field(mesh, lat, lon) {
	consts = Consts;
	eta = Eta;
	thickness = Thickness;

	UpdateMass();
};

void Mass::UpdateMass(void) {
	for (int i = 0; i < fieldLatLen; i++) {
		for (int j = 0; j < fieldLonLen; j++) {
			UpdateCellMass(i, j);
		}
	}
};

void Mass::UpdateCellMass(int i, int j) {
	solution[i][j] = pow((consts->radius.Value() + thickness->solution[i][j]), 3) - pow(consts->radius.Value(), 3);

	solution[i][j] *= (sin(lat[i]) - sin(lat[i] - dLat));

	solution[i][j] *= (lon[1] - lon[0]);

	solution[i][j] *= 1000. / 3.; //density of water

	//solution[i][j] = 1000*consts->h.Value()*pow(consts->radius.Value(),2)*cos(lat[i])*dLat*dLon;

};

void Mass::UpdateTotalMass(void) {
	totalMass = 0;

	for (int i = 0; i < fieldLatLen; i++) {
		for (int j = 0; j < fieldLonLen; j++) {
			totalMass += solution[i][j];
		}
	}
};
