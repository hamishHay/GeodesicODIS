#include "field.h"
#include "energy.h"
#include "mass.h"
#include "globals.h"

#include <math.h>

Mass::Mass(Mesh * mesh, int lat, int lon, Globals * Consts, Field * Eta) : Field(mesh, lat, lon) {
	consts = Consts;
	eta = Eta;

	UpdateMass();
};

void Mass::UpdateMass(void) {

	for (int i = 0; i < fieldLatLen-1; i++) {		
		for (int j = 0; j < fieldLonLen; j++) {
			UpdateCellMass(i, j);
		}
	}

	//UpdateTotalMass();
};

void Mass::UpdateCellMass(int i, int j) {

	solution[i][j] = pow((consts->radius.Value() + consts->h.Value()), 3) - pow(consts->radius.Value(), 3);

	solution[i][j] *= (sin(lat[i] * radConv) - sin((lat[i] - dLat)*radConv));

	solution[i][j] *= (lon[1] - lon[0])*radConv; 

	solution[i][j] *= 1000. / 3.; //density of water
};

void Mass::UpdateTotalMass(void) {
	totalMass = 0;
	
	for (int i = 0; i < fieldLatLen; i++) {
		for (int j = 0; j < fieldLonLen; j++) {
			totalMass += solution[i][j];
		}
	}
};

