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

void Mass::UpdateMass() {

	for (int i = 0; i < fieldLatLen; i++) {		
		for (int j = 0; j < fieldLonLen; j++) {
			UpdateCellMass(i, j);
		}
	}
};

void Mass::UpdateCellMass(int i, int j) {
	solution[i][j] = 1000 * (consts->h.Value() + eta->solution[i][j]) *pow(consts->radius.Value(), 2)*cos(lat[i] * radConv)*dLon*radConv*dLat*radConv;
};

void Mass::UpdateTotalMass(void) {
	totalMass = 0;
	
	//Sum over poles only once!
	//totalMass += solution[0][0];
	//totalMass += solution[fieldLatLen-1][0];
	
	for (int i = 0; i < fieldLatLen; i++) {
		for (int j = 0; j < fieldLonLen; j++) {
			totalMass += solution[i][j];
		}
	}
};

