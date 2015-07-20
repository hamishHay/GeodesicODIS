#include "vector"
#include "field.h"
#include "mesh.h"
#include "mathRoutines.h"
#include <iostream>

using namespace std;

Field::Field(Mesh *mesh, int latStagg, int lonStagg) {

	grid = mesh;

	dLat = grid->dLat*2;
	dLon = grid->dLon*2;

	fieldLatLen = 1 + grid->ReturnLatLen() / 2;
	fieldLonLen = grid->ReturnLonLen() / 2;

	if (latStagg) fieldLatLen--;
	//if (lonStagg) fieldLonLen--;
	

	solution.resize(fieldLatLen);
	for (int i = 0; i < fieldLatLen; i++) {
		solution[i].resize(fieldLonLen);
		for (int j = 0; j < fieldLonLen; j++) {
			solution[i][j] = 0; //Initial Guess
		}
	}

	for (int i = 0; i < fieldLatLen; i++) {
		if (latStagg) lat.push_back(grid->lat[i * 2 + 1]);
		else lat.push_back(grid->lat[i * 2]);
	}

	for (int j = 0; j < fieldLonLen; j++) {
		if (lonStagg) lon.push_back(grid->lon[j * 2 + 1]);
		else lon.push_back(grid->lon[j * 2]);
	}

	//populate opp pointer
	//fieldLonLen must be an even number
	opp.resize(fieldLonLen);
	for (int i = 0; i < fieldLonLen; i++) {
		if (i < fieldLonLen / 2) opp[i] = fieldLonLen / 2 + i;
		else opp[i] = i - fieldLonLen / 2; 
	}
};

int Field::ReturnFieldLatLen(){
	return fieldLatLen;
};

int Field::ReturnFieldLonLen(){
	return fieldLonLen;
};

double Field::CenterP(int lat, int lon){
	if (lat < 0) return NorthP(lat, lon);
	if (lat >= this->fieldLatLen) return SouthP(lat, lon);
	return solution[lat][lon];
};

double Field::EastP(int lat, int lon){
	if (lon >= this->ReturnFieldLonLen() - 1) return solution[lat][lon - (this->fieldLonLen - 1)];
	return solution[lat][lon + 1];
};

double Field::NorthP(int lat, int lon){
	if (lat <= 0) {
		if (lon + (this->fieldLonLen / 2) >= this->fieldLonLen) {
			return solution[0 - lat][(lon + this->fieldLonLen / 2) - this->fieldLonLen];
		}
		else {
			return solution[0 - lat][lon + this->fieldLonLen / 2];
		}
	}
	else if (lon < 0) return solution[lat - 1][(this->fieldLonLen - 1)];
	else if (lon > this->fieldLonLen - 1) return solution[lat - 1][0];
	else return solution[lat - 1][lon];
};

double Field::NorthEastP(int lat, int lon){ 
	if (lat == 0) return solution[lat][this->ReturnFieldLonLen() / 4];
	else if (lon == this->ReturnFieldLonLen()-1) return solution[lat - 1][0];
	return solution[lat - 1][lon + 1];
}

double Field::WestP(int lat, int lon){
	if (lat < 0) {
		if (lon + (this->fieldLonLen / 2) >= this->fieldLonLen) {
			return EastP(0 - lat,(lon + this->fieldLonLen / 2) - this->fieldLonLen);
		}
		else {
			return EastP(0 - lat,lon + this->fieldLonLen / 2);
		}
	}
	else if (lon <= 0) return solution[lat][(this->fieldLonLen-1)+lon];
	else return solution[lat][lon - 1];
};

double Field::SouthP(int lat, int lon){
	if (lat >= fieldLatLen - 1) {
		return solution[2 * (this->fieldLatLen - 1) - lat][this->opp[lon]];
	}
	else if (lat < 0) {
		if (lon + (this->fieldLonLen / 2) >= this->fieldLonLen) return solution[0 - lat][(lon + this->fieldLonLen / 2) - this->fieldLonLen];
		else {
			return solution[0 - lat][lon + this->fieldLonLen / 2];
		}
	}
	else if (lon < 0) {
		return solution[lat + 1][(this->fieldLonLen - 1)];
	}
	else if (lon > (this->fieldLonLen - 1)) {
		return solution[lat + 1][0];
	}
	return solution[lat + 1][lon];
};

double Field::SouthWestP(int i, int j){
	if (j == 0) return solution[i+1][fieldLonLen-1];
	else return solution[i + 1][j - 1];
};

double Field::SouthWestAvg(int i, int j) {
	//if (i == 0) {
	//	if (j == 0) {
	//		return (linearInterp2(this, i, this->fieldLonLen - 1) + linearInterp2(this, i, j) + this->SouthP(i, j) + this->SouthP(i, j - 1))*0.25;
	//	}
	//	/*else if (j == this->fieldLonLen-1) {
	//		return (linearInterp2(this, i, this->fieldLonLen - 1) + linearInterp2(this, i, j) + this->SouthP(i, j) + this->SouthP(i, j - 1))*0.25;
	//	}*/
	//	else return (linearInterp2(this, i, j-1) + linearInterp2(this, i, j) + this->SouthP(i, j) + this->SouthP(i, j - 1))*0.25;
	//}
	//else if (i == this->fieldLatLen - 2) {
	//	if (j == 0) {
	//		return (linearInterp2(this, i+1, this->fieldLonLen - 1) + linearInterp2(this, i+1, j) + this->CenterP(i, j) + this->WestP(i, j - 1))*0.25;
	//	}
	//	else return (linearInterp2(this, i+1, j - 1) + linearInterp2(this, i+1, j) + this->CenterP(i, j) + this->WestP(i, j))*0.25;
	//}
	//else 
	return (this->CenterP(i, j) + this->SouthP(i, j) + this->WestP(i, j) + this->SouthP(i, j - 1))*0.25;
};

double Field::NorthEastAvg(int i, int j) {
	if (i == 0) return (this->CenterP(i, j) + this->EastP(i, j))*0.5; //(this->CenterP(i, j) + this->EastP(i, j) + this->CenterP(i, this->opp[j]) + this->CenterP(i, this->opp[(int)j/2]))*0.25;
	else if (i == this->fieldLatLen) return (this->CenterP(i-1, j) + this->EastP(i-1, j))*0.25;
	else return (this->CenterP(i, j) + this->EastP(i, j) + this->NorthP(i, j) + this->NorthP(i, j + 1))*0.25;
};

