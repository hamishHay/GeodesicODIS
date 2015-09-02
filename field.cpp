#include "vector"
#include "field.h"
#include "mesh.h"
#include "mathRoutines.h"
#include <math.h>
#include <iostream>

using namespace std;

Field::Field(Mesh *mesh, int latStagg, int lonStagg) {

	grid = mesh;

	dLat = grid->dLat*2;
	dLon = grid->dLon*2;

	//fieldLatLen = grid->ReturnLatLen() / 2.;
	fieldLonLen = grid->ReturnLonLen() / 2;

	if (latStagg) fieldLatLen = (grid->ReturnLatLen() - 1 )/ 2;
	else fieldLatLen = 1 + (grid->ReturnLatLen() - 1)/ 2;
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

	// Make sure last value is 90 exactly
	if (!latStagg) lat[lat.size()-1] = -90.0;

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

	cosLat.resize(fieldLatLen);
	sinLat.resize(fieldLatLen);
	cos2Lat.resize(fieldLatLen);
	sin2Lat.resize(fieldLatLen);
	for (int i = 0; i < this->fieldLatLen; i++) {
		lat[i] *= radConv;
		cosLat[i] = cos(lat[i]);
		sinLat[i] = sin(lat[i]);
		cos2Lat[i] = cos(2*lat[i]);
		sin2Lat[i] = sin(2*lat[i]);
	}
	for (int j = 0; j < this->fieldLonLen; j++) lon[j] *= radConv;

	dLat *= radConv;
	dLon *= radConv;

};

int Field::ReturnFieldLatLen(){
	return fieldLatLen;
};

int Field::ReturnFieldLonLen(){
	return fieldLonLen;
};

double Field::CenterP(int lat, int lon){
	return solution[lat][lon];
};

double Field::EastP(int lat, int lon){
	if (lon >= fieldLonLen - 1) return solution[lat][lon - (fieldLonLen-1)];
	else return solution[lat][lon + 1];
};

double Field::NorthP(int lat, int lon){
	return solution[lat - 1][lon];
};

double Field::NorthEastP(int lat, int lon){ 
	if (lon >= fieldLonLen - 1) return solution[lat - 1][lon - (fieldLonLen - 1)];
	else return solution[lat - 1][lon + 1];
}

double Field::WestP(int lat, int lon){
	if (lon <= 0) return solution[lat][(fieldLonLen - 1) - lon];
	else return solution[lat][lon - 1];
};

double Field::SouthP(int lat, int lon){
	return solution[lat + 1][lon];
};

double Field::SouthWestP(int lat, int lon){
	if (lon <= 0) return solution[lat+1][(fieldLonLen-1) - lon];
	else return solution[lat + 1][lon - 1];
};

double Field::SouthWestAvg(int i, int j) {
	if (j > 0) {
		return 0.25*(solution[i][j] + solution[i + 1][j] + solution[i][j - 1] + solution[i + 1][j - 1]);
	}
	else {
		return 0.25*(solution[i][j] + solution[i + 1][j] + solution[i][fieldLonLen - 1] + solution[i + 1][fieldLonLen-1]);
	}

	//return (this->CenterP(i, j) + this->SouthP(i, j) + this->WestP(i, j) + this->SouthWestP(i, j))*0.25;
};

double Field::NorthEastAvg(int i, int j) {
	if (j < fieldLonLen-1) {
		return 0.25*(solution[i][j] + solution[i - 1][j] + solution[i][j + 1] + solution[i - 1][j + 1]);
	}
	else {
		return 0.25*(solution[i][j] + solution[i - 1][j] + solution[i][0] + solution[i - 1][0]);
	}
	//return (this->CenterP(i, j) + this->EastP(i, j) + this->NorthP(i, j) + this->NorthEastP(i, j))*0.25;
};

