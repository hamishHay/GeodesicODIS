#include "vector"
#include "depth.h"
#include "mesh.h"
#include "mathRoutines.h"
#include <math.h>
#include <iostream>

using namespace std;

Depth::Depth(Mesh * mesh) {
  bool latStagg;
	grid = mesh;

	dLat = grid->dLat;
	dLon = grid->dLon;

	fieldLonLen = grid->ReturnLonLen();

  latStagg = false;

  fieldLatLen = grid->ReturnLatLen();

	solution = new double*[fieldLatLen];
	for (int i = 0; i < fieldLatLen; i++)
  {
		solution[i] = new double[fieldLonLen];
		for (int j = 0; j < fieldLonLen; j++)
    {
			solution[i][j] = 0; //Initial Guess
		}
	}

	lat = new double[fieldLatLen];
	for (int i = 0; i < fieldLatLen; i++) lat[i] = grid->lat[i];


	// Make sure last value is 90 exactly
	if (!latStagg) lat[fieldLatLen-1] = -90.0;

	lon = new double[fieldLonLen];
	for (int j = 0; j < fieldLonLen; j++) lon[j] = grid->lon[j];

	cosLat = new double[fieldLatLen];
	sinLat = new double[fieldLatLen];
	cos2Lat = new double[fieldLatLen];
	sin2Lat = new double[fieldLatLen];
	for (int i = 0; i < fieldLatLen; i++) {
		lat[i] *= radConv;
		cosLat[i] = cos(lat[i]);
		sinLat[i] = sin(lat[i]);
		cos2Lat[i] = cos(2*lat[i]);
		sin2Lat[i] = sin(2*lat[i]);
	}
	cosLon = new double[fieldLonLen];
	sinLon = new double[fieldLonLen];
	cos2Lon = new double[fieldLonLen];
	sin2Lon = new double[fieldLonLen];
	for (int j = 0; j < fieldLonLen; j++) {
		lon[j] *= radConv;
		sinLon[j] = sin(lon[j]);
		cosLon[j] = cos(lon[j]);
		sin2Lon[j] = sin(2*lon[j]);
		cos2Lon[j] = cos(2*lon[j]);
	}
	dLat *= radConv;
	dLon *= radConv;
};

int Depth::ReturnFieldLatLen(){
	return fieldLatLen;
};

int Depth::ReturnFieldLonLen(){
	return fieldLonLen;
};

double ** Depth::MakeSolutionArrayCopy(void){
  int i,j;
  double ** solutionCopy;

  solutionCopy = new double*[fieldLatLen];
  solutionCopy[0] = new double[fieldLatLen * fieldLonLen];
  for (i=1; i < fieldLatLen; i++) {
    solutionCopy[i] = &solutionCopy[0][i*fieldLonLen];
  }

  for (i=0; i < fieldLatLen; i++) {
    for (j=0; j < fieldLonLen; j++) {
      solutionCopy[i][j] = 0.0;
    }
  }

  return solutionCopy;
  //solutionReturn = solutionCopy;

};
