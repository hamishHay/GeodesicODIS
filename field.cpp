// FILE: field.cpp
// DESCRIPTION: main file for the field class. The term field is used simply to
//              refer to a solvable quantity that is found over the numerical
//              domain. Field will generate a 2D array to store some qunatity
//              at positions that are staggered over the 2D grid, based on a mesh
//              type object. Note that field quantities are stored at half the
//              grid spacing defined in mesh.
//
//   Version              Date                   Programmer
//     0.1      |        13/09/2016        |        H. Hay        |
// ---- Initial version of ODIS, described in Hay and Matsuyama (2016)

#include "field.h"
#include "mesh.h"
#include "mathRoutines.h"
#include <math.h>
#include <iostream>

using namespace std;

Field::Field(Mesh *mesh, int latStagg, int lonStagg)
{
  // Constructor precalculates all lookup tables and finds the appropriate spatial
  // positions for the staggered field.
  // inputs: mesh: used find grid spacing/number of nodes
  //         latStagg: stagger latitde by one node south (1 yes, 0 no)
  //         lonStagg: stagger longitude by one node east (1 yes, 0 no)


  grid = mesh;

  dLat = grid->dLat*2; // staggered grid is half resolution of original grid.
  dLon = grid->dLon*2;


  // define staggered grid dimensions
  fieldLonLen = grid->ReturnLonLen() / 2;

  if (latStagg) fieldLatLen = (grid->ReturnLatLen() - 1 )/ 2;
  else fieldLatLen = (1 + grid->ReturnLatLen())/ 2;

  solution = new double*[fieldLatLen];
  solution[0] = new double[fieldLatLen * fieldLonLen];
  for (int i=1; i < fieldLatLen; i++) {
    solution[i] = &solution[0][i*fieldLonLen];
  }

 for (int i=0; i < fieldLatLen; i++) {
   for (int j=0; j < fieldLonLen; j++) {
     solution[i][j] = 0.0;
   }
 }

  // allocate and assign position values for lat and lon.
  lat = new double[fieldLatLen];
  for (int i = 0; i < fieldLatLen; i++) {
    if (latStagg) lat[i] = grid->lat[i * 2 + 1];
    else lat[i] = grid->lat[i * 2];
  }

  lon = new double[fieldLonLen];
  for (int j = 0; j < fieldLonLen; j++) {
    if (lonStagg) lon[j] = grid->lon[j * 2 + 1];
    else lon[j] = grid->lon[j * 2];
  }

  // allocate and calculate lookup tables for trig functions of position
  cosLat = new double[fieldLatLen];
  sinLat = new double[fieldLatLen];
  cos2Lat = new double[fieldLatLen];
  sin2Lat = new double[fieldLatLen];
  cosCoLat = new double[fieldLatLen];
  cosSqLat = new double[fieldLatLen];
  cosCubLat = new double[fieldLatLen];
  sinSqLat = new double[fieldLatLen];

  double colat = 0.0;
  for (int i = 0; i < fieldLatLen; i++)
  {
    colat = 90.0 - lat[i];
    lat[i] *= radConv;
    cosLat[i] = cos(lat[i]);
    sinLat[i] = sin(lat[i]);
    cos2Lat[i] = cos(2*lat[i]);
    sin2Lat[i] = sin(2*lat[i]);

    cosCoLat[i] = cos(pi*0.5 - lat[i]);
    cosSqLat[i] = pow(cos(lat[i]),2.0);
    cosCubLat[i] = pow(cos(lat[i]),3.0);

    sinSqLat[i] = pow(sin(lat[i]),2.0);

    cosCoLat[i] = cos(colat*radConv);
  }

  cosLon = new double[fieldLonLen];
  sinLon = new double[fieldLonLen];
  cos2Lon = new double[fieldLonLen];
  sin2Lon = new double[fieldLonLen];
  cos3Lon = new double[fieldLonLen];
  cosSqLon = new double[fieldLonLen];
  cosCubLon = new double[fieldLonLen];
  sinSqLon = new double[fieldLonLen];
  for (int j = 0; j < fieldLonLen; j++)
  {
    lon[j] *= radConv;
    sinLon[j] = sin(lon[j]);
    cosLon[j] = cos(lon[j]);
    sin2Lon[j] = sin(2.*lon[j]);
    cos2Lon[j] = cos(2.*lon[j]);
    cos3Lon[j] = cos(3.*lon[j]);
    cosSqLon[j] = pow(cos(lon[j]),2.0);
    cosCubLon[j] = pow(cos(lon[j]),3.0);
    sinSqLon[j] = pow(sin(lon[j]),2.0);
  }

  cosMLon = new double * [mesh->globals->l_max.Value()+1];
  sinMLon = new double * [mesh->globals->l_max.Value()+1];

  for (int m = 0; m <= mesh->globals->l_max.Value(); m++) {
    cosMLon[m] = new double[fieldLonLen];
    sinMLon[m] = new double[fieldLonLen];
    for (int j = 0; j < fieldLonLen; j++) {
      cosMLon[m][j] = cos(m*lon[j]);
      sinMLon[m][j] = sin(m*lon[j]);
    }
  }

  dLat *= radConv;
  dLon *= radConv;
};


double ** Field::MakeSolutionArrayCopy(void){
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

void Field::CalcWeights(int latStagg, int lonStagg, Field * neighbour) {
    int i,j,k,l;
    double lat1,lat2,lon1,lon2;
    double totalLength, length;

    weights = new double **[fieldLatLen];
    for (i=0; i<fieldLatLen; i++) {
        weights[i] = new double *[fieldLonLen];
        for (k=0; k<2; k++) {
            weights[i][k] = new double [2];
        }
    }


    for (i=0; i<fieldLatLen; i++) {
        totalLength = 0.0;
        j = 2;
        lat1 = lat[i];
        lon1 = lon[j];

        if (latStagg) {
            for (k=0; k<2; k++) {
                for (l=0; l<2; l++) {
                    lat2 = neighbour->lat[i + k];
                    lon2 = neighbour->lon[j + (l-1)];
                    totalLength += 1./arcLength(lat1, lat2, lon1, lon2);
                }
            }

            for (k=0; k<2; k++) {
                for (l=0; l<2; l++) {
                    lat2 = neighbour->lat[i+k];
                    // std::cout<<lat2/radConv<<std::endl;
                    lon2 = neighbour->lon[j+(l-1)];
                    length = 1./arcLength(lat1, lat2, lon1, lon2);
                    weights[i][k][l] = length/totalLength;
                }
            }

        }
        else {
            if (i > 0 && i < fieldLatLen-1) {
                for (k=0; k<2; k++) {
                    for (l=0; l<2; l++) {
                        lat2 = neighbour->lat[i + k];
                        lon2 = neighbour->lon[j + (l-1)];
                        totalLength += 1./arcLength(lat1, lat2, lon1, lon2);
                    }
                }

                for (k=0; k<2; k++) {
                    for (l=0; l<2; l++) {
                        lat2 = neighbour->lat[i+k];
                        // std::cout<<lat2/radConv<<std::endl;
                        lon2 = neighbour->lon[j+(l-1)];
                        length = 1./arcLength(lat1, lat2, lon1, lon2);
                        weights[i][k][l] = length/totalLength;
                    }
                }
            }
            else if (i == 0) {
                for (k=0; k<1; k++) {
                    for (l=0; l<2; l++) {
                        lat2 = neighbour->lat[i + k];
                        lon2 = neighbour->lon[j + (l-1)];
                        totalLength += 1./arcLength(lat1, lat2, lon1, lon2);
                    }
                }

                for (k=0; k<1; k++) {
                    for (l=0; l<2; l++) {
                        lat2 = neighbour->lat[i+k];
                        // std::cout<<lat2/radConv<<std::endl;
                        lon2 = neighbour->lon[j+(l-1)];
                        length = 1./arcLength(lat1, lat2, lon1, lon2);
                        weights[i][k][l] = length/totalLength;
                    }
                }
            }
            else if (i == fieldLatLen-1) {
                for (k=0; k<1; k++) {
                    for (l=0; l<2; l++) {
                        lat2 = neighbour->lat[i - k];
                        lon2 = neighbour->lon[j + (l-1)];
                        totalLength += 1./arcLength(lat1, lat2, lon1, lon2);
                    }
                }

                for (k=0; k<1; k++) {
                    for (l=0; l<2; l++) {
                        lat2 = neighbour->lat[i-k];
                        // std::cout<<lat2/radConv<<std::endl;
                        lon2 = neighbour->lon[j+(l-1)];
                        length = 1./arcLength(lat1, lat2, lon1, lon2);
                        weights[i][k][l] = length/totalLength;
                    }
                }
            }
        }
    }
};

int Field::ReturnFieldLatLen(){
  return fieldLatLen;
};

int Field::ReturnFieldLonLen(){
  return fieldLonLen;
};
