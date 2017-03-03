// FILE: field.h
// DESCRIPTION: header file for the field class. The term field is used simply
//              to refer to a solvable quantity that is found over the numerical
//              domain. Field will generate a 2D array to store some qunatity
//              at positions that are staggered over the 2D grid, based on a mesh
//              type object. Note that field quantities are stored at half the
//              grid spacing defined in mesh.
//
//   Version              Date                   Programmer
//     0.1      |        13/09/2016        |        H. Hay        |
// ---- Initial version of ODIS, described in Hay and Matsuyama (2016)


#ifndef FIELD_H
#define FIELD_H

#include "mesh.h"

class Field {
private:

public:
  // Number of nodes in latitude and longitude
  int fieldLatLen;
  int fieldLonLen;

  // 1D arrays of positions in latitude and longitude
  double * lat;
  double * lon;

  // 1D arrays trigonometric quantities that only need to be calculated once.
  // These are essentially lookup tables calculated upon the field initialisation.
  double * cosLat; //cos(latitude)
  double * cosSqLat; //cos^2(lat)
  double * cosCubLat; //cos^3(lat)
  double * sinLat;
  double * sinSqLat;
  double * cos2Lat; //cos(2*latitude)
  double * sin2Lat;
  double * cosLon;
  double * cosSqLon;
  double * cosCubLon;
  double * sinLon;
  double * cos2Lon;
  double * sin2Lon;
  double * sinSqLon;
  double * cos3Lon;
  double * cosCoLat;
  double ** cosMLon;
  double ** sinMLon;

  double *** weights;

  // node spacing
  double dLat;
  double dLon;

  // mesh object used to base staggered grid around.
  Mesh * grid;

  // 2D array for the numerical solution of the field
  double ** solution;

  // Constructor. First int is for latitude staggering and second is for longitude
  // staggering (1 or 0).
  Field(Mesh *, int, int);

  int ReturnFieldLatLen();
  int ReturnFieldLonLen();
  void CalcWeights(int, int, Field *);

  double ** MakeSolutionArrayCopy(void);    // Sets up double ** to be a copy of solution
};

#endif
