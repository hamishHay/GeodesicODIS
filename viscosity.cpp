#include "globals.h"
#include "field.h"

#include <math.h>

void UpdateViscosity(Field * uField, Field * vField, Globals * constants) {
  double ** u, ** v, ** uVisc, ** vVisc;
  int ** uMask, ** vMask;
  double * vCosLat, * vCosSqLat;
  double * uCosLat, * uCosSqLat;
  double dLat, vdLon, vdLat;
  int uLonLen, uLatLen, vLonLen, vLatLen;
  double viscosity, R, inv_rdlat;
  double north, south;
  double cosCoeff1, cosCoeff2;
  int i, j;

  // DEFINE CONSTANTS ----------------------------------------------------------

  u = uField->solution;
  v = vField->solution;

  uVisc = uField->viscSolution;
  vVisc = vField->viscSolution;

  uMask = uField->mask;
  vMask = vField->mask;

  uCosLat = uField->cosLat;
  uCosSqLat = uField->cosSqLat;

  vCosLat = vField->cosLat;
  vCosSqLat = vField->cosSqLat;

  dLat = uField->dLat;

  vLonLen = vField->fieldLonLen;
  vLatLen = vField->fieldLatLen;

  uLonLen = uField->fieldLonLen;
  uLatLen = uField->fieldLatLen;

  viscosity = 5e5; // Pa s
  R = constants->radius.Value();

  inv_rdlat = viscosity/(pow(R,2.0) * pow(dLat, 2.0)); // assuming dLat == dLon


  // SOLVE EAST VELOCITY -------------------------------------------------------

  for (i = 1; i < uLatLen-1; i++)
  {
    cosCoeff1 = inv_rdlat * 1./uCosLat[i];
    cosCoeff2 = inv_rdlat * 1./uCosSqLat[i];
    for (j = 0; j < uLonLen; j++)
    {
      switch (uMask[i][j])
      {
        case 1: // If cell is wet and not a boundary cell:
          north = (u[i-1][j] - u[i][j])*vCosLat[i-1];
          south = (u[i][j] - u[i+1][j])*vCosLat[i];

          uVisc[i][j] = cosCoeff1 * (north - south);

          if (j != 0 && j != uLonLen-1)           // Interior cell (most likely)
          {
            uVisc[i][j] += cosCoeff2 * (u[i][j+1] - 2.*u[i][j] + u[i][j-1]);
          }
          else if (j == 0)                        // western wrap around
          {
            uVisc[i][j] += cosCoeff2 * (u[i][j+1] - 2.*u[i][j] + u[i][uLonLen-1]);
          }
          else if (j == uLonLen-1)                // eastern wrap around
          {
            uVisc[i][j] += cosCoeff2 * (u[i][0] - 2.*u[i][j] + u[i][j-1]);
          }
          break;

        default: // Dry cell: do nothing
          break;
      }
    }
  }

  // SOLVE NORTH VELOCITY ------------------------------------------------------

  for (i = 1; i < vLatLen-1; i++)
  {
    cosCoeff1 = inv_rdlat * 1./vCosLat[i];
    cosCoeff2 = inv_rdlat * 1./vCosSqLat[i];
    for (j = 0; j < vLonLen; j++)
    {
      switch (vMask[i][j])
      {
        case 1:  // If cell is wet and not a boundary cell:
          north = (v[i-1][j] - v[i][j])*uCosLat[i-1];
          south = (v[i][j] - v[i+1][j])*uCosLat[i];

          vVisc[i][j] = cosCoeff1 * (north - south);

          if ((j != 0) && (j != vLonLen - 1))     // Interior cell (most likely)
          {
            vVisc[i][j] += cosCoeff2 * (v[i][j+1] - 2.*v[i][j] + v[i][j-1]);
          }
          else if (j == 0)                        // western wrap around
          {
            vVisc[i][j] += cosCoeff2 * (v[i][j+1] - 2.*v[i][j] + v[i][vLonLen-1]);
          }

          else if (j==vLonLen-1)                  // eastern wrap around
          {
            vVisc[i][j] += cosCoeff2 * (v[i][0] - 2.*v[i][j] + v[i][j-1]);
          }

          break;

        default:
          break;

      }
    }
  }

};
