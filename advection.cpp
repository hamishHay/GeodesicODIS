#include "globals.h"
#include "field.h"

#include <math.h>
#include <algorithm>    // std::max
#include <iostream>

void UpdateAdvection(Field * uField, Field * vField, Field * etaField, Globals * constants) {
  double ** u, ** v, ** eta, ** uAdv, ** vAdv, ** etaAdv;
  double ** vInterp, ** uInterp;
  int ** uMask, ** vMask, ** etaMask;
  double * vCosLat, * uCosLat, * etaCosLat;
  double dLat, dLon;
  int uLonLen, uLatLen, vLonLen, vLatLen, etaLonLen, etaLatLen;
  double R, inv_R;
  double uPos, uNeg, vPos, vNeg;
  double dU_forward, dU_backward;
  // double north, south;
  double cosCoeff;
  int i, j;

  // DEFINE CONSTANTS ----------------------------------------------------------

  u = uField->solution;
  v = vField->solution;
  eta = etaField->solution;

  uAdv = uField->advSolution;
  vAdv = vField->advSolution;
  etaAdv = etaField->advSolution;

  vInterp = uField->vInterp;            // v velocity interpolated to u nodes
  uInterp = vField->uInterp;            // u velocity interpolated to v nodes

  uMask = uField->mask;
  vMask = vField->mask;
  etaMask = etaField->mask;

  uCosLat = uField->cosLat;
  vCosLat = vField->cosLat;
  etaCosLat = etaField->cosLat;

  dLat = uField->dLat;
  dLon = uField->dLon;

  vLonLen = vField->fieldLonLen;
  vLatLen = vField->fieldLatLen;

  uLonLen = uField->fieldLonLen;
  uLatLen = uField->fieldLatLen;

  etaLonLen = etaField->fieldLonLen;
  etaLatLen = etaField->fieldLatLen;

  R = constants->radius.Value();
  inv_R = 1.0/R;

  // CALCULATE ADVECTION OF EASTWARD VELOCITY ----------------------------------

  for (i=1; i<uLatLen-1; i++)
  {
    cosCoeff = inv_R * 1.0/(uCosLat[i]*dLat); //assume dLat = dLon
    for (j=0; j<uLonLen; j++)
    {
      switch (uMask[i][j])
      {
        case 1:
          // COMPUTE LONGITUDINAL DERIVATIVE
          uPos = std::max(u[i][j], 0.0);
          uNeg = std::min(u[i][j], 0.0);

          if (j != 0 && j != uLonLen-1)
          {
            dU_forward = u[i][j+1] - u[i][j];
            dU_backward = u[i][j] - u[i][j-1];
          }
          else if (j == 0)
          {
            dU_forward = u[i][j+1] - u[i][j];
            dU_backward = u[i][j] - u[i][uLonLen-1];
          }
          else
          {
            dU_forward = u[i][0] - u[i][j];
            dU_backward = u[i][j] - u[i][j-1];
          }

          uAdv[i][j] = cosCoeff * (uPos*dU_backward + uNeg*dU_forward);

          // COMPUTE LATITUDINAL DERIVATIVE

          vPos = std::max(vInterp[i][j], 0.0);
          vNeg = std::min(vInterp[i][j], 0.0);

          dU_forward = u[i][j] - u[i+1][j];
          dU_backward = u[i-1][j] - u[i][j];

          // Back and forward derivatives swapped as coord system and index
          // are opposite.
          uAdv[i][j] += inv_R * (vPos*dU_forward + vNeg*dU_backward);
          break;

        default:
          break;
      }
    }
  }



  // CALCULATE ADVECTION OF NORTHWARD VELOCITY ---------------------------------

  for (i=0; i<vLatLen; i++)
  {
    cosCoeff = inv_R * 1.0/(vCosLat[i]*dLat); //assume dLat = dLon
    for (j=0; j<vLonLen; j++)
    {
      switch (vMask[i][j])
      {
        case 1:
          // COMPUTE LONGITUDINAL DERIVATIVE
          uPos = std::max(uInterp[i][j], 0.0);
          uNeg = std::min(uInterp[i][j], 0.0);

          if ((j != 0) && (j != vLonLen-1))
          {
            dU_forward = v[i][j+1] - v[i][j];
            dU_backward = v[i][j] - v[i][j-1];
          }
          else if (j == 0)
          {
            dU_forward = v[i][j+1] - v[i][0];
            dU_backward = v[i][0] - v[i][vLonLen-1];
          }
          else if (j == vLonLen-1)
          {
            dU_forward = v[i][0] - v[i][vLonLen-1];
            dU_backward = v[i][j] - v[i][j-1];
          }

          vAdv[i][j] = cosCoeff * (uPos*dU_backward + uNeg*dU_forward);

          // COMPUTE LATITUDINAL DERIVATIVE

          vPos = std::max(v[i][j], 0.0);
          vNeg = std::min(v[i][j], 0.0);

          dU_forward = v[i+1][j] - v[i][j];
          dU_backward = v[i][j] - v[i-1][j];

          // Back and forward derivatives swapped as coord system and index
          // are opposite.
          vAdv[i][j] += inv_R * (vPos*dU_backward + vNeg*dU_forward);
          break;

        default:
          break;
      }
    }
  }

};
