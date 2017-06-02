#include "field.h"

#include <math.h>

void InterpolateVel2Vel(Field * uField, Field * vField) {
  double ** u, ** v;
  double ** vInterp, ** uInterp;
  int ** uMask, ** vMask;
  int uLatLen, uLonLen;
  int vLatLen, vLonLen;
  int i, j;

  // POINT TO DESIRED ARRAYS AND SET CONSTANTS ---------------------------------

  u = uField->solution;
  v = vField->solution;

  vInterp = uField->vInterp;              // array for interpolated v at u nodes
  uInterp = vField->uInterp;              // array for interpolated u at v nodes

  uMask = uField->mask;
  vMask = vField->mask;

  uLatLen = uField->fieldLatLen;
  uLonLen = uField->fieldLonLen;

  vLatLen = vField->fieldLatLen;
  vLonLen = vField->fieldLonLen;


  // INTERPOLATE NORTH VELOCITY TO EAST VELOCITY NODE LOCATIONS ----------------

  for (i=1; i<uLatLen-1; i++)
  {
    for (j=0; j<uLonLen; j++)
    {
      switch (uMask[i][j])
      {
        case 1:
          if (j != uLonLen-1) vInterp[i][j] = 0.25*(v[i][j] + v[i-1][j] + v[i-1][j+1] + v[i][j+1]);
          else vInterp[i][uLonLen-1] = 0.25*(v[i][uLonLen-1] + v[i-1][uLonLen-1] + v[i-1][0] + v[i][0]);
          break;
        default:
          break;
      }
    }
  }


  // INTERPOLATE EAST VELOCITY TO NORTH VELOCITY NODE LOCATIONS ----------------

  for (i=0; i<vLatLen; i++)
  {
    for (j=0; j<vLonLen; j++)
    {
      switch (vMask[i][j])
      {
        case 1:
           if (j != 0) uInterp[i][j] = 0.25*(u[i][j] + u[i+1][j] + u[i][j-1] + u[i+1][j-1]);
           else uInterp[i][0] = 0.25*(u[i][0] + u[i+1][0] + u[i][uLonLen-1] + u[i+1][uLonLen-1]);
          break;
        default:
          break;
      }
    }
  }
};

void InterpolateDisp2Vel(Field * uField, Field * vField, Field * etaField) {
  double ** eta;
  double ** etaInterp;
  int ** uMask, ** vMask;
  int uLatLen, uLonLen;
  int vLatLen, vLonLen;
  int i, j;

  // POINT TO DESIRED ARRAYS AND SET CONSTANTS ---------------------------------

  eta = etaField->solution;

  etaInterp = uField->etaInterp;              // array for interpolated eta at u nodes

  uMask = uField->mask;
  vMask = vField->mask;

  uLatLen = uField->fieldLatLen;
  uLonLen = uField->fieldLonLen;

  vLatLen = vField->fieldLatLen;
  vLonLen = vField->fieldLonLen;


  // INTERPOLATE DISPLACEMENT TO EAST VELOCITY NODE LOCATIONS ------------------

  for (i=0; i<uLatLen; i++)
  {
    for (j=0; j<uLonLen; j++)
    {
      switch (uMask[i][j])
      {
        case 1:
          if (j != uLonLen-1) etaInterp[i][j] = 0.5*(eta[i][j] + eta[i][j+1]);
          else etaInterp[i][uLonLen-1] = 0.5*(eta[i][uLonLen-1] + eta[i][0]);
          break;
        default:
          break;
      }
    }
  }


  // INTERPOLATE EAST VELOCITY TO NORTH VELOCITY NODE LOCATIONS ----------------

  etaInterp = vField->etaInterp;

  for (i=0; i<vLatLen; i++)
  {
    for (j=0; j<vLonLen; j++)
    {
      switch (vMask[i][j])
      {
        case 1:
           etaInterp[i][j] = 0.5*(eta[i][j] + eta[i+1][j]);
          break;
        default:
          break;
      }
    }
  }
};
