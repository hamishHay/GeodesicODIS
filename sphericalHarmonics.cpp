#include "mesh.h"
#include "globals.h"
#include "array3d.h"
#include "array2d.h"
#include "array1d.h"

extern "C"
{
    void extractshcoeff_(double *, double*, double *, int *, int *, double *);
    // void legendrederiv_(double * P, double * dP, int * lmax, double * cosColat);
    // void legendre_(double * P, int * lmax, double * cosColat);
}

void getSHCoeffs(Mesh * grid, Array1D<double> & data, Array2D<double> & sh_coeffs)
{
    
}
