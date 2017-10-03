#include "mesh.h"
#include "sphericalHarmonics.h"
#include "globals.h"
#include "array3d.h"
#include "array2d.h"
#include "array1d.h"

extern "C"
{
    void extractshcoeffgg_(double *, double*, double *, int *, int *, double *);
    void extractshcoeffll_(double *, int *, int *, double *);
    void legendrederiv_(double * P, double * dP, int * lmax, double * cosColat);
    void legendre_(double * P, int * lmax, double * cosColat);
}

void getSHCoeffsGG(Array2D<double> & coords, Array1D<double> & data, Array3D<double> & sh_coeffs, int node_num, int l_max)
{
    int i,l,m,index;

    double * data_fort;
    double * sh_coeffs_fort;
    double * lats_fort;
    double * lons_fort;

    double rad2deg;

    rad2deg = 1.0/radConv;

    data_fort = new double[node_num];
    sh_coeffs_fort = new double[ 2*(l_max + 1)*(l_max + 1)];

    lats_fort = new double[node_num];
    lons_fort = new double[node_num];

    for (i=0; i<node_num; i++)
    {
        lats_fort[i] = coords(i,0)*rad2deg;
        lons_fort[i] = coords(i,1)*rad2deg;
        data_fort[i] = data(i);
    }

    extractshcoeffgg_(data_fort, lats_fort, lons_fort, &node_num, &l_max, sh_coeffs_fort);

    index = 0;

    for (m=0; m<l_max+1; m++)
    {
        for (l=0; l<l_max+1; l++)
        {
                if (m<=l)
                {
                    sh_coeffs(l, m, 0) = sh_coeffs_fort[index];
                    sh_coeffs(l, m, 1) = sh_coeffs_fort[index + 1];
                }

                index += 2;
            }
        }

    delete[] lats_fort;
    delete[] lons_fort;
    delete[] data_fort;
    delete[] sh_coeffs_fort;

};

void getSHCoeffsLL(Array2D<double> & data,
                   Array3D<double> & sh_coeffs,
                   int N_ll,
                   int l_max)
{
    int i,j,l,m,index,count,lat_len;

    double * data_fort;
    double * sh_coeffs_fort;

    lat_len = 180/N_ll;

    data_fort = new double[360/N_ll * 180/N_ll];
    sh_coeffs_fort = new double[ 2*(l_max + 1)*(l_max + 1)];

    // THIS IS AN EVIL AND SLOW LOOP HELP
    count = 0;
    for (j = 0; j < 360/N_ll; j++) {
        for (i = 0; i < 180/N_ll; i++) {
            data_fort[count] = data(i,j);
            count++;
        }
    }

    extractshcoeffll_(data_fort, &lat_len, &l_max, sh_coeffs_fort);

    index = 0;

    for (m=0; m<l_max+1; m++)
    {
        for (l=0; l<l_max+1; l++)
        {
                if (m<=l)
                {
                    // if (fabs(sh_coeffs_fort[index]) > 1e-16) sh_coeffs(l, m, 0) = sh_coeffs_fort[index];
                    // else sh_coeffs(l, m, 0) = 0.0;
                    // if (fabs(sh_coeffs_fort[index+1]) > 1e-16) sh_coeffs(l, m, 1) = sh_coeffs_fort[index+1];
                    // else sh_coeffs(l, m, 1) = 0.0;

                    sh_coeffs(l, m, 0) = sh_coeffs_fort[index];
                    sh_coeffs(l, m, 1) = sh_coeffs_fort[index+1];
                }

                index += 2;
            }
        }

    delete[] data_fort;
    delete[] sh_coeffs_fort;

};

void getLegendreFuncs(double cosCoLat, Array2D<double> & legendreFuncs, int l_max)
{
    int l, m, index;
    double * legendre_fort;

    legendre_fort = new double[(l_max + 1)*(l_max + 2)/2];

    legendre_(legendre_fort, &l_max, &cosCoLat);

    for (l=0; l<l_max+1; l++)
    {
        for (m=0; m<=l; m++)
        {
            index = l*(l+1)/2 + m;
            legendreFuncs(l, m) = legendre_fort[index];

        }
    }

    delete[] legendre_fort;
};

void getLegendreFuncsDeriv(double cosCoLat, Array2D<double> & legendreFuncsDeriv, int l_max)
{
    int l, m, index;
    double * legendre_fort, * dlegendre_fort;

    legendre_fort = new double[(l_max + 1)*(l_max + 2)/2];
    dlegendre_fort = new double[(l_max + 1)*(l_max + 2)/2];

    if (fabs(cosCoLat) < 0.9999999) legendrederiv_(legendre_fort, dlegendre_fort, &l_max, &cosCoLat);

    for (l=0; l<l_max+1; l++)
    {
        for (m=0; m<=l; m++)
        {
            index = l*(l+1)/2 + m;
            if (fabs(cosCoLat) < 0.9999999) legendreFuncsDeriv(l, m) = dlegendre_fort[index];
            else legendreFuncsDeriv(l, m) = 0.0;
        }
    }

    delete[] legendre_fort;
    delete[] dlegendre_fort;
};
