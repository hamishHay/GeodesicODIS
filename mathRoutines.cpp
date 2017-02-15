#include "mathRoutines.h"
#include "field.h"
#include "globals.h"
#include "solver.h"

double linearInterp1(Field * field, int i, int j) {
    //first order accurate forward difference method
    if (i < field->fieldLatLen / 2) {
        return field->solution[i + 1][j] + (field->solution[i + 1][j] - field->solution[i + 2][j]);
    }
    else {
        return field->solution[i - 1][j] + (field->solution[i - 1][j] - field->solution[i - 2][j]);
    }
};

double linearInterp2(Field * field, int i, int j) {
    //second order accurate forward difference method
    if (i < field->fieldLatLen / 2) {
        return field->solution[i + 1][j] + (1.5* field->solution[i + 1][j] - 2 * field->solution[i + 2][j] + 0.5* field->solution[i + 3][j]);
    }
    else {
        return field->solution[i - 1][j] + (1.5* field->solution[i - 1][j] - 2 * field->solution[i - 2][j] + 0.5* field->solution[i - 3][j]);

    }
};

double linearInterp4(Field * field, int i, int j) {
    //fourth order accurate forward difference method
    if (i < field->fieldLatLen / 2) {
        return field->solution[i + 1][j] + (25. / 12.* field->solution[i + 1][j] - 4 * field->solution[i + 2][j] + 3 * field->solution[i + 3][j] - 4 / 3. * field->solution[i + 4][j] + .25 * field->solution[i + 5][j]);
    }
    else {
        return field->solution[i - 1][j] + (25. / 12.* field->solution[i - 1][j] - 4 * field->solution[i - 2][j] + 3 * field->solution[i - 3][j] - 4 / 3. * field->solution[i - 4][j] + .25 * field->solution[i - 5][j]);
    }
};

double lagrangeInterp2(Field * field, int i, int j) {
    double L[3];
    double x[4];
    int inc = 1;

    if (i > field->fieldLatLen / 2) inc = -inc;

    x[0] = field->lat[i];
    x[1] = field->lat[i + inc];
    x[2] = field->lat[i + inc*2];
    x[3] = field->lat[i + inc*3];

    //solve for each quadratic term
    L[0] = (x[0] - x[2])*(x[0] - x[3]) / ((x[1] - x[2])*(x[1] - x[3])) * field->solution[i+inc][j]; //L1
    L[1] = (x[0] - x[1])*(x[0] - x[3]) / ((x[2] - x[1])*(x[2] - x[3])) * field->solution[i + inc*2][j]; //L2
    L[2] = (x[0] - x[1])*(x[0] - x[2]) / ((x[3] - x[1])*(x[3] - x[2]))* field->solution[i + inc * 3][j]; //L3

    return L[0]+L[1]+L[2];
};

double lagrangeInterp3(Field * field, int i, int j) {
    double L[4];
    double x[5];
    int inc = 1;

    if (i > field->fieldLatLen / 2) inc = -inc;

    //x[0] = field->lat[i];
    //x[1] = field->lat[i + inc * 2];
    //x[2] = field->lat[i + inc];
    //x[3] = x[0] + field->dLat*radConv;  //field->lat[i + inc];
    //x[4] = x[0] + 2 * field->dLat*radConv;//field->lat[i + inc * 2];

    ////solve for each quadratic term
    //L[0] = (x[0] - x[2])*(x[0] - x[3])*(x[0] - x[4]) / ((x[1] - x[2])*(x[1] - x[3])*(x[1] - x[4])) * field->solution[i + inc * 2][j]; //L1
    //L[1] = (x[0] - x[1])*(x[0] - x[3])*(x[0] - x[4]) / ((x[2] - x[1])*(x[2] - x[3])*(x[2] - x[4])) * field->solution[i + inc][j]; //L2
    //L[2] = (x[0] - x[1])*(x[0] - x[2])*(x[0] - x[4]) / ((x[3] - x[1])*(x[3] - x[2])*(x[3] - x[4])) * field->solution[field->opp[i + inc]][j]; //L3
    //L[3] = (x[0] - x[1])*(x[0] - x[2])*(x[0] - x[3]) / ((x[4] - x[1])*(x[4] - x[2])*(x[4] - x[3])) * field->solution[field->opp[i + inc] + 1][j];

    //return L[0] + L[1] + L[2] + L[3];

    x[0] = field->lat[i];
    x[1] = field->lat[i + inc];
    x[2] = field->lat[i + inc*2];
    x[3] = field->lat[i + inc*3];
    x[4] = field->lat[i + inc*4];

    //solve for each quadratic term
    L[0] = (x[0] - x[2])*(x[0] - x[3])*(x[0] - x[4]) / ((x[1] - x[2])*(x[1] - x[3])*(x[1] - x[4])) * field->solution[i + inc][j]; //L1
    L[1] = (x[0] - x[1])*(x[0] - x[3])*(x[0] - x[4]) / ((x[2] - x[1])*(x[2] - x[3])*(x[2] - x[4])) * field->solution[i + inc*2][j]; //L2
    L[2] = (x[0] - x[1])*(x[0] - x[2])*(x[0] - x[4]) / ((x[3] - x[1])*(x[3] - x[2])*(x[3] - x[4])) * field->solution[i + inc*3][j]; //L3
    L[3] = (x[0] - x[1])*(x[0] - x[2])*(x[0] - x[3]) / ((x[4] - x[1])*(x[4] - x[2])*(x[4] - x[3])) * field->solution[i + inc*4][j];

    return L[0] + L[1] + L[2] + L[3];
};

/**
 *  Derivative of the associated Legendre poynomial at latidude lat
 *
 *  inputs: int l:          degree of the polynomial
 *          int m:          order of the polynomial
 *          double lat:     latitude position (in radians)
 *
 *  returns: double val:    dP_lm @ lat
 *
 */
double assLegendreDLat(int l, int m, double lat)
{
    // double dPlm[5][5];
    double val;
    double colat;

    colat = pi*0.5 - lat;

    // Make sure l doesn't exceed number of polynomials coded
    if (l > 4) {
        printf("WARNING: selected l_max too high. Maximum is 5.\n\n");
    }


    // Nested switches select appropriate function for d(Plm)/dlat
    switch (l) {

    // DEGREE l=0
    case 0:
        val = 0.0;
        break;

    // DEGREE l=1
    case 1:
        switch (m) {
        case 0:
            val = -sin(colat);
            break;
        case 1:
            val = -cos(colat);
            break;
        }
        break;

    // DEGREE l=2
    case 2:
        switch (m) {
        case 0:
            val = -3.*sin(colat)*cos(colat);
            break;
        case 1:
            val = 3.*(pow(sin(colat),2.0) - pow(cos(colat),2.0));
            break;
        case 2:
            val = -6.*sin(colat)*cos(colat);
            break;
        }
        break;

    // DEGREE l=3
    case 3:
        switch (m) {
        case 0:
            val = -0.375*(sin(colat) + 5.*sin(3.*colat));
            break;
        case 1:
            val = 3.*cos(colat)*(5.*pow(sin(colat),2.0) - 0.5*(5.*pow(cos(colat),2.0) - 1.0));
            break;
        case 2:
            val = 15.*sin(colat)*(3.*pow(cos(colat),2.0) - 1.0);
            break;
        case 3:
            val = -45.*pow(sin(colat),2.0) * cos(colat);
            break;
        }
        break;

    // DEGREE l=4
    case 4:
        switch (m) {
        case 0:
            val = -0.3125*(2.0*sin(2.0*colat) + 7.*sin(4.*colat));
            break;
        case 1:
            val = -2.5*(pow(sin(colat),2.0)*(3.-21.*pow(cos(colat),2.0)) + pow(cos(colat),2.0)*(7.*pow(cos(colat),2.0) -3.));
            break;
        case 2:
            val = 30.*sin(colat)*sin(colat)*(7.*pow(cos(colat),2.0) - 4.);
            break;
        case 3:
            val = 105.*pow(sin(colat),2.0)*(1. - 4.*pow(cos(colat),2.0));
            break;
        case 4:
            val = 420.*pow(sin(colat),3.0)*cos(colat);
            break;
        }
        break;
    }

    return val;

};
