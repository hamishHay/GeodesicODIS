#include "mathRoutines.h"
#include "field.h"
#include "globals.h"
#include "solver.h"

double arcLength(double lat1, double lat2, double lon1, double lon2) {
    // return fabs(acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(fabs(lon1-lon2))));
    return fabs(acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*(cos(lon1)*cos(lon2) + sin(lon1)*sin(lon2))));

}

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
