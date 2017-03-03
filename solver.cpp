#ifdef _WIN32
#include <Windows.h>
#define mkdir CreateDirectory

#elif _WIN64
#include <Windows.h>

#elif __linux__
#include <unistd.h>

#else
#error "OS not supported!"
#endif

#include "solver.h"
#include "globals.h"
#include "field.h"
#include "depth.h"
#include "mesh.h"
#include "mathRoutines.h"
#include "tidalPotentials.h"
#include "energy.h"
#include "outFiles.h"
#include "vector"

#include <math.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cstdio>
#include <sys/stat.h>
#include <errno.h>
#include "H5Cpp.h"
#include <signal.h>
#include <typeinfo>

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

extern "C"
{
    void extractshcoeff_(double *, int *, int *, double *);
    void legendrederiv_(double * P, double * dP, int * lmax, double * cosColat);
    void legendre_(double * P, int * lmax, double * cosColat);
}

int flag = 0;

void CatchExit(int sig) {
    printf("%s\n", "Caught Terminate Signal...");
    flag = 1;
}

Solver::Solver(int type, int dump, Globals * Consts, Mesh * Grid, Field * UGradLon, Field * UGradLat, Field * VelU, Field * VelV, Field * Eta, Energy * EnergyField, Depth * Depth_h) {
    solverType = type;
    dumpTime = dump;
    Out = &Consts->Output;

    //Assign member pointers to passed in pointers from main.cpp
    consts = Consts;
    grid = Grid;

    u = VelU;
    v = VelV;
    eta = Eta;
    energy = EnergyField;

    loading = false;

    radius = consts->radius.Value();
    angVel = consts->angVel.Value();
    theta = consts->theta.Value();
    ecc = consts->e.Value();
    smAxis = consts->a.Value();

    etaOldArray = eta->solution;
    etaNewArray = eta->MakeSolutionArrayCopy();
    oceanLoadingArrayEta = eta->MakeSolutionArrayCopy();

    uOldArray = u->solution;
    uNewArray = u->MakeSolutionArrayCopy();
    uDissArray = u->MakeSolutionArrayCopy();
    vNEAvgArray = u->MakeSolutionArrayCopy();
    oceanLoadingArrayU = u->MakeSolutionArrayCopy();

    vOldArray = v->solution;
    vNewArray = v->MakeSolutionArrayCopy();
    vDissArray = v->MakeSolutionArrayCopy();
    uSWAvgArray = v->MakeSolutionArrayCopy();
    oceanLoadingArrayV = v->MakeSolutionArrayCopy();

    dUlat = UGradLat;
    dUlatArray = dUlat->solution;

    dUlon = UGradLon;
    dUlonArray = dUlon->solution;

    depth = Depth_h;
    depthArray = Depth_h->solution;

    l_max = consts->l_max.Value();

    SH_cos_coeff = new double*[l_max+1];
    SH_sin_coeff = new double*[l_max+1];
    for (int i=0; i<l_max+1; i++) {
        SH_cos_coeff[i] = new double[l_max+1];
        SH_sin_coeff[i] = new double[l_max+1];
    }

    // SH_cos_coeff = new double[(l_max+1)*(l_max+1)];
    // SH_sin_coeff = new double[(l_max+1)*(l_max+1)];

    uLatLen = u->fieldLatLen;
    uLonLen = u->fieldLonLen;
    udLon = u->dLon;
    udLat = u->dLat;

    vLatLen = v->fieldLatLen;
    vLonLen = v->fieldLonLen;
    vdLon = v->dLon;
    vdLat = v->dLat;

    etaLatLen = eta->fieldLatLen;
    etaLonLen = eta->fieldLonLen;
    etadLon = eta->dLon;
    etadLat = eta->dLat;

    cellArea = grid->cellArea;

    dt = consts->timeStep.Value();

    // etaLegendreArray = new double**[etaLatLen];
    // for (int i=0; i<etaLatLen; i++) {
    //     etaLegendreArray[i] = new double*[l_max+1];
    //     for (int l=0; l<l_max+1; l++) {
    //         etaLegendreArray[i][l] = new double[l_max+1];
    //         for (int m=0; m <= l; m++) {
    //             etaLegendreArray[i][l][m] = 0.0; //assLegendre(l, m, eta->cosCoLat[i]);
    //         }
    //     }
    // }

    //matrix size is x*y*z
  // double ***etaLegendreArray;
  // etaLegendreArray = (double***) malloc(etaLatLen*sizeof(double**));
  // etaLegendreArray[0] = (double**)malloc(etaLatLen*(l_max+1)*sizeof(double*));
  // etaLegendreArray[0][0] = (double*)malloc(etaLatLen*(l_max+1)*(l_max+1)*sizeof(double));
  //
  // for(int j=1; j<(l_max+1); j++)
  //   etaLegendreArray[0][j]=etaLegendreArray[0][j-1]+(l_max+1);
  //  for(int i=1; i<etaLatLen; i++){
  //   etaLegendreArray[i]=etaLegendreArray[i-1]+(l_max+1);
  //   etaLegendreArray[i][0]=etaLegendreArray[i-1][0]+(l_max+1)*(l_max+1);
  //   for(int j=1; j<(l_max+1); j++){
  //     etaLegendreArray[i][j]=etaLegendreArray[i][j-1]+etaLatLen;
  //   }
  //  }

    etaLegendreArray = new double[etaLatLen*(l_max+1)*(l_max+1)];

    // etaLegendreArray = new double**[etaLatLen];
    // // etaLegendreArray[0] = new double*[etaLatLen*(l_max+1)]
    // for (int i=1; i<etaLatLen; i++) {
    //     *(etaLegendreArray[i]) = new double*[l_max+1];
    //     for (int l=0; l<l_max+1; l++) {
    //         etaLegendreArray[i][l] = new double[l_max+1];
    //         for (int m=0; m <= l; m++) {
    //             etaLegendreArray[i][l][m] = 0.0; //assLegendre(l, m, eta->cosCoLat[i]);
    //         }
    //     }
    // }

    // etaLegendreArray = new double[etaLatLen*(l_max+1)*(l_max+1)];


    uLegendreArray = new double**[uLatLen];
    for (int i=0; i<uLatLen; i++) {
        uLegendreArray[i] = new double*[l_max+1];
        for (int l=0; l<l_max+1; l++) {
            uLegendreArray[i][l] = new double[l_max+1];
            for (int m=0; m <= l; m++) {
                uLegendreArray[i][l][m] = 0.0; //assLegendre(l, m, u->cosCoLat[i]);
            }
        }
    }

    Legendre();

    vdLegendreArray = new double**[vLatLen];
    for (int i=0; i<vLatLen; i++) {
        vdLegendreArray[i] = new double*[l_max+1];
        for (int l=0; l<l_max+1; l++) {
            vdLegendreArray[i][l] = new double[l_max+1];
            for (int m=0; m <= l; m++) {
                vdLegendreArray[i][l][m] = 0.0;
            }
        }
    }

    LegendreDeriv();

    // etaCosMLon = new double*[etaLonLen];
    // etaSinMLon = new double*[etaLonLen];
    // for (int j=0; j<etaLonLen; j++) {
    //     etaCosMLon[j] = new double[l_max+1];
    //     etaSinMLon[j] = new double[l_max+1];
    //     for (int m=0; m<l_max+1; m++) {
    //         etaCosMLon[j][m] = eta->cosMLon[m][j];
    //         etaSinMLon[j][m] = eta->sinMLon[m][j];
    //         // std::cout<<&etaSinMLon[j][m]<<std::endl;
    //     }
    // }
    // etaSinMLon = new double*[etaLonLen];
    // etaSinMLon[0] = new double[etaLonLen*(l_max+1)];
    // for (int j=1; j<etaLonLen; j++) {
    //     etaSinMLon[j] = &etaSinMLon[0][j*(l_max+1)];
    //     for (int m=0; m<l_max+1; m++) {
    //         etaSinMLon[j][m] = eta->sinMLon[m][j];
    //         // std::cout<<&etaSinMLon[j][m]<<std::endl;
    //     }
    // }

    // std::cout<<"BREAK\n\n"<<std::endl;
    // etaCosMLon = new double*[etaLonLen];
    // etaCosMLon[0] = new double[etaLonLen*(l_max+1)];
    // for (int j=1; j<etaLonLen; j++) {
    //     etaCosMLon[j] = &etaCosMLon[0][j*(l_max+1)];
    //     for (int m=0; m<l_max+1; m++) {
    //         etaCosMLon[j][m] = eta->cosMLon[m][j];
    //         // std::cout<<&etaCosMLon[j][m]<<std::endl;
    //     }
    // }

    etaSinMLon = new double[etaLonLen*(l_max+1)];
    for (int j=0; j<etaLonLen; j++) {
        for (int m=0; m<l_max+1; m++) {
            etaSinMLon[j*(l_max+1) + m] = eta->sinMLon[m][j];
        }
    }


    etaCosMLon = new double[etaLonLen*(l_max+1)];
    for (int j=0; j<etaLonLen; j++) {
        for (int m=0; m<l_max+1; m++) {
            etaCosMLon[j*(l_max+1) + m] = eta->cosMLon[m][j];
        }
    }


    vCosMLon = new double*[vLonLen];
    vSinMLon = new double*[vLonLen];
    for (int j=0; j<vLonLen; j++) {
        vCosMLon[j] = new double[l_max+1];
        vSinMLon[j] = new double[l_max+1];
        for (int m=0; m<l_max+1; m++) {
            vCosMLon[j][m] = v->cosMLon[m][j];
            vSinMLon[j][m] = v->sinMLon[m][j];
        }
    }

    uCosMLon = new double*[uLonLen];
    uSinMLon = new double*[uLonLen];
    for (int j=0; j<uLonLen; j++) {
        uCosMLon[j] = new double[l_max+1];
        uSinMLon[j] = new double[l_max+1];
        for (int m=0; m<l_max+1; m++) {
            uCosMLon[j][m] = u->cosMLon[m][j];
            uSinMLon[j][m] = u->sinMLon[m][j];
        }
    }

    loadK = new double[l_max+1];
    loadH = new double[l_max+1];
    gammaFactor = new double[l_max+1];

    double mu_bar = 0.0;

    for (int l=0; l<l_max+1; l++)
    {
        mu_bar = 6.78e9/(1610.0*consts->g.Value()*consts->radius.Value());
        mu_bar *= (double)(2*l*l + 4*l +3)/((double)l);

        loadK[l] = -1.0 / (1.0 + mu_bar);
        loadH[l] = -1.0 / (1.0 + mu_bar) * (2.0*(double)l + 1.0)/3.0;

        gammaFactor[l] = 1.0 - (1.0 + loadK[l] - loadH[l]) * 3000.0 / ((2*l + 1) * 1610.0);
    }

    if (consts->potential.Value() == "ECC_RAD")         tide = ECC_RAD;
    else if (consts->potential.Value() == "ECC_LIB")    tide = ECC_LIB;
    else if (consts->potential.Value() == "ECC")        tide = ECC;
    else if (consts->potential.Value() == "OBLIQ")      tide = OBLIQ;
    else if (consts->potential.Value() == "FULL")       tide = FULL;
    else if (consts->potential.Value() == "TOTAL")      tide = TOTAL;
    else if (consts->potential.Value() == "ECC_W3")     tide = ECC_W3;
    else if (consts->potential.Value() == "OBLIQ_W3")   tide = OBLIQ_W3;
    else {
        outstring << "No potential forcing found." << std::endl;
        Out->Write(ERR_MESSAGE, &outstring);
        Out->TerminateODIS();
    }

    CreateHDF5FrameWork();

    signal(SIGINT, CatchExit);

};

int Solver::Legendre(void) {
    double * fort_dplm_dtheta;
    double * fort_plm;
    int n,i,j,l,m,index;
    double val;
    double val2;

    n = (l_max + 1)*(l_max + 2)/2;
    // Create 1d arrays for for fortran to return values with
    fort_plm = new double[n];
    fort_dplm_dtheta = new double[n];


    for (i=0; i<etaLatLen; i++) {
        val = eta->cosCoLat[i];
        legendre_(fort_plm, &l_max, &val);

        for (l=0; l<l_max+1; l++) {
            for (m=0; m<=l; m++) {
                index = l*(l+1)/2 + m;
                val2 = fort_plm[index];
                etaLegendreArray[(i*(l_max+1) + l)*(l_max+1) + m] = val2;
                // etaLegendreArray[i][l][m] = val2; //FOR SOME REASON ARRAY MUST BE SET WITH VAL2??????
                // etaLegendreArray[i*(l_max+1)*(l_max+1) + l*(l_max+1) + m] = val2;
                uLegendreArray[i][l][m] = val2;

                // if (l==2 && m==0) std::cout<< std::setprecision(8)<<i<<'\t'<<l<<'\t'<<m<<'\t'<<uLegendreArray[i][l][m]<<'\t'<<std::endl;
                // std::cout<<i<<'\t'<<l<<'\t'<<m<<'\t'<<vdLegendreArray[i][l][m]<<'\t'<<val<<std::endl;
            }
        }
    }
    // std::cout<<std::endl<<std::endl;

    delete[] fort_plm;
    delete[] fort_dplm_dtheta;

};

int Solver::LegendreDeriv(void) {
    double * fort_dplm_dtheta;
    double * fort_plm;
    int n,i,j,l,m,index;
    double val;
    double val2;

    n = (l_max + 1)*(l_max + 2)/2;
    // Create 1d arrays for for fortran to return values with
    fort_plm = new double[n];
    fort_dplm_dtheta = new double[n];


    for (i=0; i<vLatLen; i++) {
        val = v->cosCoLat[i];
        legendrederiv_(fort_plm, fort_dplm_dtheta, &l_max, &val);

        for (l=0; l<l_max+1; l++) {
            for (m=0; m<=l; m++) {
                index = l*(l+1)/2 + m;
                val2 = fort_dplm_dtheta[index];

                vdLegendreArray[i][l][m] = val2*sin(pi*0.5 - v->lat[i]); //FOR SOME REASON ARRAY MUST BE SET WITH VAL2??????


                // if (l==2 && m==0) std::cout<<i<<'\t'<<l<<'\t'<<m<<'\t'<<vdLegendreArray[i][l][m]<<'\t'<<std::endl;
            }
        }

    }

    // Out->TerminateODIS();

    delete[] fort_plm;
    delete[] fort_dplm_dtheta;

};

int Solver::InitialConditions(void) {
    bool action = consts->init.Value();

    outstring << "Use initial conditions: ";

    if (action) {
        outstring << "Yes." << std::endl;
        Out->Write(OUT_MESSAGE, &outstring);
        ReadInitialConditions(action);
        return 1;
    }
    else {
        ReadInitialConditions(action);
        outstring << "No." << std::endl;
        Out->Write(OUT_MESSAGE, &outstring);
    }
    return 0;
};

void Solver::Solve() {
    outstring << "Entering solver: ";
    switch (solverType) {
    case 0:
        outstring << "Explicit\n";
        Explicit();
        break;
    default:
        outstring << "No solver type selected.\n Terminating execution.";
    }
};

void Solver::UpdatePotential() {
    switch (tide) {
    case ECC_RAD:
        deg2EccRad(dUlat, dUlon, simulationTime, radius, angVel, ecc);
        break;

    case ECC_LIB:
        deg2EccLib(dUlat, dUlon, simulationTime, radius, angVel, ecc);
        break;

    case ECC:
        deg2Ecc(dUlat, dUlon, simulationTime, radius, angVel, ecc);
        break;

    case OBLIQ:
        deg2Obliq(dUlat, dUlon, simulationTime, radius, angVel, theta);
        break;

    // case FULL:
    //     UpdateFullPotential();
    //     break;
    //
    // case TOTAL:
    //     UpdateTotalPotential();
    //     break;

    case ECC_W3:
        deg3Ecc(dUlat, dUlon, simulationTime, radius, smAxis, angVel, ecc);
        break;

    case OBLIQ_W3:
        deg3Obliq(dUlat, dUlon, simulationTime, radius, smAxis, angVel, theta);
        break;
    }
};


void Solver::UpdateEastVel(){
    double coriolis = 0.0;//0;
    double tidalForce = 0;

    double dSurfLon = 0;
    double surfHeight = 0;
    double eastEta = 0;
    double westEta = 0;
    double oceanLoadingTerm = 0;

    double alpha = consts->alpha.Value();
    double angVel = consts->angVel.Value();

    double a,b,c,d,sum;

    double old;

    int i, j, i_a, j_a;

    for (i = 1; i < uLatLen - 1; i++) {
        for (j = 0; j < uLonLen; j++) {
            i_a = i*2;
            j_a = j*2;

            if (j < vLonLen - 1) {
                vNEAvgArray[i][j] = 0.25*vOldArray[i][j]
                                    + 0.25*vOldArray[i - 1][j]
                                    + 0.25*vOldArray[i][j + 1]
                                    + 0.25*vOldArray[i - 1][j + 1];
            }
            else {
                vNEAvgArray[i][j] = 0.25*vOldArray[i][j]
                                    + 0.25*vOldArray[i - 1][j]
                                    + 0.25*vOldArray[i][0]
                                    + 0.25*vOldArray[i - 1][0];
            }
        }
    }


    switch (consts->fric_type) {
    case LINEAR:
        for (int i = 1; i < uLatLen - 1; i++) {
            for (int j = 0; j < uLonLen; j++) {
                uDissArray[i][j] = uOldArray[i][j] * alpha;
            }
        }
        break;

    case QUADRATIC:
        double alphah = 0.0;
        int i_h = 0;
        int j_h = 0;

        for (int i = 1; i < uLatLen - 1; i++) {
            i_h = i*2;
            for (int j = 0; j < uLonLen; j++) {
                j_h = j*2;
                alphah = alpha / (depthArray[i_h][j_h+1]);
                uDissArray[i][j] = alphah * uOldArray[i][j] * sqrt(uOldArray[i][j]*uOldArray[i][j] + vNEAvgArray[i][j]*vNEAvgArray[i][j]);
            }
        }
        break;
    }

    double loveRadius = consts->loveReduct.Value() / consts->radius.Value();
    double gRadius = consts->g.Value() / consts->radius.Value();
    double coriolisFactor = 0;
    double tidalFactor = 0;
    double surfFactor = 0;

    double * uCosLat = u->cosLat;
    double * uSinLat = u->sinLat;

    double g = consts->g.Value();
    double r = consts->radius.Value();

    int i_h = 0;
    int j_h = 0;

    for (int i = 1; i < uLatLen - 1; i++) {
        i_h = i*2;

        coriolisFactor =  2. * angVel*uSinLat[i];
        tidalFactor = loveRadius/uCosLat[i];
        surfFactor = gRadius/uCosLat[i];
        for (int j = 0; j < uLonLen; j++) {
            j_h = j*2;
            if (j != uLonLen - 1) eastEta = depthArray[i_h][j_h + 2] + etaOldArray[i][j+1];
            else eastEta = depthArray[i_h][0] + etaOldArray[i][0];

            westEta = depthArray[i_h][j_h] + etaOldArray[i][j];

            dSurfLon = (eastEta - westEta) / (etadLon);

            coriolis = coriolisFactor * vNEAvgArray[i][j];
            tidalForce = tidalFactor * dUlonArray[i][j];


            if (!loading) {
                surfHeight = surfFactor*dSurfLon;
                oceanLoadingTerm = 0.0;
            }
            else {
                oceanLoadingTerm = surfFactor*oceanLoadingArrayU[i][j];
                surfHeight = 0.0;
            }

            uNewArray[i][j] = (coriolis - surfHeight - oceanLoadingTerm + tidalForce - uDissArray[i][j])*dt + uOldArray[i][j];
        }
    }


    double du_east, du_west, du_centre, du_avg, u_avg;
    for (int j = 0; j < uLonLen; j++) {
        uNewArray[0][j] = linearInterp1Array(u,uNewArray, 0, j);
        uNewArray[uLatLen - 1][j] = linearInterp1Array(u,uNewArray, uLatLen - 1, j);
        // uNewArray[0][j] = lagrangeInterp4ArrayCenter(u,uNewArray, 0, j);
        // uNewArray[uLatLen - 1][j] = lagrangeInterp4ArrayCenter(u,uNewArray, uLatLen - 1, j);
        // uNewArray[0][j] = 0.0;//uNewArray[1][j];
        // uNewArray[uLatLen - 1][j] = 0.0;//uNewArray[uLatLen - 2][j];

    }

}

int Solver::UpdateNorthVel(){
    double coriolis = 0;
    double tidalForce = 0;
    double dSurfLat = 0;
    double surfHeight = 0;

    double northEta = 0;
    double southEta = 0;

    double oceanLoadingTerm = 0;
    double loveRadius = consts->loveReduct.Value() / consts->radius.Value();
    double gRadius = consts->g.Value() / consts->radius.Value();
    double coriolisFactor = 0;

    double g = consts->g.Value();
    double r = consts->radius.Value();

    double alpha = consts->alpha.Value();
    double angVel = consts->angVel.Value();

    double a,b,c,d,sum;
    int i_a;
    int j_a;
    double leftU;
    double rightU;

    for (int i = 0; i < vLatLen; i++) {
        i_a = i*2;
        for (int j = 0; j < vLonLen; j++) {
            j_a = j*2;
            if (j > 0) {
                uSWAvgArray[i][j] = 0.25*(uOldArray[i][j]
                                    + uOldArray[i + 1][j]
                                    + uOldArray[i][j - 1]
                                    + uOldArray[i + 1][j - 1]);
            }
            else {
                uSWAvgArray[i][j] = 0.25*(uOldArray[i][j]
                                    + uOldArray[i + 1][j]
                                    + uOldArray[i][uLonLen - 1]
                                    + uOldArray[i + 1][uLonLen - 1]);
            }
        }



    }

    switch (consts->fric_type) {
    case LINEAR:
        for (int i = 0; i < vLatLen; i++) {
            for (int j = 0; j < vLonLen; j++) {
                vDissArray[i][j] = vOldArray[i][j] * alpha;
            }
        }
        break;

    case QUADRATIC:
        double alphah = 0.0;
        int i_h = 0;
        int j_h = 0;
        for (int i = 0; i < vLatLen; i++) {
            i_h = i*2;
            for (int j = 0; j < vLonLen; j++) {
                j_h = j*2;
                //alphah = alpha / (depthArray[i_h+1][j_h] + etaVAvgArray[i][j]);
                alphah = alpha / (depthArray[i_h+1][j_h]);
                vDissArray[i][j] = alphah * vOldArray[i][j] * sqrt(vOldArray[i][j] * vOldArray[i][j] + uSWAvgArray[i][j]*uSWAvgArray[i][j]);
            }
        }
        break;
    }


    int i_h = 0;
    int j_h = 0;

    for (int i = 0; i < vLatLen; i++) {
        coriolisFactor = 2. * angVel * v->sinLat[i];
        i_h = i*2;

        for (int j = 0; j < vLonLen; j++) {
            j_h = j*2;

            northEta = depthArray[i_h][j_h] + etaOldArray[i][j];
            southEta = depthArray[i_h+2][j_h] + etaOldArray[i+1][j];

            dSurfLat = (northEta - southEta) / etadLat;

            if (!loading) {
                surfHeight = gRadius*dSurfLat;
                oceanLoadingTerm = 0.0;
            }
            else {
                oceanLoadingTerm = gRadius*oceanLoadingArrayV[i][j];
                surfHeight = 0.0;
            }

            coriolis = coriolisFactor*uSWAvgArray[i][j];

            tidalForce = loveRadius * dUlatArray[i][j];

            vNewArray[i][j] = (-coriolis - surfHeight - oceanLoadingTerm + tidalForce - vDissArray[i][j])*dt + vOldArray[i][j];

        }
    }

}

void Solver::UpdateSurfaceHeight(){
    double vGrad = 0;
    double uGrad = 0;
    double northv = 0;
    double southv = 0;
    double eastu = 0;
    double westu = 0;

    double cosLat;
    double vdLat = v->dLat;
    double vdLon = v->dLon;
    double r = consts->radius.Value();
    double hRadius = 0.0;
    double h = consts->h.Value();
    int i_h = 0;
    int j_h = 0;


    for (int i = 1; i < etaLatLen-1; i++) {
        for (int j = 0; j < etaLonLen; j++) {
            cosLat = eta->cosLat[i];
            i_h = i*2;

            j_h = j*2;

            northv = ( depthArray[i_h - 1][j_h]) * vNewArray[i - 1][j] * v->cosLat[i - 1];
            southv = ( depthArray[i_h + 1][j_h]) * vNewArray[i][j] * v->cosLat[i];
            //
            vGrad = (northv - southv) / vdLat;

            //northv = vNewArray[i - 1][j];// * v->cosLat[i - 1];
            //southv = vNewArray[i][j];// * v->cosLat[i];

            //vGrad = h*(eta->cosLat[i]*(northv - southv) / vdLat
            //        - eta->sinLat[i]*(vNewArray[i - 1][j] + vNewArray[i][j])*0.5);


            if (j > 0) {
                eastu = (depthArray[i_h][j_h + 1]) * uNewArray[i][j];
                westu = (depthArray[i_h][j_h - 1]) * uNewArray[i][j - 1];
            }
            else {
                eastu = (depthArray[i_h][j_h + 1]) * uNewArray[i][j];
                westu = (depthArray[i_h][uLonLen*2 - 1]) * uNewArray[i][uLonLen - 1];
            }

            uGrad = (eastu - westu) / vdLon;

            hRadius = 1.0 / r;

            etaNewArray[i][j] = hRadius/cosLat*(-vGrad - uGrad)*dt + etaOldArray[i][j];
        }
    }

    // if (!loading) {
        // double deta_east, deta_west, deta_centre, deta_avg, eta_avg;
        // for (int j = 0; j < etaLonLen; j++) {
        //     // if (loading) {
        //     //     etaNewArray[0][j] =  etaNewArray[1][j];
        //     //     etaNewArray[etaLatLen - 1][j] = etaNewArray[etaLatLen - 2][j];
        //     // }
        //     //  etaNewArray[0][j] =  etaNewArray[1][j];
        //     //  etaNewArray[etaLatLen - 1][j] = etaNewArray[etaLatLen - 2][j];
        //     // etaNewArray[0][j] = lagrangeInterp2Array(eta,etaNewArray, 0, j);
        //     // etaNewArray[etaLatLen - 1][j] = lagrangeInterp2Array(eta,etaNewArray, etaLatLen - 1, j);
        //     // else {
        //     // etaNewArray[0][j] = lagrangeInterp4ArrayCenter(eta,etaNewArray, 0, j);
        //     // etaNewArray[etaLatLen - 1][j] = lagrangeInterp4ArrayCenter(eta,etaNewArray, etaLatLen - 1, j);
        //     etaNewArray[0][j] = linearInterp4Array(eta,etaNewArray, 0, j);
        //     etaNewArray[etaLatLen - 1][j] = linearInterp4Array(eta,etaNewArray, etaLatLen - 1, j);
        //     // }
        //
        //     // eta_avg = 0.5*(etaNewArray[1][j] + etaNewArray[2][j]);
        //     // deta_centre = (etaNewArray[1][j] - etaNewArray[2][j])/etadLat;
        //     // if (j > 0 && j < etaLonLen-1) {
        //     //     deta_east = (etaNewArray[1][j+1] - etaNewArray[2][j+1])/etadLat;
        //     //     deta_west = (etaNewArray[1][j-1] - etaNewArray[2][j-1])/etadLat;
        //     // }
        //     // else if (j == 0) {
        //     //     deta_east = (etaNewArray[1][j+1] - etaNewArray[2][j+1])/etadLat;
        //     //     deta_west = (etaNewArray[1][etaLonLen-1] - etaNewArray[2][etaLonLen-1])/etadLat;
        //     // }
        //     // else if (j == etaLonLen-1) {
        //     //     deta_east = (etaNewArray[1][0] - etaNewArray[2][0])/etadLat;
        //     //     deta_west = (etaNewArray[1][j-1] - etaNewArray[2][j-1])/etadLat;
        //     // }
        //     //
        //     // deta_avg = 0.5*deta_centre + 0.25*(deta_east + deta_west);
        //     //
        //     // etaNewArray[0][j] = eta_avg - etadLat*deta_avg;
        // }
        //
        // double npoleSum = 0;
        // double spoleSum = 0;
        // for (int j = 0; j < etaLonLen; j++) {
        //     npoleSum += etaNewArray[0][j];
        //     spoleSum += etaNewArray[etaLatLen - 1][j];
        // }
        // npoleSum = npoleSum / etaLonLen;
        // spoleSum = spoleSum / etaLonLen;

};

int Solver::FindAverages(void) {

}

int Solver::ExtractSHCoeff(void) {
    double * fort_array;
    double * fort_harm_coeff;

    int coeff_num = 2*(l_max + 1)*(l_max + 1);

    int i_len = eta->ReturnFieldLatLen()-1;
    int j_len = eta->ReturnFieldLonLen();

    // if (i_len%2 != 0) i_len -= 1; //minus 1 required as SHTOOLS requires even samples in latitude

    int n = i_len*j_len;

    int i,j,k,l,m;

    fort_array = new double[n];
    fort_harm_coeff = new double[coeff_num];

    //THIS LOOP IS SLOW
    int count = 0;
    for (j = 0; j<j_len; j++) {
        for (i = 0; i<i_len; i++) {
            fort_array[count] = etaNewArray[i][j];
            count++;
        }
    }

    extractshcoeff_(fort_array, &i_len, &l_max, fort_harm_coeff);

    // TODO - Speed up loop
    count = 0;
    for (j=0; j<l_max+1; j++) {
        for (k=0; k<l_max+1; k++) {
            if (fabs(fort_harm_coeff[count]) < 1e-10) SH_cos_coeff[k][j] = 0.0;
            else SH_cos_coeff[k][j] = fort_harm_coeff[count];
            // SH_cos_coeff[k][j] = fort_harm_coeff[count];

            if (fabs(fort_harm_coeff[count+1]) < 1e-10) SH_sin_coeff[k][j] = 0.0;
            else SH_sin_coeff[k][j] = fort_harm_coeff[count+1];
            // // SH_sin_coeff[k][j] = fort_harm_coeff[count+1];

            // if (fabs(fort_harm_coeff[count]) < 1e-10) SH_cos_coeff[k*(l_max+1) + j] = 0.0;
            // else SH_cos_coeff[k*(l_max+1) + j] = fort_harm_coeff[count];
            // // SH_cos_coeff[k][j] = fort_harm_coeff[count];
            //
            // if (fabs(fort_harm_coeff[count+1]) < 1e-10) SH_sin_coeff[k*(l_max+1) + j] = 0.0;
            // else SH_sin_coeff[k*(l_max+1) + j] = fort_harm_coeff[count+1];

            count+=2;
        }
    }

    delete[] fort_array;
    delete[] fort_harm_coeff;

    return 1;
}

int Solver::UpdateLoading(void) {
    double loading = 0.0;
    double loadingTotal = 0.0;
    double loadingEta = 0.0;
    double loadingEtaTotal = 0.0;
    double loadingFactor = 0.0;
    double omega = consts->angVel.Value();
    double omegaTime = omega*simulationTime;
    double lon;
    double normalise;
    double loadingDLat, loadingDLon;
    double loadingDLatTotal, loadingDLonTotal;
    double loadingEtaNorth, loadingEtaSouth;
    int i,j,l,m,k, degree;

    ExtractSHCoeff();

    // NORTH VELOCITY
    // for (i = 0; i < vLatLen; i++) {
    //     for (j = 0; j < vLonLen; j++) {
    //         loadingDLatTotal = 0.;
    //         oceanLoadingArrayV[i][j] = 0.0;
    //         for (l=1; l<l_max+1; l++) {
    //             loadingDLat = 0.0;
    //             for (m=0; m<=l; m++) {
    //                 loadingDLat += vdLegendreArray[i][l][m] * (SH_cos_coeff[l*(l_max+1) + m]*vCosMLon[j][m] + SH_sin_coeff[l*(l_max+1) + m]*vSinMLon[j][m]);
    //             }
    //             loadingDLat *= gammaFactor[l];
    //             loadingDLatTotal += loadingDLat;
    //         }
    //         oceanLoadingArrayV[i][j] = loadingDLatTotal;
    //     }
    // }
    //
    // // EAST VELOCITY
    // for (i = 0; i < uLatLen; i++) {
    //     for (j = 0; j < uLonLen; j++) {
    //         loadingDLonTotal = 0.;
    //
    //         for (l=1; l<l_max+1; l++) {
    //             loadingDLon = 0.0;
    //             for (m=0; m<=l; m++) {
    //                 loadingDLon += uLegendreArray[i][l][m] * (-SH_cos_coeff[l*(l_max+1) + m]*(double)m*uSinMLon[j][m] + SH_sin_coeff[l*(l_max+1) + m]*(double)m*uCosMLon[j][m]);
    //             }
    //             loadingDLon *= gammaFactor[l];
    //             loadingDLonTotal += loadingDLon;
    //
    //         }
    //         oceanLoadingArrayU[i][j] = loadingDLonTotal;
    //     }
    // }

    return 1;

}

int Solver::InterpPoles() {
    double loadingEta;
    int i,j,l,m,k, degree;

    if (!loading) ExtractSHCoeff();

    for (i=0; i<etaLatLen; i++) {
        for (j = 0; j < etaLonLen; j++) {
            loadingEta = 0.0;
            for (l=0; l<l_max+1; l++) {
                for (m=0; m<=l; m++) {
                    loadingEta += etaLegendreArray[(i*(l_max+1) + l)*(l_max+1) + m]//etaLegendreArray[i][l][m]
                                  * (SH_cos_coeff[l][m]
                                  * etaCosMLon[j*(l_max+1) + m] //etaCosMLon[j][m]
                                  + SH_sin_coeff[l][m]
                                  * etaSinMLon[j*(l_max+1) + m]);//etaSinMLon[j][m]);
                }
            }
            etaNewArray[i][j] = loadingEta;
        }
    }

    return 1;
}


int Solver::Explicit() {

    int outCount = 0;
    double timeStepCount = 0;
    int inc = (int) (consts->period.Value()/dt);

    InitialConditions(); // SHOULD BE IN CONSTRUCTOR

    outstring << "Entering time loop:\n\n";
    outstring << "\t\t End time: \t" << consts->endTime.Value() / 86400.0 << " days\n";
    outstring << "\t\t Time step: \t" << dt << " seconds\n\n";
    Out->Write(OUT_MESSAGE, &outstring);

    DumpSolutions(-1,timeStepCount);


    energy->mass->UpdateMass();

    //Update cell energies and globally averaged energy
    energy->UpdateKinE(uNewArray,vNewArray);

    while (simulationTime <= consts->endTime.Value() && !energy->converged) {

        if (flag) {
            return -1;
        }

        timeStepCount+=dt;
        outputCount+=dt;

        UpdatePotential();

        simulationTime += dt;

        UpdateNorthVel();
        UpdateEastVel();

        UpdateSurfaceHeight();

        if (!loading) {
            if (simulationTime > 1.1*consts->endTime.Value()) {
                printf("Kicking in ocean loading at orbit num %d\n", output);
                loading = true;
            }
        }

        if (loading) {
            UpdateLoading();
            // output++;
            // DumpFields(output);
            // flag = true;
        }
        // else InterpPoles();
        InterpPoles();

        for (int i = 0; i < vLatLen; i++) {
            for (int j = 0; j < vLonLen; j++) {
                vOldArray[i][j] = vNewArray[i][j];
            }
        }

        for (int i = 0; i < uLatLen; i++) {
            for (int j = 0; j < uLonLen; j++) {
                uOldArray[i][j] = uNewArray[i][j];
            }
        }

        for (int i = 0; i < etaLatLen; i++) {
            for (int j = 0; j < etaLonLen; j++) {
                etaOldArray[i][j] = etaNewArray[i][j];
            }
        }

        energy->UpdateKinE(uNewArray,vNewArray);

        energy->UpdateDtKinEAvg();

        //Check for output
        if (timeStepCount >= consts->period.Value()) {
            orbitNumber++;
            timeStepCount -= consts->period.Value();

            printf("TIME: %.2f, number: %d, dissipation: %.3f, h_0: %.3f\n", simulationTime/consts->period.Value(), output, energy->currentDissEAvg, consts->h.Value());

            energy->UpdateOrbitalKinEAvg(inc);

            // Check for convergence

            if (energy->converged) convergeCount++;

            outstring << std::fixed << std::setprecision(2) << simulationTime / 86400.0 << " days: \t" << 100 * (simulationTime / consts->endTime.Value()) << "%\t" << output;
            Out->Write(OUT_MESSAGE, &outstring);

            output++;
            outCount++;


            outCount = 1;

            energy->timePos = 0; //Reset time position after energy data output




            DumpFields(output);
        }
        else if (timeStepCount >= consts->period.Value()*consts->outputTime.Value()*outCount) {
            output++;

            outCount++;


            DumpFields(output);

        }
        iteration++;

    }

    outstring << "\nSimulation complete.\n\nEnd time: \t\t" << simulationTime / 86400.0 << "\n";
    outstring << "Total iterations: \t" << iteration;
    Out->Write(OUT_MESSAGE, &outstring);

    return 1;

};

void Solver::DumpSolutions(int out_num, double time) {

    if (out_num == -1) {

        //Only for windows

        //remove(&(consts->path + SEP + "diss.txt")[0]);

        std::ofstream uLat(consts->path + SEP + "Grid" + SEP + "u_lat.txt", std::ofstream::out);
        std::ofstream uLon(consts->path + SEP + "Grid" + SEP + "u_lon.txt", std::ofstream::out);
        std::ofstream vLat(consts->path + SEP + "Grid" + SEP + "v_lat.txt", std::ofstream::out);
        std::ofstream vLon(consts->path + SEP + "Grid" + SEP + "v_lon.txt", std::ofstream::out);

        for (int j = 0; j < u->ReturnFieldLonLen(); j++) uLon << u->lon[j] * 1 / radConv << '\t';
        for (int i = 0; i < u->ReturnFieldLatLen(); i++) uLat << u->lat[i] * 1 / radConv << '\t';

        for (int j = 0; j < v->ReturnFieldLonLen(); j++) vLon << v->lon[j] * 1 / radConv << '\t';
        for (int i = 0; i < v->ReturnFieldLatLen(); i++) vLat << v->lat[i] * 1 / radConv << '\t';
    }
    else if (out_num == -2) {
        //WHAT IS THIS FOR!?
    }
    else {
        if (energy->converged) DumpFields(out_num);
        // else if (out_num % 5 == 0) DumpFields(out_num); //dump every 5 orbits

        // DumpFields(out_num); //FOR DEBUGGING
    }
};

void Solver::ReadInitialConditions(bool yes) {

    if (yes) {
        hid_t init_file = H5Fopen("InitialConditions/initial_conditions.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

        hid_t displacement_h5 = H5Dopen(init_file, "displacement", H5P_DEFAULT);
        hid_t northVel_h5 = H5Dopen(init_file, "north velocity", H5P_DEFAULT);
        hid_t eastVel_h5 = H5Dopen(init_file, "east velocity", H5P_DEFAULT);

        H5Dread(displacement_h5, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, eta_1D);
        H5Dread(northVel_h5, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, v_1D);
        H5Dread(eastVel_h5, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, u_1D);
    }

    for (int i=0; i< etaLatLen; i++) {
        for (int j=0; j<etaLonLen; j++) {
            if (yes) etaOldArray[i][j] = eta_1D[i*etaLonLen + j];
            else etaOldArray[i][j] = 0.0;
        }
    }

    for (int i=0; i< vLatLen; i++) {
        for (int j=0; j<vLonLen; j++) {
            if (yes) vOldArray[i][j] = v_1D[i*vLonLen + j];
            else vOldArray[i][j] = 0.0;
        }
    }

    for (int i=0; i< uLatLen; i++) {
        for (int j=0; j<uLonLen; j++) {
            if (yes) uOldArray[i][j] = u_1D[i*uLonLen + j];
            else uOldArray[i][j] = 0.0;
        }
    }

    double h_thick = consts->h.Value();

    for (int i = 0; i < depth->ReturnFieldLatLen(); i++) {
        for (int j = 0; j < depth->ReturnFieldLonLen(); j++) {
            depthArray[i][j] = h_thick;
        }
    }
};

void Solver::DumpFields(int output_num) {
    std::string out = std::to_string(output_num);



    double factor = 0;

    switch (consts->fric_type) {
    case LINEAR:
        factor = consts->alpha.Value();
        break;

    case QUADRATIC:
        factor = consts->alpha.Value()/consts->h.Value();
        break;
    }


    //-----------------Unravel solutions into 1D arrays---------------------------

    for (int i=0; i<etaLatLen; i++) {
        for (int j=0; j<etaLonLen; j++) {
            eta_1D[i*etaLonLen + j] = (float)etaNewArray[i][j];
        }
    }

    for (int i=0; i<uLatLen; i++) {
        for (int j=0; j<uLonLen; j++) {
            u_1D[i*uLonLen + j] = (float)uNewArray[i][j];
        }
    }

    for (int i=0; i<vLatLen; i++) {
        for (int j=0; j<vLonLen; j++) {
            v_1D[i*vLonLen + j] = (float)vNewArray[i][j];
        }
    }

    for (int i=0; i<etaLatLen; i++) {
        for (int j=0; j<etaLonLen; j++) {
            diss_1D[i*etaLonLen + j] = (float)energy->solution[i][j]*2.0*factor;
        }
    }


    for (int i=0; i<l_max+1; i++)
    {
        for (int j=0; j<l_max+1; j++)
        {
            harm_coeff_1D[i*harm_cols + j] = (float)SH_cos_coeff[i][j];
            // harm_coeff_1D[i*harm_cols + j] = (float)SH_cos_coeff[i*(l_max+1) + j];
        }
    }

    for (int i=0; i<l_max+1; i++)
    {
        for (int j=0; j<l_max+1; j++)
        {
            harm_coeff_1D[(l_max+1)*(l_max+1) + i*harm_cols + j] = (float)SH_sin_coeff[i][j];
            // harm_coeff_1D[(l_max+1)*(l_max+1) + i*harm_cols + j] = (float)SH_sin_coeff[i*(l_max+1) + j];
        }
    }


    diss_avg_1D[0] = energy->currentDissEAvg;

    start[0] = output_num-1;
    start[1] = 0;
    start[2] = 0;

    count[0] = 1;

    start_1D[0] = output_num - 1;

    start_harm[0] = output_num - 1;

    // count 1D is already set to 1.

    //----------------------------------------------------------------------------
    //------------------------WRITE ARRAYS TO DATA FILE---------------------------
    //----------------------------------------------------------------------------


    // ----------------------- Write displacement field --------------------------

    if (consts->field_displacement_output.Value()) {
        count[1] = eta_rows;
        count[2] = eta_cols;

        H5Sselect_hyperslab(data_space_eta, H5S_SELECT_SET, start, NULL, count, NULL);

        H5Dwrite(data_set_eta, H5T_NATIVE_FLOAT, mem_space_eta, data_space_eta, H5P_DEFAULT, eta_1D);
    }


    // ----------------------- Write velocity fields -----------------------------

    if (consts->field_velocity_output.Value()) {

        // ----------------------- Write north velocity field ------------------------

        count[1] = v_rows;
        count[2] = v_cols;

        H5Sselect_hyperslab(data_space_v, H5S_SELECT_SET, start, NULL, count, NULL);

        H5Dwrite(data_set_v, H5T_NATIVE_FLOAT, mem_space_v, data_space_v, H5P_DEFAULT, v_1D);

        // ----------------------- Write east velocity field -------------------------

        count[1] = u_rows;
        count[2] = u_cols;

        H5Sselect_hyperslab(data_space_u, H5S_SELECT_SET, start, NULL, count, NULL);

        H5Dwrite(data_set_u, H5T_NATIVE_FLOAT, mem_space_u, data_space_u, H5P_DEFAULT, u_1D);
    }

    // ----------------------- Write dissipated energy field ---------------------
    if (consts->field_diss_output.Value()) {

        count[1] = eta_rows;
        count[2] = eta_cols;

        H5Sselect_hyperslab(data_space_diss, H5S_SELECT_SET, start, NULL, count, NULL);

        H5Dwrite(data_set_diss, H5T_NATIVE_FLOAT, mem_space_diss, data_space_diss, H5P_DEFAULT, diss_1D);
    }

    // ------------------- Write global average dissipation ----------------------
    if (consts->diss.Value()) {

        H5Sselect_hyperslab(data_space_1D_avg, H5S_SELECT_SET, start_1D, NULL, count_1D, NULL);

        H5Dwrite(data_set_1D_avg, H5T_NATIVE_FLOAT, mem_space_1D_avg, data_space_1D_avg, H5P_DEFAULT, diss_avg_1D);
    }

    // ----------------------- Write SH coefficients ------------------------
    if (consts->sh_coeff_output.Value()) {
        H5Sselect_hyperslab(data_space_harm_coeff, H5S_SELECT_SET, start_harm, NULL, count_harm, NULL);

        H5Dwrite(data_set_harm_coeff, H5T_NATIVE_FLOAT, mem_space_harm_coeff, data_space_harm_coeff, H5P_DEFAULT, harm_coeff_1D);
    }

};

void Solver::CreateHDF5FrameWork(void) {
  #if _WIN32
    mkdir(&(consts->path + SEP + "Grid" + SEP)[0], NULL);

    mkdir(&(consts->path +  SEP + "DATA" + SEP)[0], NULL);

  #elif __linux__
    mkdir(&(consts->path +  SEP + "Grid" + SEP)[0], S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP);

    mkdir(&(consts->path +  SEP + "DATA" + SEP)[0], S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP);

  #endif


    eta_1D = new float[etaLatLen*etaLonLen];
    v_1D = new float[vLatLen*vLonLen];
    u_1D = new float[uLatLen*uLonLen];
    diss_1D = new float[etaLatLen*etaLonLen];
    diss_avg_1D = new float[1];
    harm_coeff_1D = new float[2*(l_max+1)*(l_max+1)];

    eta_rows = etaLatLen;
    eta_cols = etaLonLen;
    v_rows = vLatLen;
    v_cols = vLonLen;
    u_rows = uLatLen;
    u_cols = uLonLen;
    harm_rows = l_max+1;
    harm_cols = l_max+1;

    time_slices = (consts->endTime.Value()/consts->period.Value())/consts->outputTime.Value() + 1;

    rank_field = 3;
    hsize_t rank_size = 1;
    rank_1D = 1;
    rank_harm = 4;

    // char dataFile[] = "DATA/data.h5";
    char dataFile[1024];
    std::strcpy(dataFile, (Out->dataPath).c_str());
    // std::cout << "HERE" << '\n';


    start = new hsize_t[3];
    count = new hsize_t[3];

    start_1D = new hsize_t[1];
    count_1D = new hsize_t[1];

    start_harm = new hsize_t[4];
    count_harm = new hsize_t[4];

    // Create HDF5 file
    file = H5Fcreate(dataFile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // std::cout << "HERE" << '\n';


    dims_eta = new hsize_t[3];
    dims_eta[0] = 1;
    dims_eta[1] = eta_rows;
    dims_eta[2] = eta_cols;

    dims_u = new hsize_t[3];
    dims_u[0] = 1;
    dims_u[1] = u_rows;
    dims_u[2] = u_cols;

    dims_v = new hsize_t[3];
    dims_v[0] = 1;
    dims_v[1] = v_rows;
    dims_v[2] = v_cols;

    dims_1D_avg = new hsize_t[1];
    dims_1D_avg[0] = 1;

    dims_harm_coeff = new hsize_t[4];
    dims_harm_coeff[0] = 1;
    dims_harm_coeff[1] = 2;
    dims_harm_coeff[2] = harm_rows;
    dims_harm_coeff[3] = harm_cols;

    max_dims_eta = new hsize_t[3];
    max_dims_eta[0] = time_slices;
    max_dims_eta[1] = eta_rows;
    max_dims_eta[2] = eta_cols;

    max_dims_u = new hsize_t[3];
    max_dims_u[0] = time_slices;
    max_dims_u[1] = u_rows;
    max_dims_u[2] = u_cols;

    max_dims_v = new hsize_t[3];
    max_dims_v[0] = time_slices;
    max_dims_v[1] = v_rows;
    max_dims_v[2] = v_cols;

    max_dims_1D_avg = new hsize_t[1];
    max_dims_1D_avg[0] = time_slices;

    max_dims_harm_coeff = new hsize_t[4];
    max_dims_harm_coeff[0] = time_slices;
    max_dims_harm_coeff[1] = 2;
    max_dims_harm_coeff[2] = harm_rows;
    max_dims_harm_coeff[3] = harm_cols;

    if (consts->field_displacement_output.Value()) {
        data_space_eta = H5Screate_simple(rank_field, max_dims_eta, NULL); // 3D data space
        data_set_eta = H5Dcreate(file, "displacement", H5T_NATIVE_FLOAT, data_space_eta, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        mem_space_eta = H5Screate_simple(rank_field, dims_eta, NULL);
    }

    if (consts->field_velocity_output.Value()) {
        data_space_u = H5Screate_simple(rank_field, max_dims_u, NULL);
        data_space_v = H5Screate_simple(rank_field, max_dims_v, NULL);
        data_set_u = H5Dcreate(file, "east velocity", H5T_NATIVE_FLOAT, data_space_u, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        data_set_v = H5Dcreate(file, "north velocity", H5T_NATIVE_FLOAT, data_space_v, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        mem_space_u = H5Screate_simple(rank_field, dims_u, NULL);
        mem_space_v = H5Screate_simple(rank_field, dims_v, NULL);
    }

    if (consts->field_diss_output.Value()) {
        data_space_diss = H5Screate_simple(rank_field, max_dims_eta, NULL);
        data_set_diss = H5Dcreate(file, "dissipated energy", H5T_NATIVE_FLOAT, data_space_diss, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        mem_space_diss = H5Screate_simple(rank_field, dims_eta, NULL);
    }

    if (consts->diss.Value()) {
        data_space_1D_avg = H5Screate_simple(rank_1D, max_dims_1D_avg, NULL);
        data_set_1D_avg = H5Dcreate(file, "dissipated energy avg", H5T_NATIVE_FLOAT, data_space_1D_avg, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        mem_space_1D_avg = H5Screate_simple(rank_1D, dims_1D_avg, NULL);
    }

    if (consts->sh_coeff_output.Value()) {
        data_space_harm_coeff = H5Screate_simple(rank_harm, max_dims_harm_coeff, NULL);
        data_set_harm_coeff = H5Dcreate(file, "sh coefficients", H5T_NATIVE_FLOAT, data_space_harm_coeff, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        mem_space_harm_coeff = H5Screate_simple(rank_harm, dims_harm_coeff, NULL);
    }

    start[0] = 0;
    start[1] = 0;
    start[2] = 0;

    count[0] = 1;
    count[1] = eta_rows;
    count[2] = eta_cols;

    start_1D[0] = 0;
    count_1D[0] = 1;

    start_harm[0] = 0;
    start_harm[1] = 0;
    start_harm[2] = 0;
    start_harm[3] = 0;

    count_harm[0] = 1;
    count_harm[1] = 2;
    count_harm[2] = harm_rows;
    count_harm[3] = harm_cols;

    //-----------------------Write dataset attributes---------------------------

    hsize_t dims_coord[1];
    dims_coord[0] = 2;

    hid_t attr_space;
    hid_t attr;

    int data_set_size[2];
    float data_set_dx[2];


    //---------------------------Displacement-----------------------------------
    if (consts->field_displacement_output.Value()) {
        attr_space = H5Screate_simple(rank_size, dims_coord, NULL);

        attr = H5Acreate(data_set_eta, "nlat; nlon", H5T_NATIVE_INT, attr_space, H5P_DEFAULT, H5P_DEFAULT);

        data_set_size[0] = etaLatLen;
        data_set_size[1] = etaLonLen;

        H5Awrite(attr, H5T_NATIVE_INT, data_set_size);

        attr = H5Acreate(data_set_eta, "dlat; dlon", H5T_NATIVE_FLOAT, attr_space, H5P_DEFAULT, H5P_DEFAULT);

        data_set_dx[0] = (float)etadLat*1.0/radConv;
        data_set_dx[1] = (float)etadLon*1.0/radConv;

        H5Awrite(attr, H5T_NATIVE_FLOAT, data_set_dx);
    }

    //---------------------------Velocity---------------------------------------
    if (consts->field_velocity_output.Value()) {
        attr = H5Acreate(data_set_v, "nlat; nlon", H5T_NATIVE_INT, attr_space, H5P_DEFAULT, H5P_DEFAULT);

        data_set_size[0] = vLatLen;
        data_set_size[1] = vLonLen;

        H5Awrite(attr, H5T_NATIVE_INT, data_set_size);

        attr = H5Acreate(data_set_v, "dlat; dlon", H5T_NATIVE_FLOAT, attr_space, H5P_DEFAULT, H5P_DEFAULT);

        data_set_dx[0] = (float)vdLat*1.0/radConv;
        data_set_dx[1] = (float)vdLon*1.0/radConv;

        H5Awrite(attr, H5T_NATIVE_FLOAT, data_set_dx);

        attr = H5Acreate(data_set_u, "nlat; nlon", H5T_NATIVE_INT, attr_space, H5P_DEFAULT, H5P_DEFAULT);

        data_set_size[0] = uLatLen;
        data_set_size[1] = uLonLen;

        H5Awrite(attr, H5T_NATIVE_INT, data_set_size);

        attr = H5Acreate(data_set_u, "dlat; dlon", H5T_NATIVE_FLOAT, attr_space, H5P_DEFAULT, H5P_DEFAULT);

        data_set_dx[0] = (float)udLat*1.0/radConv;
        data_set_dx[1] = (float)udLon*1.0/radConv;

        H5Awrite(attr, H5T_NATIVE_FLOAT, data_set_dx);
    }

    //-------------------------Dissipated energy----------------------------------
    if (consts->field_diss_output.Value()) {
        attr = H5Acreate(data_set_diss, "nlat; nlon", H5T_NATIVE_INT, attr_space, H5P_DEFAULT, H5P_DEFAULT);

        data_set_size[0] = etaLatLen;
        data_set_size[1] = etaLonLen;

        H5Awrite(attr, H5T_NATIVE_INT, data_set_size);

        attr = H5Acreate(data_set_diss, "dlat; dlon", H5T_NATIVE_FLOAT, attr_space, H5P_DEFAULT, H5P_DEFAULT);

        data_set_dx[0] = (float)etadLat*1.0/radConv;
        data_set_dx[1] = (float)etadLon*1.0/radConv;

        H5Awrite(attr, H5T_NATIVE_FLOAT, data_set_dx);
    }

    // delete[] dataFile;

}
