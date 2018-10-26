#include "mesh.h"
#include "globals.h"
#include "array3d.h"
#include "array2d.h"
#include "array1d.h"
#include "spatialOperators.h"
#include <math.h>
#include <iostream>

#include <mkl.h>
#include <mkl_spblas.h>

void setBeta(Array1D<double> & beta, int m, int n, int node_num, Array3D<double> & trigMLon, Array3D<double> & trigNLat)
{
    int i;

    for (i = 0 ; i < node_num; i++)
    {
        beta(i) = 1e6*trigMLon(i, m, 0) * pow(trigNLat(i, n, 0), 4.0);
    }
};

void setU(Array2D<double> & u, int m, int n, int node_num, Array3D<double> & trigMLon, Array3D<double> & trigNLat, Array2D<double> & trigLon, Array2D<double> & trigLat)
{
    int i;

    for (i = 0 ; i < node_num; i++)
    {
        u(i, 0) = -4. * 1e6* n * trigLon(i, 1) * trigMLon(i, m, 0);
        u(i, 0) *= trigNLat(i, n, 1);
        u(i, 0) *= trigNLat(i, n, 0)*trigNLat(i, n, 0)*trigNLat(i, n, 0);


        u(i, 1) = -m * 1e6*trigLon(i, 1) * trigMLon(i, m, 1);
        u(i, 1) *= trigNLat(i, n, 0)*trigNLat(i, n, 0)
                   *trigNLat(i, n, 0)*trigNLat(i, n, 0)/trigLat(i, 0);

    }
};

void setDivU(Array1D<double> & divU, int m, int n, int node_num, double r, Array3D<double> & trigMLon, Array3D<double> & trigNLat, Array2D<double> & trigLon, Array2D<double> & trigLat)
{
    int i;

    for (i = 0 ; i < node_num; i++)
    {
        divU(i)  = 2.*m*trigLon(i, 1)*trigMLon(i, m, 1);
        divU(i) -= trigLon(i, 0)*trigMLon(i, m, 0);
        divU(i) *= trigNLat(i, n, 1);
        divU(i) *= trigNLat(i, n, 0)*trigNLat(i, n, 0)*trigNLat(i, n, 0);
        divU(i) *= 1e6*4.*n/(r * trigLat(i, 0));
    }
};

void setLaplaceBeta(Array1D<double> & lapBeta, int m, int n, int node_num, double r, Array3D<double> & trigMLon, Array3D<double> & trigNLat, Array2D<double> & trigLon, Array2D<double> & trigLat)
{
    int i;

    double dn = (double)n;
    double dm = (double)m;
    double sinNLat, sinLat;
    double sin2NLat;
    double tanNLat;
    double cosNLat, cosLat;
    double cosMLon, sinMLon;


    for (i = 0 ; i < node_num; i++)
    {
        cosNLat = trigNLat(i,n,0);
        sinNLat = trigNLat(i,n,1);
        sinNLat = trigNLat(i,2*n,1);
        cosLat  = trigNLat(i,1,0);
        sinLat  = trigNLat(i,1,1);
        cosMLon = trigMLon(i,m,0);
        sinMLon = trigMLon(i,m,1);
        tanNLat = sinNLat/cosNLat;

        // lapBeta(i)  = trigNLat(i, 1, 0)*(dn - 2*dn * trigNLat(i, 2*n, 0));
        // lapBeta(i) += trigNLat(i, 1, 1)*trigNLat(i, n, 1)*trigNLat(i, n, 0);
        // lapBeta(i) *= 4*dn*trigNLat(i, n, 0)*trigNLat(i, n, 0)*trigMLon(i, m, 0);
        // lapBeta(i) *= 1e12/(r*r*trigNLat(i, 1, 0));

        // lapBeta(i)  = dn*cosLat*pow(cosNLat, 4.0) - sinLat*sinNLat*pow(cosNLat, 3.0);
        // lapBeta(i) -= 3*dn*cosLat*pow(sinNLat, 2.0)*pow(cosNLat, 2.0);
        //
        // lapBeta(i) *= -1e6*4*n*cosMLon/(r*r*cosLat);

        lapBeta(i)  = dm*dm/cosLat + 4.*dn*dn*cosLat - 12.*dn*dn*cosLat*pow(tanNLat, 2.0) - 4.*dn*sinLat*tanNLat;
        lapBeta(i) *= -1e6*cosMLon*pow(cosNLat, 4.0)/(r*r*cosLat);

        // lapBeta(i) = trigNLat(i, 1, 1)*trigNLat(i, n, 1)*trigNLat(i, n, 0);
        // lapBeta(i) -= n*trigNLat(i, 1, 0)*(pow(trigNLat(i, n, 0), 2.0) - 3*pow(trigNLat(i, n, 1), 2.0));
        // lapBeta(i) *= 4 * n * pow(trigNLat(i, n, 0), 2.0) * trigMLon(i, m, 0);


        // lapBeta(i) += -1e12/(r*r*trigNLat(i, 1, 0)*trigNLat(i, 1, 0)) * (dm*dm*trigMLon(i, m, 0)*trigMLon(i, m, 0)*pow(trigNLat(i, n, 0), 4.0));
        // lapBeta(i) *= -1e6/(r*r*trigNLat(i, 1, 0)*trigNLat(i, 1, 0));
    }
};

void setGradBeta(Array2D<double> & gradBeta, int m, int n, int node_num, double r, Array3D<double> & trigMLon, Array3D<double> & trigNLat, Array2D<double> & trigLat)
{
    int i;

    for (i = 0 ; i < node_num; i++)
    {
        gradBeta(i, 1) =  -1e6 * 4.0 * n * trigMLon(i, m, 0);
        gradBeta(i, 1) *= trigNLat(i, n, 1);
        gradBeta(i, 1) *= trigNLat(i, n, 0)*trigNLat(i, n, 0)*trigNLat(i, n, 0);
        gradBeta(i, 1) /= r;

        gradBeta(i, 0) = -m * 1e6 * trigMLon(i, m, 1);
        gradBeta(i, 0) *= trigNLat(i, n, 0)*trigNLat(i, n, 0)
                   *trigNLat(i, n, 0)*trigNLat(i, n, 0);
        gradBeta(i, 0) /= r*trigLat(i, 0);
    }
};

void runOperatorTests(Globals * globals, Mesh * mesh)
{
    int node_num;
    int i, m, n;
    int mMax = 1;
    int mMin = 1;
    int nMax = 1;
    int nMin = 1;

    double r;

    node_num = mesh->node_num;
    r = globals->radius.Value();

    // Define alpha and beta test functions
    Array1D<double> * alpha;
    Array1D<double> * beta;
    Array2D<double> * u;

    Array2D<double> * gradBeta;
    Array2D<double> * gradTest;
    Array1D<double> * laplaceBeta;
    Array1D<double> * laplaceTest;
    Array1D<double> * divU;
    Array1D<double> * divTest;

    Array2D<double> * trigLon;
    Array2D<double> * trigLat;
    Array2D<double> * trigSqLon;
    Array2D<double> * trigSqLat;
    Array3D<double> * trigMLon;
    Array3D<double> * trigNLat;

    alpha       = new Array1D<double>(node_num);
    beta        = new Array1D<double>(node_num);
    u           = new Array2D<double>(node_num, 2);

    double * beta2;
    beta2 = new double[node_num];
    double * ux;
    ux = new double[node_num];
    double * uy;
    uy = new double[node_num];
    double * laplace2;
    laplace2 = new double[node_num];
    double * laplacex;
    laplacex = new double[node_num];
    double * laplacey;
    laplacey = new double[node_num];

    gradBeta    = new Array2D<double>(node_num, 2);
    gradTest    = new Array2D<double>(node_num, 2);
    laplaceBeta = new Array1D<double>(node_num);
    laplaceTest = new Array1D<double>(node_num);
    divU        = new Array1D<double>(node_num);
    divTest     = new Array1D<double>(node_num);


    trigLon     = &(mesh->trigLon);
    trigLat     = &(mesh->trigLat);
    trigSqLon   = &(mesh->trigSqLon);
    trigSqLat   = &(mesh->trigSqLat);
    trigMLon    = &(mesh->trigMLon);
    trigNLat    = new Array3D<double>(node_num, 2*(nMax+1), 2);

    for (i=0; i < node_num; i++)
    {
        for (n=0; n < 2*(nMax+1); n++)
        {
            (*trigNLat)(i,n,0) = cos(mesh->node_pos_sph(i,0));
            (*trigNLat)(i,n,1) = sin(mesh->node_pos_sph(i,0));
        }
    }

    int a = 0;
    int error;
    double sum = -1.0;
    for (n=nMin; n < nMax + 1; n++)
    {
        for (m = mMin; m < mMax + 1; m++)
        {
            setBeta(*beta, m, n, node_num, *trigMLon, *trigNLat);
            for (i=0; i < node_num; i++) beta2[i] = (*beta)(i);

            setU(*u, m, n, node_num, *trigMLon, *trigNLat, *trigLon, *trigLat);
            for (i=0; i < node_num; i++)
            {
                ux[i] = (*u)(i,0);
                uy[i] = (*u)(i,1);
            }

            setDivU(*divU, m, n, node_num, r, *trigMLon, *trigNLat, *trigLon, *trigLat);

            // setDivU(*laplaceBeta, m, n, node_num, r, *trigMLon, *trigNLat, *trigLon, *trigLat);

            setGradBeta(*gradBeta, m, n, node_num, r, *trigMLon, *trigNLat, *trigLat);

            setLaplaceBeta(*laplaceBeta, m, n, node_num, r, *trigMLon, *trigNLat, *trigLon, *trigLat);

            // velocityDivergence(mesh, *divTest, *u, sum, -1.0);
            //
            // pressureGradient(mesh, *gradTest, *beta, node_num, -1.0);

            sparse_operation_t operation = SPARSE_OPERATION_NON_TRANSPOSE;
            matrix_descr descript;
            // struct matrix_descr {
            //     sparse_matrix_type_t type;
            // } descript;
            descript.type = SPARSE_MATRIX_TYPE_GENERAL;
            double alpham = 1.0;
            double betam = 0.0;
            error = mkl_sparse_d_mv(operation, alpham, *(mesh->operatorLaplacian), descript, beta2, betam, laplace2);
            // error = mkl_sparse_d_mv(operation, alpham, *(mesh->operatorTest), descript, ux, betam, laplacex);
            // error = mkl_sparse_d_mv(operation, alpham, *(mesh->operatorLaplacian), descript, uy, betam, laplacey);
            std::cout<<"AND THE ERROR IS.... "<<error<<std::endl;

            // pressureGradient(mesh, *gradTest, *beta, node_num, -1.0);
            //
            // for (i=0; i < node_num; i++) {
            //         (*gradTest)(i, 0) *= (*trigMLon)(i, 1, 1);
            //         (*gradTest)(i, 1) *= (*trigMLon)(i, 1, 1);
            // }
            //
            // velocityDivergence(mesh, *laplaceTest, *u, sum, -1.0);
            double mean_err = 0.0;
            for (i = 0; i < node_num; i++)
            {
                // std::cout<<(*gradBeta)(i,0)<<'\t';
                // // std::cout<<(*u)(i,a)/(r*(*trigLon)(i,1))<<'\t';
                // std::cout<<(*gradTest)(i,0)<<'\t';
                // std::cout<<i<<'\t'<<(*gradBeta)(i,1)<<'\t';
                // std::cout<<laplace2[i]<<std::endl;
                // std::cout<<(*gradTest)(i,0)<<std::endl;
                // mean_err += fabs((*gradBeta)(i,0) - (*gradTest)(i,0));

                // std::cout<<(*gradBeta)(i,1)<<'\t';
                // std::cout<<(*gradTest)(i,1)<<std::endl;
                // mean_err += fabs((*gradBeta)(i,1) - (*gradTest)(i,1));

                // std::cout<<i<<'\t'<<(*divU)(i)<<'\t';
                // std::cout<<laplacex[i]+laplacey[i]<<std::endl;
                // std::cout<<(*laplaceTest)(i)<<std::endl;
                // mean_err += fabs((*divU)(i) - (*divTest)(i));

                std::cout<<i<<'\t'<<(*laplaceBeta)(i)<<'\t';
                std::cout<<laplace2[i]<<'\t'<<(*laplaceBeta)(i) - laplace2[i]<<std::endl;
                // mean_err += fabs((*laplaceBeta)(i) - laplace2[i]);
                (*gradTest)(i,0) = 0.0;
                (*gradTest)(i,1) = 0.0;
                (*divTest)(i) = 0.0;
            }
            // std::cout<<mean_err/node_num<<std::endl;

                // set u vector
                // set laplacian solution
                // set div solution
                // set grad solution
        }
    }
};
