#include "mesh.h"
#include "globals.h"
#include "array3d.h"
#include "array2d.h"
#include "array1d.h"
#include "spatialOperators.h"
#include <math.h>
#include <iostream>
#include <mkl.h>

void setBeta(Array1D<double> & beta, int m, int n, int node_num, Array3D<double> & trigMLon, Array3D<double> & trigNLat)
{
    int i;

    for (i = 0 ; i < node_num; i++)
    {
        beta(i) = 1e6*trigMLon(i, m, 0) * trigNLat(i, n, 0) * trigNLat(i, n, 0)
                  * trigNLat(i, n, 0) * trigNLat(i, n, 0);
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

    gradBeta    = new Array2D<double>(node_num, 2);
    gradTest    = new Array2D<double>(node_num, 2);
    laplaceBeta = new Array1D<double>(node_num);
    divU        = new Array1D<double>(node_num);
    divTest     = new Array1D<double>(node_num);


    trigLon     = &(mesh->trigLon);
    trigLat     = &(mesh->trigLat);
    trigSqLon   = &(mesh->trigSqLon);
    trigSqLat   = &(mesh->trigSqLat);
    trigMLon    = &(mesh->trigMLon);
    trigNLat    = new Array3D<double>(node_num, nMax+1, 2);

    for (i=0; i < node_num; i++)
    {
        for (n=0; n < nMax+1; n++)
        {
            (*trigNLat)(i,n,0) = cos(mesh->node_pos_sph(i,0));
            (*trigNLat)(i,n,1) = sin(mesh->node_pos_sph(i,0));
        }
    }

    int a = 0;
    double sum = -1.0;
    for (n=nMin; n < nMax + 1; n++)
    {
        for (m = mMin; m < mMax + 1; m++)
        {
            setBeta(*beta, m, n, node_num, *trigMLon, *trigNLat);

            setU(*u, m, n, node_num, *trigMLon, *trigNLat, *trigLon, *trigLat);

            setDivU(*divU, m, n, node_num, r, *trigMLon, *trigNLat, *trigLon, *trigLat);

            setGradBeta(*gradBeta, m, n, node_num, r, *trigMLon, *trigNLat, *trigLat);

            velocityDivergence(mesh, *divTest, *u, sum, -1.0);

            // pressureGradient(mesh, *gradTest, *beta, node_num, -1.0);

            sparse_operation_t operation = SPARSE_OPERATION_NON_TRANSPOSE;
            matrix_descr descript;
            // struct matrix_descr {
            //     sparse_matrix_type_t type;
            // } descript;
            descript.type = SPARSE_MATRIX_TYPE_GENERAL;
            double alpham = 1.0;
            double betam = 1.0;
            int error;

            // double * u_vel, * v_vel;
            // u_vel = new double[node_num];
            // v_vel = new double[node_num];
            double * du_vel, * dv_vel;
            du_vel = new double[2*node_num];
            // dv_vel = new double[node_num];
            double * vec = new double[3*node_num];
            for (i=0; i<node_num; i++)
            {
                // dvdt(i,0) = 0.0;
                // dvdt(i,1) = 0.0;
                // pressure(i) = 0.001*i;
                // u_vel[i] = velocity(i,0);
                // v_vel[i] = velocity(i,1);
                du_vel[2*i] = 0.0;
                du_vel[2*i+1] = 0.0;
                // dv_vel[i] = dvdt(i,1);
                vec[3*i] = 0.0;
                vec[3*i+1] = 0.0;
                vec[3*i+2] = (*beta)(i);
            }

            error = mkl_sparse_d_mv(operation, -1.0, *(mesh->operatorGradientX), descript, vec, betam, du_vel);

            for (i = 0; i < node_num; i++)
            {
                std::cout<<i<<' ';
                std::cout<<(*gradBeta)(i,0)<<'\t';
                // std::cout<<(*gradBeta)(i,1)<<'\t';
                // std::cout<<(*u)(i,a)/(r*(*trigLon)(i,1))<<'\t';
                // std::cout<<(*gradTest)(i,0)<<'\t';
                // std::cout<<(*gradTest)(i,1)<<'\t';
                std::cout<<du_vel[2*i]<<'\t';
                // std::cout<<(*gradTest)(i,1)<<'\t';
                // std::cout<<fabs((*gradBeta)(i,a) - (*gradTest)(i,a))/(*gradBeta)(i,a)*100.0<<std::endl;

                std::cout<<std::endl;
                // std::cout<<(*divU)(i)<<'\t';
                // std::cout<<(*divTest)(i)<<std::endl;
                // std::cout<<(*divU)(i) - (*divTest)(i)<<std::endl;


                (*gradTest)(i,0) = 0.0;
                (*gradTest)(i,1) = 0.0;
                (*divTest)(i) = 0.0;
            }

                // set u vector
                // set laplacian solution
                // set div solution
                // set grad solution
        }
    }


    // Define vector u

    // Define analytical solutions to grad beta, laplacian beta, and div u
};
