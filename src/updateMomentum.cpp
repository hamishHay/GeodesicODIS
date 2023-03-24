#include "mesh.h"
#include "globals.h"
#include "outFiles.h"
#include "solver.h"
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include "initialConditions.h"
#include "momAdvection.h"
#include "gridConstants.h"
#include <math.h>
#include <iostream>
#include <sstream>
#include <iomanip>

int updateMomentum(Globals * constants,
                   Mesh * grid, 
                   Array1D<double> & dvdt, 
                   Array1D<double> & v_tm1, 
                   Array1D<double> & p_tm1, 
                   Array1D<double> & h_total, 
                   Array1D<double> & ekin, 
                   double GAMMA=0.0, 
                   double IMPLICIT=0.0)
{
    double g = constants->g.Value();

    Eigen::Map<Eigen::VectorXd> etaEigen(&p_tm1(0), NODE_NUM);
    Eigen::Map<Eigen::VectorXd> velEigen(&v_tm1(0), FACE_NUM);
    Eigen::Map<Eigen::VectorXd> dvdtEigen(&dvdt(0), FACE_NUM);

    // Perform sparse matrix * vector operation
    

    if (constants->advection.Value()) 
    {
        dvdtEigen = -g*(1-GAMMA*IMPLICIT)*grid->operatorGradient * etaEigen;
        calculateMomentumAdvection(constants, grid, dvdt, v_tm1, h_total, ekin);

    }
    else {
        dvdtEigen = -g*(1-GAMMA*IMPLICIT)*grid->operatorGradient * etaEigen + grid->operatorCoriolis * velEigen;
    }
    
    return 1;

};