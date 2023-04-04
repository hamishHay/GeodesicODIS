#include "mesh.h"
#include "globals.h"
#include "array1d.h"
#include "interpolation.h"
#include "gridConstants.h"

int updateEta(Globals * globals, 
              Mesh * grid, 
              Array1D<double> & deta_dt, 
              Array1D<double> & v_t0, 
              Array1D<double> & eta, 
              Array1D<double> &h_total)
{

    double h = globals->h.Value();
    int advection = globals->advection.Value();

    
    Eigen::Map<Eigen::VectorXd> detadtEigen(&deta_dt(0), grid->node_num_ng);

    
    if (globals->advection.Value()) 
    {
        Array1D<double> hv(FACE_NUM);
        Eigen::Map<Eigen::VectorXd> hvEigen(&hv(0), FACE_NUM);


        // Get thickness flux normal to each face
        interpolateLSQFlux(globals, grid, hv, v_t0, h_total);

        // Get divergence of thickness flux to find dEta/dt
        detadtEigen = grid->operatorDivergence * hvEigen;
    }
    else 
    {
        Eigen::Map<Eigen::VectorXd> hvEigen(&v_t0(0), grid->face_num);

        // Get divergence of velocity field to find dEta/dt
        detadtEigen = h * grid->operatorDivergence * hvEigen;   
    }
     

    return 1;
};
