#include "mesh.h"
#include "globals.h"
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include <math.h>
#include <iostream>
#include "spatialOperators.h"

// #include <mkl_spblas.h>
// #include <mkl.h>

/**
*   @purpose    Function interpolates data from the geodesic grid to the regular
*               lat-lon grid. Interpolation is biquadratic and takes the form
*               Vc = d, where V is a 6x6 vandermonde matrix, c is a 6x1 column vector
*               of the interpolating function coefficients, and d is a 6x1 column
*               vector of the geodesic grid data to be interpolated. See Lee and
*               Macdonald (2008).
*
*   @params
*
*   *globals        pointer to the Globals object in ODIS
*   *mesh           pointer to the Mesh object in ODIS
*   ll_data         reference to the initialised 2d array for data on the
*                   lat-lon grid
*   gg_data         reference to the initialised geodesic grid data
*   ll_map_coords   reference to the coords of lat-lon nodes in mapping space
*   V_inv           reference to the inverse Vandermonde matrix for every cell
*   cell_ID         reference to the geodesic cell ID for each lat-lon node
*
*   @returns
*
*   ll_data         array with data that has been interpolated from the geodesic
*                   grid to the 2D lat-lon grid
*
*   @author         Hamish Hay
*
*   @history        15/09/2017 - created
*                   05/01/2018 - switched to conservative interpolation
*
*/

int interpolateGG2LLConservative(Globals * globals,
                                 Mesh * mesh,
                                 Array2D<double> & ll_data,
                                 Array1D<double> & gg_data)
{
    double * gg_data_1D;
    double * ll_data_1D;

    int dLat;
    dLat = globals->dLat.Value();

    double r;
    r = globals->radius.Value();

    // No. of nodes on the geodesic and lat-lon grids.
    int node_num_gg = mesh->node_num;
    int node_num_ll = (int)(360/dLat)*(int)(180/dLat);

    // geodesic grid data (0: val, 1: lat gradient, 2: lon gradient)
    gg_data_1D = new double[node_num_gg * 3];

    // lat-lon interpolated solution
    ll_data_1D = new double[node_num_ll] (); // <--- initialize to zero!!!!


    // -------------------------------------------------------------------------
    // Compute gradient of the solution to be interpolated
    // -------------------------------------------------------------------------

    Array2D<double> * gradient;
    gradient = new Array2D<double>(node_num_gg, 2);

    pressureGradient(mesh, *gradient, gg_data, node_num_gg, -1.0);


    // -------------------------------------------------------------------------
    // Fill the data vector with the geodesic grid solution and its gradients
    // -------------------------------------------------------------------------

    double * cosLat = &(mesh->trigLat(0,0));
    for (int i=0; i<node_num_gg; i++)
    {
        gg_data_1D[3*i]     = gg_data(i);
        gg_data_1D[3*i + 1] = r*(*gradient)(i, 1);
        gg_data_1D[3*i + 2] = r*(*gradient)(i, 0);
    }


    // -------------------------------------------------------------------------
    // Define variables for matrix-vector multiplication
    // -------------------------------------------------------------------------

    double alpha = 1.0;
    double beta = 1.0;

    matrix_descr descrp;
    descrp.type = SPARSE_MATRIX_TYPE_GENERAL;

    sparse_operation_t operation = SPARSE_OPERATION_NON_TRANSPOSE;

    sparse_status_t err;

    sparse_matrix_t dest;
    sparse_index_base_t index_type = SPARSE_INDEX_BASE_ZERO;
    int nrows = (360/dLat)*(180/dLat);
    int ncols = 3*node_num_gg;

    // -------------------------------------------------------------------------
    // Multiply the sparse interpolation matrix by the geodesic grid data vector
    // and return the interpolated solution in ll_data_1D
    // -------------------------------------------------------------------------

    err = mkl_sparse_d_mv (operation, alpha, *(mesh->interpMatrix), descrp, gg_data_1D, beta, ll_data_1D);


    // -------------------------------------------------------------------------
    // Convert the 1D interpolated array to 2D. This is actually not necessary,
    // and could be avoided in the future if the SH routines are modified.
    // -------------------------------------------------------------------------

    int count = 0;
    // double lat, lon;
    for (int i=0; i<(int)(180/dLat); i++)
    {
        for (int j=0; j<(int)(360/dLat); j++)
        {
            ll_data(i, j) = ll_data_1D[count];

            count++;
        }
    }

    // globals->Output->TerminateODIS();

    delete[] ll_data_1D;
    delete[] gg_data_1D;
    delete gradient;

    return 1;
}
