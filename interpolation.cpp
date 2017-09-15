#include "mesh.h"
#include "globals.h"
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include <math.h>
#include <iostream>

/**
*   @purpose    Function interpolates data from the geodesic grid to the regular
*               lat-lon grid. Interpolation is biquadratic and takes the form
*               Vc = d, where V is a 6x6 vandemonde matrix, c is a 6x1 column vector
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
*
*/
int interpolateGG2LL(Globals * globals,
                     Mesh * mesh,
                     Array2D<double> & ll_data,
                     Array1D<double> & gg_data,
                     Array3D<double> & ll_map_coords,
                     Array3D<double> & V_inv,
                     Array2D<int> & cell_ID)

{
    int i, j, k, l, f, N_ll, gg_ID, f_num;
    double x, y;
    double d_vec[6];            // static data vector
    double c_vec[6];            // static coefficient vector
    double V_inv_mat[6][6];     // static Vandermonde inverse matrix

    Array2D<int> * friend_list;
    friend_list = &(mesh->node_friends);

    N_ll = globals->dLat.Value();

    for (i = 0; i < 180/N_ll; i++)
    {
        for (j = 0; j < 360/N_ll; j++)
        {
            /*------------------------------------------------------------------
            * STEP 1: FIND THE GEODESIC CELL THAT LAT-LON NODE i,j BELONGS TO
            * ------------------------------------------------------------------
            */

            // get geodesic cell ID of lat-lon node i,j
            gg_ID = cell_ID( i, j );


            /*------------------------------------------------------------------
            * STEP 2: FILL DATA VECTOR WITH DATA FROM THE GEODESIC GRID NODES
            * ------------------------------------------------------------------
            */

            // find number of friends
            f_num = 6;
            if ((*friend_list)(gg_ID, 5) == -1)
            {
                f_num = 5;

                // If cell is pentagon, use middle node for interpolation
                d_vec[5] = gg_data(gg_ID);
            }

            // loop through friend nodes,
            for (k = 0; k < f_num; k++)
            {
                // get friend ID
                f = (*friend_list)(gg_ID, k);

                // assign data vector to values of surrounding friends
                d_vec[k] = gg_data(f);
            }


            /*------------------------------------------------------------------
            * STEP 3: FILL VANDERMONDE INVERSE MATRIX
            * ------------------------------------------------------------------
            */

            for (k = 0; k < 6; k++)
            {
                // Unrolling inner loop. Potential speed benefit?
                V_inv_mat[k][0] = V_inv(gg_ID, k, 0);
                V_inv_mat[k][1] = V_inv(gg_ID, k, 1);
                V_inv_mat[k][2] = V_inv(gg_ID, k, 2);
                V_inv_mat[k][3] = V_inv(gg_ID, k, 3);
                V_inv_mat[k][4] = V_inv(gg_ID, k, 4);
                V_inv_mat[k][5] = V_inv(gg_ID, k, 5);
            }

            /*------------------------------------------------------------------
            * STEP 4: PERFORM MATRIX VECTOR DOT PRODUCT
            * ------------------------------------------------------------------
            */

            // This naive implementation of matrix multiplication may be
            // sufficiently fast for the small 6x6 matricies used here
            for (k = 0; k < 6; k++)
            {
                c_vec[k] = 0.0;
                for (l = 0; l < 6; l++)
                {
                    c_vec[k] += V_inv_mat[k][l] * d_vec[l];
                }
            }

            /*------------------------------------------------------------------
            * STEP 5: ASSEMBLE THE INTERPOLATED DATA FROM C_VEC
            * ------------------------------------------------------------------
            */

            // get map coords of current lat-lon node relative to geodesic cell
            x = ll_map_coords(i, j, 0);
            y = ll_map_coords(i, j, 1);

            ll_data(i, j) = c_vec[0]
                            + c_vec[1]*x
                            + c_vec[2]*x*x
                            + c_vec[3]*y
                            + c_vec[4]*x*y
                            + c_vec[5]*y*y;
            //
            // std::cout<<gg_ID<<' '<<x<<' '<<y<<std::endl;
            // std::cout<<c_vec[0]<<' '<<c_vec[1]<<' '<<c_vec[2]<<' '<<c_vec[3]<<' '<<c_vec[4]<<' '<<c_vec[5]<<std::endl;
            // std::cout<<V_inv_mat[0][0]<<' '<<V_inv_mat[0][1]<<' '<<V_inv_mat[0][2]<<' '<<V_inv_mat[0][3]<<' '<<V_inv_mat[0][4]<<' '<<V_inv_mat[0][5]<<std::endl;
            // std::cout<<V_inv_mat[1][0]<<' '<<V_inv_mat[1][1]<<' '<<V_inv_mat[1][2]<<' '<<V_inv_mat[1][3]<<' '<<V_inv_mat[1][4]<<' '<<V_inv_mat[1][5]<<std::endl;
            // std::cout<<V_inv_mat[2][0]<<' '<<V_inv_mat[2][1]<<' '<<V_inv_mat[2][2]<<' '<<V_inv_mat[2][3]<<' '<<V_inv_mat[2][4]<<' '<<V_inv_mat[2][5]<<std::endl;
            // std::cout<<V_inv_mat[3][0]<<' '<<V_inv_mat[3][1]<<' '<<V_inv_mat[3][2]<<' '<<V_inv_mat[3][3]<<' '<<V_inv_mat[3][4]<<' '<<V_inv_mat[3][5]<<std::endl;
            // std::cout<<V_inv_mat[4][0]<<' '<<V_inv_mat[4][1]<<' '<<V_inv_mat[4][2]<<' '<<V_inv_mat[4][3]<<' '<<V_inv_mat[4][4]<<' '<<V_inv_mat[4][5]<<std::endl;
            // std::cout<<V_inv_mat[5][0]<<' '<<V_inv_mat[5][1]<<' '<<V_inv_mat[5][2]<<' '<<V_inv_mat[5][3]<<' '<<V_inv_mat[5][4]<<' '<<V_inv_mat[5][5]<<std::endl;
            // //
            // char c;
            //
            // std::cin>>c;

            // dgemv_(&TRANS, &M, &N, &ALPHA, V_inv_mat, &LDA, d_vec, &INCX, &BETA, c_vec, &INCY);

        }
    }

    return 1;
};
