#include "mesh.h"
#include "globals.h"
#include "math.h"
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include "mathRoutines.h"
#include "gridConstants.h"
// #include "sphericalHarmonics.h"

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cstring>

#include "H5Cpp.h"
#include "mpi.h"

#include <Eigen/Sparse>

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

// Constructor. Calls the series of tasks needed to contruct all information 
// relevant to the grid.
Mesh::Mesh(Globals *Globals, int N, int face_N, int vertex_N, int N_ll, int l_max)
{
    globals = Globals; // define reference to all constants

    // Read in grid file
    // ReadMeshFile(); // Legacy, no longer used.
    
    ReadGridFile();

    CalcTrigFunctions();
    // CalcLegendreFuncs();
    CalcMappingCoords();
    CalcVelocityTransformFactors();
    CalcTangentialVelocityWeights();
    CalcCartesianComponents();

    CalcGradOperatorCoeffs();
    CalcDivOperatorCoeffs();
    CalcCurlOperatorCoeffs();
    CalcCoriolisOperatorCoeffs();
    CalcLinearDragOperatorCoeffs();

    CalcAdjacencyMatrix();
    CalcControlVolumeInterpMatrix();

    CalcRBFInterpMatrix();
    CalcRBFInterpMatrix2();

    CalcFreeSurfaceSolver();

    CalcMaxTimeStep();

    globals->Output->DumpGridData(this); // ?
};

// Function to convert and store normal vectors to cartesian components
int Mesh::CalcCartesianComponents(void)
{
    double normal_vec_xyz[3];
    double normal_vec_sph[2];
    double sph_pos[2];
    
    for (unsigned i=0; i<FACE_NUM; i++) {
        sph_pos[0] = face_centre_pos_sph(i, 0); // latitude 
        sph_pos[1] = face_centre_pos_sph(i, 1); // longitude

        // Note the swapped index here!
        normal_vec_sph[0] = face_normal_vec_map(i, 1);  // Get latitude component
        normal_vec_sph[1] = face_normal_vec_map(i, 0);  // Get longitude component

        // Convert the spherical vector to cartesian
        normalVectorSph2Cart(sph_pos, normal_vec_sph, normal_vec_xyz); // --> overwrites normal_vec_xyz

        face_normal_vec_xyz(i, 0) = normal_vec_xyz[0];
        face_normal_vec_xyz(i, 1) = normal_vec_xyz[1];
        face_normal_vec_xyz(i, 2) = normal_vec_xyz[2];
    }
    return 1;
}

// Function to calculate the w_ee' weights from Eq. 33 Ringler et al (2009).
// This function requires that the face_friends array has ordered the friends
// list in clockwise or anti-clockwise fashion.
int Mesh::CalcTangentialVelocityWeights(void)
{
    double vertex_sub_area = 0.0;
    double weight;
    int fnum;
    int tev, ne;

    int vshared_ID;
    int fstart_ID, face_ID, face_ID2, adj_ID;
    int node_ID;

    // For each face, loop over each neighbouring node. This will give the corresponding weights.
    for (unsigned i=0; i<FACE_NUM; i++)
    {
        // Loop over each node adjacent to face i
        for (unsigned k=0; k<2; k++)
        {
            node_ID = face_nodes(i, k);
            fnum = node_fnum(node_ID);
            tev = 0;

            // Get indicator tev for the vertex shared between 
            // face i and its adjacent face

            adj_ID = face_friends(i, k, 0); // First friend

            // find shared vertex between face i and face adj_ID
            if (face_vertexes(adj_ID, 0) == face_vertexes(i,0)) 
                vshared_ID = face_vertexes(i, 0);
            else if (face_vertexes(adj_ID, 0) == face_vertexes(i,1)) 
                vshared_ID = face_vertexes(i, 1);
            else // face_vertexes(adj_ID, 0) is not shared with face_vertexes(i,:), therefore
                vshared_ID = face_vertexes(adj_ID, 1);

            // Find tev of vertex vshared_ID with face i
            for (unsigned j3=0; j3<3; j3++)
            {
                if (vertex_faces(vshared_ID, j3) == i) 
                {
                    tev = vertex_face_dir(vshared_ID, j3);
                    break;
                }
            }

            // Loop over each face e' to calculate the weights wee'
            for (int j=0; j<fnum-1; j++)
            {
                weight = 0;     // weight w_ee'
            
                fstart_ID = face_friends(i, k, j);

                // This statement assumes that the node order in face_nodes gives
                // the upwind node at face_nodes(:, 0), and the downwind 
                // node at face_nodes(:, 1)
                if (face_nodes(fstart_ID,0) == node_ID) ne = 1; // Node is upwind, face points outwards
                else ne = -1;                                   // Node is downwind of face fstart_ID, face points inwards 
                
                // Now sum up all node->vertex areas, starting from face_ID, to the parent face
                for (int j2 = j; j2>=0; j2--)
                {
                    // Get the current face e' ID
                    face_ID = face_friends(i, k, j2);   // ID of current face friend

                    // Get ID of the face that is next closest to the parent face, i
                    if (j2==0)  face_ID2 = i;                           // Next closest friend is the parent face, i
                    else        face_ID2 = face_friends(i, k, j2 - 1);  // ID of next closest face friend

                    // Find the vertex shared by the pair of faces.
                    if (face_vertexes(face_ID, 0) == face_vertexes(face_ID2,0)) 
                        vshared_ID = face_vertexes(face_ID2, 0);
                    else if (face_vertexes(face_ID, 0) == face_vertexes(face_ID2,1)) 
                        vshared_ID = face_vertexes(face_ID2, 1);
                    else // face_vertexes(face_ID, 0) is not shared with face_vertexes(face_ID2,:), therefore
                        vshared_ID = face_vertexes(face_ID, 1);

                    // Find the area that vertex vshared_ID shares with node node_ID.
                    vertex_sub_area = 0.0;
                    for (unsigned j3=0; j3<fnum; j3++) {
                        // Get the correct area by finding the index j3 where
                        // vshared_ID lives at in the array vertexes (and 
                        // therefore node_vertex_area)
                        
                        if (vertexes(node_ID, j3) == vshared_ID) 
                        { 
                            vertex_sub_area = node_vertex_area(node_ID, j3);
                            break;
                        }
                    }

                    weight += vertex_sub_area*cv_area_sph_r(node_ID);  
                }

                weight = weight - 0.5;  // Eq. 33 Ringler et al
                
                face_interp_weights(i, k, j) = weight * -tev * ne;
            }
        }     
    }

    return 1;
};

int Mesh::CalcControlVolumeInterpMatrix(void)
{
    double r = globals->radius.Value();
    double r_recip = 1.0/r;

    // interpSolvers = new Eigen::ColPivHouseholderQR<Eigen::MatrixXd>[NODE_NUM];
    interpSolvers = new Eigen::FullPivLU<Eigen::MatrixXd>[NODE_NUM];
    // interpVandermonde = new Eigen::MatrixXd[NODE_NUM];

    // Eigen::VectorX< > interpVectors =
    // interpSolvers
    // vandermondeInv = SpMat((NODE_NUM-12)*7 + 12*6, 2);
    interpSolvers2 = new Eigen::LLT<Eigen::MatrixXd>[NODE_NUM];

    // test = Eigen::MatrixXd<Eigen::ColPivHouseholderQR<Eigen::MatrixXd>

    for (int i = 0; i < NODE_NUM; i++)
    {

        int f_num = node_fnum(i);

        Vandermonde V(f_num + 1, 6);

        int j = 0;
        for (j = 0; j < f_num + 1; j++)
        {
            double x, y;

            x = node_pos_map(i, j, 0) * r_recip;
            y = node_pos_map(i, j, 1) * r_recip;

            V(j, 0) = 1.0;
            V(j, 1) = x;
            V(j, 2) = y;
            V(j, 3) = pow(x, 2.0);
            V(j, 4) = x * y;
            V(j, 5) = pow(y, 2.0);
        }

        interpSolvers[i].compute(V);

        // Eigen::HouseholderQR<Eigen::MatrixXd> qr(V.transpose());
        // Eigen::MatrixXd pinv;
        // pinv.setIdentity(V.cols(), V.rows());
        // pinv = qr.householderQ() * pinv;
        // pinv = qr.matrixQR().topLeftCorner(V.rows(),V.rows()).triangularView<Eigen::Upper>().transpose().solve<Eigen::OnTheRight>(pinv);

        // Eigen::MatrixXd B = (V.transpose() * V).completeOrthogonalDecomposition().pseudoInverse() * V.transpose();
        // // interpVandermonde[i] = V;
        // // interpVandermonde.push_back((V.transpose() * V).completeOrthogonalDecomposition().pseudoInverse() * V.transpose() );
        // interpVandermonde.push_back(pinv);
        // if (i< 200) std::cout<<i<<' '<<sizeof(V)<<std::endl;
    }

    {
        int friend_num, f_IDj, f_IDi;
        double sphn[2];
        double sphfi[2], sphfj[2];
        double arc;
        double eps = 0.125; // If you change this, make sure you also change it in CalcRBFInterpMatrix2!!!!
        double x1, x2, y1, y2;
        for (int i = 0; i < NODE_NUM; i++)
        {
            // check if node is pentagon
            friend_num = node_fnum(i);

            DenseMat matrix(friend_num + 1, friend_num + 1);

            sphn[0] = node_pos_sph(i, 0); // Node lat
            sphn[1] = node_pos_sph(i, 1); // Node Lon

            arc = distanceBetweenSph2(sphn, sphn);

            matrix(0, 0) = 1.0;

            // loop through all faces again
            for (int j1 = 0; j1 < friend_num; j1++)
            {
                // Get face ID
                f_IDi = node_friends(i, j1);

                sphfi[0] = node_pos_sph(f_IDi, 0); // Face 2 lat
                sphfi[1] = node_pos_sph(f_IDi, 1); // Face 2 lon

                // Get radial basis function between face centre j1 and j1
                // using constant shapre parameter eps=2
                arc = distanceBetweenSph2(sphn, sphfi);

                x1 = node_pos_map(i, j1 + 1, 0) * r_recip;
                y1 = node_pos_map(i, j1 + 1, 1) * r_recip;
                arc = sqrt(x1 * x1 + y1 * y1);

                matrix(0, j1 + 1) = RBF(arc, eps);
            }

            // loop through all j faces
            for (int j1 = 0; j1 < friend_num; j1++)
            {
                f_IDj = node_friends(i, j1);

                sphfj[0] = node_pos_sph(f_IDj, 0); // Face 1 lat
                sphfj[1] = node_pos_sph(f_IDj, 1); // Face 1 lon

                arc = distanceBetweenSph2(sphfj, sphn);

                x1 = node_pos_map(i, j1 + 1, 0) * r_recip;
                y1 = node_pos_map(i, j1 + 1, 1) * r_recip;
                arc = sqrt(x1 * x1 + y1 * y1);

                matrix(j1 + 1, 0) = RBF(arc, eps);

                node_node_RBF(i, j1 + 1) = RBF(arc, eps);

                // loop through all faces again
                for (int j2 = 0; j2 < friend_num; j2++)
                {
                    // Get face ID
                    f_IDi = node_friends(i, j2);

                    sphfi[0] = node_pos_sph(f_IDi, 0); // Face 2 lat
                    sphfi[1] = node_pos_sph(f_IDi, 1); // Face 2 lon

                    // Get radial basis function between face centre j1 and j2
                    // using constant shapre parameter eps=2
                    arc = distanceBetweenSph2(sphfj, sphfi);

                    x2 = node_pos_map(i, j2 + 1, 0) * r_recip;
                    y2 = node_pos_map(i, j2 + 1, 1) * r_recip;
                    arc = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));

                    matrix(j1 + 1, j2 + 1) = RBF(arc, eps);
                }

                // arc = distanceBetweenSph2(sphn, sphfj);
                // node_node_RBF(i, j1+1) = RBF(arc, eps);
            }
            node_node_RBF(i, 0) = 1.0;

            // // Compute the Cholesky decomposition of the interpolation matrix
            // // interpSolvers[i].compute(matrix);

            // interpVandermondeScalar.push_back(matrix.completeOrthogonalDecomposition().pseudoInverse());
            interpVandermondeScalar.push_back(matrix);
            // std::cout<<matrix<<std::endl<<std::endl;
        }
    }

    // Get solvers for RBF matrix
    r = globals->radius.Value();
    {
        int friend_num, f_IDj, f_IDi;
        double sphn[2];
        double sphfi[2], sphfj[2];
        double njx, njy, njz;
        double nix, niy, niz;
        double arc, phi, aij;
        for (int i = 0; i < NODE_NUM; i++)
        {
            friend_num = node_fnum(i);

            DenseMat matrix(friend_num, friend_num);

            sphn[0] = node_pos_sph(i, 0); // Node lat
            sphn[1] = node_pos_sph(i, 1); // Node Lon

            // loop through all j faces
            for (int j1 = 0; j1 < friend_num; j1++)
            {
                f_IDj = faces(i, j1);

                sphfj[0] = face_centre_pos_sph(f_IDj, 0); // Face 1 lat
                sphfj[1] = face_centre_pos_sph(f_IDj, 1); // Face 1 lon

                njx = face_normal_vec_xyz(f_IDj, 0);
                njy = face_normal_vec_xyz(f_IDj, 1);
                njz = face_normal_vec_xyz(f_IDj, 2);

                // loop through all faces again

                for (int j2 = 0; j2 < friend_num; j2++)
                {
                    // Get face ID
                    f_IDi = faces(i, j2);

                    sphfi[0] = face_centre_pos_sph(f_IDi, 0); // Face 2 lat
                    sphfi[1] = face_centre_pos_sph(f_IDi, 1); // Face 2 lon

                    // normal vector to face j2
                    nix = face_normal_vec_xyz(f_IDi, 0);
                    niy = face_normal_vec_xyz(f_IDi, 1);
                    niz = face_normal_vec_xyz(f_IDi, 2);

                    // Get radial basis function between face centre j1 and j2
                    // using constant shapre parameter eps=2
                    arc = distanceBetweenSph2(sphfj, sphfi);
                    phi = RBF(arc, 1.0);

                    // coefficient is radial-basis-function (phi) x nf1 \cdot nf2
                    aij = phi * (njx * nix + njy * niy + njz * niz);

                    matrix(j1, j2) = aij;
                }

                arc = distanceBetweenSph2(sphn, sphfj);
                node_face_RBF(i, j1) = RBF(arc, 1.0);
            }

            // Compute the Cholesky decomposition of the interpolation matrix
            interpSolvers2[i].compute(matrix);

            interpVandermonde.push_back(matrix.inverse());
        }
    }

    return 1;
};

/*
// int Mesh::DefineBoundaryCells(void)
// {
//     // simply define if a control volume is either interior to a boundary (0),
//     // shares a face with a boundary (1), or is completely outside of the
//     // boundary (2).
//     int i, j, f, friend_num;
//     double lat1, lat2, lon1, lon2, dist;
//     double r = 1.0;

//     // r = globals->radius.Value();

//     for (i=0; i<NODE_NUM; i++)
//     {
//         f = node_friends(i,5);

//         friend_num = 6;                                     // Assume hexagon (6 centroids)
//         if (f == -1) {
//             friend_num = 5;                                 // Check if pentagon (5 centroids)
//             node_dists(i,5) = -1.0;
//         }

//         lat1 = 0.0;                     // set map coords for first
//         lon1 = pi;

//         lat2 = node_pos_sph(i, 0);                     // set map coords for second centroid
//         lon2 = node_pos_sph(i, 1);

//         distanceBetweenSph(dist, lat1, lat2, lon1, lon2, r);
//         // dist /= r;

//         // std::cout<<dist*180./pi<<std::endl;
//         if (dist*180./pi <= 90) cell_is_boundary(i) = 0;
//         else cell_is_boundary(i) = 1;

//         // if (cell_is_boundary(i) == 0)
//         // {
//         //   for (j=0; j<friend_num; j++)                        // Loop through all centroids in the control volume
//         //   {
//         //     int f_ID = faces(i,j);
//         //     lat2 = face_intercept_pos_sph(f_ID, 0);                     // set map coords for second centroid
//         //     lon2 = face_intercept_pos_sph(f_ID, 1);
//         //
//         //     distanceBetweenSph(dist, lat1, lat2, lon1, lon2, r);   //
//         //     // dist /= r;
//         //     // std::cout<<i<<'\t'<<dist*180./pi<<std::endl;
//         //     if (dist*180./pi < 20) {
//         //       // cell_is_boundary(i) = 1;
//         //       face_is_boundary(f_ID) = 1;
//         //       // break;
//         //     }
//         //   }
//         // }
//         // if (cell_is_boundary(i) == 2) std::cout<<i<<'\t'<<cell_is_boundary(i)<<std::endl;

//         cell_is_boundary(i) = 0;
//     }

// }*/

// Function to calculate the cosine(alpha) and sine(alpha) velocity tranform
// factors from Lee and Macdonald (2009). These factors are constant in time,
// and so are pre-calculated here to be easily accessed during the spatial
// integration.
int Mesh::CalcVelocityTransformFactors(void)
{
    int i, j, f;
    double *cos_a, *sin_a;
    double lat1, lat2, lon1, lon2;

    for (i = 0; i < NODE_NUM; i++)
    {
        lat1 = node_pos_sph(i, 0);
        lon1 = node_pos_sph(i, 1);

        // Set pointers to address of variables we want to change
        cos_a = &node_vel_trans(i, 0, 0);
        sin_a = &node_vel_trans(i, 0, 1);

        // Pass pointers by reference
        velTransform(*cos_a, *sin_a, lat1, lat1, lon1, lon1);

        // Loop through all friends (len is 7 as first element is parent node)
        // Assign node mapping factors and coordinates
        for (j = 1; j < 7; j++)
        {
            cos_a = &node_vel_trans(i, j, 0);
            sin_a = &node_vel_trans(i, j, 1);

            f = node_friends(i, j - 1); // j-1 because len node_friends[i] = 6;
            switch (f)
            {
            case -1: // if node is pentagon center
                *cos_a = -1.0;
                *sin_a = -1.0;
                break;
            default: // if node is a hexagon center
                lat2 = node_pos_sph(f, 0);
                lon2 = node_pos_sph(f, 1);

                // assign transform factors values to arrays
                velTransform(*cos_a, *sin_a, lat1, lat2, lon1, lon2);

                break;
            }
        }
    }

    for (i = 0; i < FACE_NUM; i++)
    {
        lat2 = face_intercept_pos_sph(i, 0);
        lon2 = face_intercept_pos_sph(i, 1);

        int n1 = face_nodes(i, 0);
        int n2 = face_nodes(i, 1);

        // lat1 = node_pos_sph(n1, 0);
        // lon1 = node_pos_sph(n1, 1);

        // double cos_a, sin_a;

        // // find which friend number face i is associated with for node n1
        // int face_friend_ID;
        // int friend_num = 6;
        // if ( node_friends(n1, 5) == -1 ) friend_num--;
        // for (j = 0; j<friend_num; j++) {
        //     if (node_friends(n1, j) == n2) {
        //         face_friend_ID = (j+friend_num-1)%friend_num;
        //         break;
        //     }
        // }

        // // Get the face normal direction in mapped coords with node n1 at
        // // the coordinate centre
        // double face_nx = control_vol_edge_normal_map(n1, face_friend_ID, 0);
        // double face_ny = control_vol_edge_normal_map(n1, face_friend_ID, 1);

        // // Make sure normal vector points in the positive direction
        // face_nx *= face_dir(i, 0);
        // face_ny *= face_dir(i, 0);

        // // Find angle between the mapping coordinate unit vectors and the
        // // face normal vector

        // // unit vector in x-direction = (1.0, 0.0)^T
        // double node_x_nx = 1.0;
        // double node_x_ny = 0.0;

        // double dot = face_nx*node_x_nx + face_ny*node_x_ny;      // dot product between [face_nx, face_ny] and [node_x_nx, node_x_ny]
        // double det = face_nx*node_x_ny - face_ny*node_x_nx ;     // determinant
        // double angle = atan2(det, dot);

        // cos_a = cos(angle);//face_nx*node_x_nx + face_ny*node_x_ny;
        // sin_a = sin(angle);//face_ny*node_x_nx - face_nx*node_x_ny;

        // face_node_vel_trans(i, 0, 0) = cos_a;
        // face_node_vel_trans(i, 0, 1) = sin_a;

        // friend_num = 6;
        // if ( node_friends(n2, 5) == -1 ) friend_num--;
        // for (j = 0; j<friend_num; j++) {
        //     if (node_friends(n2, j) == n1) {
        //         face_friend_ID = (j+friend_num-1)%friend_num;
        //         break;
        //     }
        // }

        // // Get the face normal direction in mapped coords with node n1 at
        // // the coordinate centre
        // face_nx = control_vol_edge_normal_map(n2, face_friend_ID, 0);
        // face_ny = control_vol_edge_normal_map(n2, face_friend_ID, 1);

        // // Make sure normal vector points in the positive direction
        // face_nx *= face_dir(i, 1);
        // face_ny *= face_dir(i, 1);

        // // Find angle between the mapping coordinate unit vectors and the
        // // face normal vector

        // dot = face_nx*node_x_nx + face_ny*node_x_ny;      // dot product between [face_nx, face_ny] and [node_x_nx, node_x_ny]
        // det = face_nx*node_x_ny - face_ny*node_x_nx;      // determinant
        // angle = atan2(det, dot);

        // cos_a = cos(angle);//face_nx*node_x_nx + face_ny*node_x_ny;
        // sin_a = sin(angle);

        // // cos_a = face_nx*node_x_nx + face_ny*node_x_ny;
        // // sin_a = face_ny*node_x_nx - face_nx*node_x_ny;

        // face_node_vel_trans(i, 1, 0) = cos_a;
        // face_node_vel_trans(i, 1, 1) = sin_a;

        // // std::cout<<cos_a<<' '<<sin_a<<' '<<face_ny<<std::endl;

        lat1 = face_intercept_pos_sph(i, 0);
        lon1 = face_intercept_pos_sph(i, 1);

        lat2 = node_pos_sph(n1, 0);
        lon2 = node_pos_sph(n1, 1);

        // Set pointers to address of variables we want to change
        cos_a = &face_node_vel_trans(i, 0, 0);
        sin_a = &face_node_vel_trans(i, 0, 1);

        // Pass pointers by reference
        velTransform(*cos_a, *sin_a, lat1, lat2, lon1, lon2);

        lat2 = node_pos_sph(n2, 0);
        lon2 = node_pos_sph(n2, 1);

        // Set pointers to address of variables we want to change
        cos_a = &face_node_vel_trans(i, 1, 0);
        sin_a = &face_node_vel_trans(i, 1, 1);

        // Pass pointers by reference
        velTransform(*cos_a, *sin_a, lat1, lat2, lon1, lon2);
    }
    return 1;
};

// Function to convert all spherical coorinate quantities to mapping coordinates
// x and y from Lee and Macdonald (2009). Mapping factors are also calculated
// and stored. These quantities are stored for node positions, centroid positions,
// node mapping factors, and centroid mapping factors.
int Mesh::CalcMappingCoords(void)
{
    int i, j, f;
    double *x, *y, *m, r;
    double lat1, lat2, lon1, lon2;

    m = new double;

    r = globals->radius.Value();

    for (i = 0; i < NODE_NUM; i++)
    {
        lat1 = node_pos_sph(i, 0);
        lon1 = node_pos_sph(i, 1);

        // Set pointers to address of variables we want to change
        x = &node_pos_map(i, 0, 0);
        y = &node_pos_map(i, 0, 1);

        // Pass pointers by reference
        mapAtPoint(*m, *x, *y, lat1, lat1, lon1, lon1, r);

        // Loop through all friends (len is 7 as first element is parent node)
        // Assign node mapping factors and coordinates
        for (j = 1; j < 7; j++)
        {
            x = &node_pos_map(i, j, 0);
            y = &node_pos_map(i, j, 1);

            f = node_friends(i, j - 1); // j-1 because len node_friends[i] = 6;
            switch (f)
            {
            case -1: // if node is pentagon center
                *m = -1.0;
                *x = -1.0;
                *y = -1.0;
                break;
            default: // if node is a hexagon center
                lat2 = node_pos_sph(f, 0);
                lon2 = node_pos_sph(f, 1);

                // assign mapped values to arrays
                mapAtPoint(*m, *x, *y, lat1, lat2, lon1, lon2, r);

                break;
            }
        }
    }

    for (i = 0; i < NODE_NUM; i++)
    {
        lat1 = node_pos_sph(i, 0);
        lon1 = node_pos_sph(i, 1);

        // Assign centroid mapping factors and coordinates
        for (j = 0; j < 6; j++)
        {
            x = &centroid_pos_map(i, j, 0);
            y = &centroid_pos_map(i, j, 1);

            f = node_friends(i, j);
            switch (f)
            {
            case -1: // if node is pentagon center
                *m = -1.0;
                *x = -1.0;
                *y = -1.0;
                break;
            default: // if node is a hexagon center
                lat2 = centroid_pos_sph(i, j, 0);
                lon2 = centroid_pos_sph(i, j, 1);

                // assign mapped values to arrays
                mapAtPoint(*m, *x, *y, lat1, lat2, lon1, lon2, r);
                break;
            }
        }
    }

    delete m;

    return 1;
};

/*
// Function to find the length of control volume edges for each node, in mapping
// coordinates. Values are stored in control_vol_edge_len 2D array.
// int Mesh::CalcControlVolumeEdgeLengths(void)
// {
//     int i, j, f, friend_num;
//     double * x1, * y1, * x2, * y2, * edge_len;

//     for (i=0; i<NODE_NUM; i++)
//     {

//         f = node_friends(i,5);

//         friend_num = 6;                                     // Assume hexagon (6 centroids)
//         if (f == -1) {
//             friend_num = 5;                                   // Check if pentagon (5 centroids)
//             control_vol_edge_len(i,5) = -1.0;
//         }

//         for (j=0; j<friend_num; j++)                        // Loop through all centroids in the control volume
//         {
//             edge_len = &control_vol_edge_len(i,j);            // set pointer to edge length array

//             x1 = &centroid_pos_map(i, j, 0);                  // set map coords for first centroid
//             y1 = &centroid_pos_map(i, j, 1);

//             x2 = &centroid_pos_map(i, (j+1)%friend_num, 0);   // set map coords for second centroid
//             y2 = &centroid_pos_map(i, (j+1)%friend_num, 1);   // automatically loops around using %

//             distanceBetween(*edge_len, *x1, *x2, *y1, *y2);   // calculate distance between the two centroids.
//             // Edge_len automatically assigned the length
//         }
//     }

//     return 1;
// };
*/


int Mesh::CalcMaxTimeStep(void)
{
    double dt, dt_target;
    int dt_num = 100;

    dt_target = globals->timeStep.Value();

    dt = globals->period.Value();

    std::cout << "TARGET DT: " << dt_target << std::endl;
    while (dt > dt_target)
    {
        dt = globals->period.Value() / (double)dt_num;

        dt_num += 100;

        std::cout << dt << ' ' << dt_num << std::endl;
    }
    dt_num -= 100;

    globals->timeStep.SetValue(dt);
    globals->totalIter.SetValue(dt_num);

    return 1;
}

// Function to evaluate trigonometric functions over latitude and longitude for
// every node.
int Mesh::CalcTrigFunctions(void)
{
    int i;
    double lat, lon;

    for (i = 0; i < NODE_NUM; i++)
    {
        lat = node_pos_sph(i, 0);
        lon = node_pos_sph(i, 1);

        trigLat(i, 0) = cos(lat);
        trigLat(i, 1) = sin(lat);

        trigLon(i, 0) = cos(lon);
        trigLon(i, 1) = sin(lon);

        trig2Lat(i, 0) = cos(2.0 * lat);
        trig2Lat(i, 1) = sin(2.0 * lat);

        trig2Lon(i, 0) = cos(2.0 * lon);
        trig2Lon(i, 1) = sin(2.0 * lon);

        trigSqLat(i, 0) = pow(cos(lat), 2.0);
        trigSqLat(i, 1) = pow(sin(lat), 2.0);

        trigSqLon(i, 0) = pow(cos(lon), 2.0);
        trigSqLon(i, 1) = pow(sin(lon), 2.0);
    }

    for (int i = 0; i < VERTEX_NUM; i++)
    {
        vertex_sinlat(i) = sin(vertex_pos_sph(i, 0));
    }

    return 1;
};

/*int Mesh::CalcLegendreFuncs(void)
{
    // int i, l, m, l_max;
    // double cosCoLat;

    // Array2D<double> * temp_legendre;    // temp array for legendre polynomials
    // Array2D<double> * temp_dlegendre;   // temp array for legendre polynomial derivs

    // //-----------------------------

    // l_max = globals->l_max.Value();

    // temp_legendre = new Array2D<double>(2 * (l_max+1), 2 * (l_max+1));
    // temp_dlegendre = new Array2D<double>(2 * (l_max+1), 2 * (l_max+1));

    // Array3D<double> Pbar_lm(NODE_NUM, l_max+1, l_max+1);
    // Array3D<double> Pbar_cosMLon(NODE_NUM, l_max+1, l_max+1);
    // Array3D<double> Pbar_sinMLon(NODE_NUM, l_max+1, l_max+1);

    // for (i=0; i<NODE_NUM; ++i)
    // {
    //     cosCoLat = cos(pi*0.5 - node_pos_sph(i,0));     // cos(colatitude) of node i

    //     getLegendreFuncs(cosCoLat, *temp_legendre, l_max);
    //     getLegendreFuncsDeriv(cosCoLat, *temp_dlegendre, l_max);



    //     for (l = 0; l < l_max+1; l++)
    //     {
    //         for (m = 0; m <= l; m++)
    //         {

    //             // assign legendre pol at node i for degree l and order m
    //             Pbar_lm(i, l, m) = (*temp_legendre)(l, m);

    //             // assign legendre pol derivative at node i for degree l and order m
    //             // Pbar_lm_deriv(i, l, m) = (*temp_dlegendre)(l, m)*sin(pi*0.5 - node_pos_sph(i, 0));

    //             if ((l==2) && (m==0)) Pbar_20(i) = (*temp_legendre)(l, m);
    //             if ((l==2) && (m==2)) Pbar_22(i) = (*temp_legendre)(l, m);

    //         }
    //     }

    //     // calculate cos(m*longitude) and sin(m*longitude) for node i
    //     for (m=0; m < l_max+1; m++)
    //     {
    //         trigMLon(i, m, 0) = cos( (double)m * node_pos_sph( i, 1 ) );
    //         trigMLon(i, m, 1) = sin( (double)m * node_pos_sph( i, 1 ) );
    //     }

    //     // int count = 0;
    //     for (l = 0; l < l_max+1; l++)
    //     {
    //         for (m = 0; m <= l; m++)
    //         {
    //             Pbar_cosMLon(i, l, m) = Pbar_lm(i, l, m)*trigMLon(i, m, 0);
    //             Pbar_sinMLon(i, l, m) = Pbar_lm(i, l, m)*trigMLon(i, m, 1);

    //             // Pbar_deriv_cosMLon(i, l, m) = Pbar_lm_deriv(i, l, m)*trigMLon(i, m, 0);
    //             // Pbar_deriv_sinMLon(i, l, m) = Pbar_lm_deriv(i, l, m)*trigMLon(i, m, 1);

    //         }
    //     }

    //     int count = 0;
    //     for (l = 2; l < l_max+1; l++)
    //     {
    //         for (m = 0; m <= l; m++)
    //         {
    //             sh_matrix(i, count) = Pbar_cosMLon(i, l, m);
    //             if (globals->surface_type == FREE_LOADING) sh_matrix(i, count) *= globals->loading_factor[l];
    //             else if (globals->surface_type == LID_MEMBR ||
    //                      globals->surface_type == LID_LOVE) sh_matrix(i, count) *= globals->shell_factor_beta[l];


    //             count++;

    //             sh_matrix(i, count) = Pbar_sinMLon(i, l, m);
    //             if (globals->surface_type == FREE_LOADING) sh_matrix(i, count) *= globals->loading_factor[l];
    //             else if (globals->surface_type == LID_MEMBR ||
    //                      globals->surface_type == LID_LOVE) sh_matrix(i, count) *= globals->shell_factor_beta[l];

    //             count++;
    //         }
    //     }
    // }

    // int count = 0;
    // for (m = 0; m < (l_max+1)*(l_max+2) - 6; m++)
    // {
    //     for (i=0; i<NODE_NUM; i++)
    //     {
    //         sh_matrix_fort(count) = sh_matrix(i, m);
    //         // std::cout<<sh_matrix_fort(count)<<std::endl;
    //         count++;
    //     }
    // }

    // // free up memory!
    // delete temp_legendre;
    // delete temp_dlegendre;

    // return 1;
}
*/

// Function to create the adjaceny matrix that maps the neighbouring faces to each node
int Mesh::CalcRBFInterpMatrix(void)
{
    //--------------- Construct sparse matrix in CSR format --------------------

    //--------------- Define objects to construct matrix -----------------------

    // There are 12 pentagons on the geodesic grid, each with 5 faces.
    // The other cells are hexagons with 6
    // faces. Therefore, the total number
    // of non-zero coefficients in our sparse matrix is given as:
    int nNonZero = 5 * 5 * 12 + (NODE_NUM - 12) * 6 * 6;

    SpMat Vinv = SpMat(12 * 5 + (NODE_NUM - 12) * 6, 12 * 5 + (NODE_NUM - 12) * 6);

    // Allocate memory for non-zero entries
    Vinv.reserve(nNonZero);

    std::vector<Triplet> coefficients;
    coefficients.reserve(nNonZero);

    // assign the non-zero coefficients and their column indexes in CSR format
    int friend_num, face_ID;
    double coeff;
    int count = 0;
    int nh = 0;
    int np = 0;

    for (int i = 0; i < NODE_NUM; ++i)
    {
        friend_num = node_fnum(i);

        for (int x = 0; x < friend_num; x++)
        {
            int row = x + nh * 6 + np * 5;
            for (int y = 0; y < friend_num; y++)
            {

                int col = y + nh * 6 + np * 5;

                coeff = interpVandermonde[i].coeff(x, y);

                coefficients.push_back(Triplet(row, col, interpVandermonde[i].coeff(x, y)));
            }
            count++;
        }

        if (friend_num == 5)
            np++;
        else
            nh++;
    }

    // Fill the sparse operator with the list of triplets
    Vinv.setFromTriplets(coefficients.begin(), coefficients.end());
    Vinv.makeCompressed();

    nNonZero = 3 * (12 * 5 + (NODE_NUM - 12) * 6);
    SpMat RBF = SpMat(3 * NODE_NUM, 12 * 5 + (NODE_NUM - 12) * 6);
    RBF.reserve(nNonZero);

    std::vector<Triplet> coefficients2;
    coefficients2.reserve(nNonZero);

    count = 0;
    // int face_ID;
    for (int i = 0; i < NODE_NUM; ++i)
    {
        friend_num = node_fnum(i);

        for (int j = 0; j < friend_num; j++)
        {
            face_ID = faces(i, j);

            coeff = node_face_RBF(i, j) * face_normal_vec_xyz(face_ID, 0); // x component row
            coefficients2.push_back(Triplet(3 * i, count, coeff));

            coeff = node_face_RBF(i, j) * face_normal_vec_xyz(face_ID, 1); // y component row
            coefficients2.push_back(Triplet(3 * i + 1, count, coeff));

            coeff = node_face_RBF(i, j) * face_normal_vec_xyz(face_ID, 2); // z component row
            coefficients2.push_back(Triplet(3 * i + 2, count, coeff));
            count++;
        }
    }

    RBF.setFromTriplets(coefficients2.begin(), coefficients2.end());
    RBF.makeCompressed();

    operatorRBFinterp = RBF * (Vinv * node2faceAdj);
    operatorRBFinterp.makeCompressed();

    return 1;
}

// Function to create the adjaceny matrix that maps the neighbouring faces to each node
int Mesh::CalcRBFInterpMatrix2(void)
{
    //--------------- Define objects to construct matrix -----------------------

    // There are 12 pentagons on the geodesic grid, each with 5 faces.
    // The other cells are hexagons with 6
    // faces. Therefore, the total number
    // of non-zero coefficients in our sparse matrix is given as:
    int nNonZero = 6 * 6 * 12 + (NODE_NUM - 12) * 7 * 7;

    vandermondeInv = SpMat(12 * 6 + (NODE_NUM - 12) * 7, 12 * 6 + (NODE_NUM - 12) * 7);
    // vandermondeInv = SpMat(12*6 + (NODE_NUM - 12)*7, NODE_NUM);

    // Allocate memory for non-zero entries
    vandermondeInv.reserve(nNonZero);

    std::vector<Triplet> coefficients;
    coefficients.reserve(nNonZero);

    // assign the non-zero coefficients and their column indexes in CSR format
    int friend_num;
    double coeff;
    // int count=0;
    int nh = 0;
    int np = 0;

    for (int i = 0; i < NODE_NUM; ++i)
    {
        friend_num = node_fnum(i)+1;

        Eigen::MatrixXd mat = interpVandermondeScalar[i].inverse();

        for (int x = 0; x < friend_num; x++)
        {
            int row = x + nh * 7 + np * 6;
            for (int y = 0; y < friend_num; y++)
            {

                int col = y + nh * 7 + np * 6;

                coeff = mat.coeff(x, y);

                coefficients.push_back(Triplet(row, col, coeff));
            }
        }

        if (friend_num == 6)
            np++;
        else
            nh++;
    }

    // Fill the sparse operator with the list of triplets
    vandermondeInv.setFromTriplets(coefficients.begin(), coefficients.end());
    vandermondeInv.makeCompressed();

    nNonZero = (6 * 12 + (NODE_NUM - 12) * 7) * 3;

    SpMat rbfDeriv2(3 * NODE_NUM, 12 * 6 + (NODE_NUM - 12) * 7);
    // vandermondeInv = SpMat(12*6 + (NODE_NUM - 12)*7, NODE_NUM);

    // Allocate memory for non-zero entries
    rbfDeriv2.reserve(nNonZero);

    std::vector<Triplet> coefficients2;
    coefficients2.reserve(nNonZero);

    {
        double arc;
        double eps = 0.125; // If you change this, make sure you also change it in CalcControlVolumeInterpMatrix!!!!
        double x1, y1;
        double xx, yy, xy, arc2, arc3, dphidr, d2phidr2;
        int row, col;
        double r_recip = 1.0 / globals->radius.Value();
        // RBFSUM = DenseMat(3*NODE_NUM, 12*6 + (NODE_NUM - 12)*7);
        np = 0;
        nh = 0;
        for (int i = 0; i < NODE_NUM; i++)
        {
            // check if node is pentagon
            friend_num = node_fnum(i);

            coeff = -2.0 * eps * eps;

            row = 3 * i;
            col = (6 * np + 7 * nh);
            coefficients2.push_back(Triplet(row, col, coeff));

            row = 3 * i + 1;
            // coefficients2.push_back( Triplet(row, col, 0.0));

            row = 3 * i + 2;
            coefficients2.push_back(Triplet(row, col, coeff));

            // loop through all j faces
            for (int j = 0; j < friend_num; j++)
            {
                x1 = node_pos_map(i, j + 1, 0) * r_recip;
                y1 = node_pos_map(i, j + 1, 1) * r_recip;
                arc = sqrt(x1 * x1 + y1 * y1);

                xx = x1 * x1;
                yy = y1 * y1;
                xy = x1 * y1;
                arc2 = 1.0 / (arc * arc);
                arc3 = 1.0 / (arc * arc * arc);

                double rbf = node_node_RBF(i, j + 1);
                d2phidr2 = 2 * eps * eps * rbf * (2 * arc * arc * eps * eps - 1);
                dphidr = -2 * arc * eps * eps * rbf;

                row = 3 * i;
                col = (6 * np + 7 * nh) + j + 1;
                coeff = (yy * arc3 * dphidr + xx * arc2 * d2phidr2);
                coefficients2.push_back(Triplet(row, col, coeff));

                row = 3 * i + 1;
                coeff = (-xy * arc3 * dphidr + xy * arc2 * d2phidr2);
                coefficients2.push_back(Triplet(row, col, coeff));

                row = 3 * i + 2;
                coeff = (xx * arc3 * dphidr + yy * arc2 * d2phidr2);
                coefficients2.push_back(Triplet(row, col, coeff));
            }

            if (friend_num == 6)
                nh++;
            else
                np++;
        }

        rbfDeriv2.setFromTriplets(coefficients2.begin(), coefficients2.end());
        rbfDeriv2.makeCompressed();

        // std::cout<<rbfDeriv2.nonZeros()<<std::endl;
        operatorSecondDeriv = r_recip * r_recip * rbfDeriv2 * vandermondeInv * node2nodeAdj;
    }

    nNonZero = 2 * 9 * FACE_NUM;

    SpMat R(3 * 2 * FACE_NUM, 3 * NODE_NUM);

    // Allocate memory for non-zero entries
    R.reserve(nNonZero);

    std::vector<Triplet> coefficients3;
    coefficients3.reserve(nNonZero);

    nNonZero = 2 * 3 * FACE_NUM;

    // SpMat operatorDirectionalSecondDeriv(2*FACE_NUM, 3*NODE_NUM);

    SpMat N(2 * FACE_NUM, 2 * 3 * FACE_NUM);

    // Allocate memory for non-zero entries
    N.reserve(nNonZero);

    std::vector<Triplet> coefficients4;
    coefficients4.reserve(nNonZero);

    {
        int inner_ID;
        int outer_ID;
        double nx, ny;
        double cosa, sina;
        int row, col;
        for (int i = 0; i < FACE_NUM; i++)
        {
            inner_ID = face_nodes(i, 0);

            cosa = face_node_vel_trans(i, 0, 0);
            sina = face_node_vel_trans(i, 0, 1);

            // fxx = 2*(c3*cosa*cosa - c4*cosa*sina + c5*sina*sina);
            // fyy = 2*(c3*sina*sina + c4*cosa*sina + c5*cosa*cosa);
            // fxy = 2*(c3 - c5)*cosa*sina + c4*(cosa*cosa - sina*sina);

            row = 6 * i;
            col = 3 * inner_ID;
            coeff = cosa * cosa;
            coefficients3.push_back(Triplet(row, col, coeff));

            col = 3 * inner_ID + 1;
            coeff = -2 * cosa * sina;
            coefficients3.push_back(Triplet(row, col, coeff));

            col = 3 * inner_ID + 2;
            coeff = sina * sina;
            coefficients3.push_back(Triplet(row, col, coeff));

            row = 6 * i + 1;
            col = 3 * inner_ID;
            coeff = cosa * sina;
            coefficients3.push_back(Triplet(row, col, coeff));

            col = 3 * inner_ID + 1;
            coeff = cosa * cosa - sina * sina;
            coefficients3.push_back(Triplet(row, col, coeff));

            col = 3 * inner_ID + 2;
            coeff = -sina * cosa;
            coefficients3.push_back(Triplet(row, col, coeff));

            row = 6 * i + 2;
            col = 3 * inner_ID;
            coeff = sina * sina;
            coefficients3.push_back(Triplet(row, col, coeff));

            col = 3 * inner_ID + 1;
            coeff = 2 * cosa * sina;
            coefficients3.push_back(Triplet(row, col, coeff));

            col = 3 * inner_ID + 2;
            coeff = cosa * cosa;
            coefficients3.push_back(Triplet(row, col, coeff));

            outer_ID = face_nodes(i, 1);

            cosa = face_node_vel_trans(i, 1, 0);
            sina = face_node_vel_trans(i, 1, 1);

            row = 6 * i + 3;
            col = 3 * outer_ID;
            coeff = cosa * cosa;
            coefficients3.push_back(Triplet(row, col, coeff));

            col = 3 * outer_ID + 1;
            coeff = -2 * cosa * sina;
            coefficients3.push_back(Triplet(row, col, coeff));

            col = 3 * outer_ID + 2;
            coeff = sina * sina;
            coefficients3.push_back(Triplet(row, col, coeff));

            row = 6 * i + 4;
            col = 3 * outer_ID;
            coeff = cosa * sina;
            coefficients3.push_back(Triplet(row, col, coeff));

            col = 3 * outer_ID + 1;
            coeff = cosa * cosa - sina * sina;
            coefficients3.push_back(Triplet(row, col, coeff));

            col = 3 * outer_ID + 2;
            coeff = -sina * cosa;
            coefficients3.push_back(Triplet(row, col, coeff));

            row = 6 * i + 5;
            col = 3 * outer_ID;
            coeff = sina * sina;
            coefficients3.push_back(Triplet(row, col, coeff));

            col = 3 * outer_ID + 1;
            coeff = 2 * cosa * sina;
            coefficients3.push_back(Triplet(row, col, coeff));

            col = 3 * outer_ID + 2;
            coeff = cosa * cosa;
            coefficients3.push_back(Triplet(row, col, coeff));

            nx = face_normal_vec_map(i, 0);
            ny = face_normal_vec_map(i, 1);

            row = 2 * i;
            coefficients4.push_back(Triplet(row, 6 * i, nx * nx));
            coefficients4.push_back(Triplet(row, 6 * i + 1, 2 * nx * ny));
            coefficients4.push_back(Triplet(row, 6 * i + 2, ny * ny));

            row = 2 * i + 1;
            coefficients4.push_back(Triplet(row, 6 * i + 3, nx * nx));
            coefficients4.push_back(Triplet(row, 6 * i + 4, 2 * nx * ny));
            coefficients4.push_back(Triplet(row, 6 * i + 5, ny * ny));
        }
    }

    R.setFromTriplets(coefficients3.begin(), coefficients3.end());
    N.setFromTriplets(coefficients4.begin(), coefficients4.end());

    R.makeCompressed();
    N.makeCompressed();

    // operatorDirectionalSecondDeriv = SpMat(2*FACE_NUM, 3*NODE_NUM);

    // operatorDirectionalSecondDeriv = N*R;
    // operatorDirectionalSecondDeriv.makeCompressed();

    operatorDirectionalSecondDeriv = SpMat(2 * FACE_NUM, NODE_NUM);

    operatorDirectionalSecondDeriv = N * R * operatorSecondDeriv;
    operatorDirectionalSecondDeriv.makeCompressed();
    operatorDirectionalSecondDeriv.prune(0.0);

    return 1;
}

// Function to create the adjaceny matrix that maps the neighbouring faces to each node
int Mesh::CalcAdjacencyMatrix(void)
{
    int count = 0, friend_num, face_ID, node_ID;
    double coeff = 1.0; // matrix coeff is always 1

    // There are 12 pentagons on the geodesic grid, each with 5 faces.
    // The other cells are hexagons with 6
    // faces. Therefore, the total number
    // of non-zero coefficients in our sparse matrix is given as:
    int nNonZero = 12 * 5 + (NODE_NUM - 12) * 6;

    node2faceAdj = SpMat(nNonZero, FACE_NUM);

    // Allocate memory for non-zero entries
    node2faceAdj.reserve(nNonZero);

    std::vector<Triplet> coefficients;
    coefficients.reserve(nNonZero);


    count = 0;
    for (unsigned i = 0; i < NODE_NUM; ++i)
    {
        friend_num = node_fnum(i);

        for (unsigned j = 0; j < friend_num; j++)
        {
            face_ID = faces(i, j);

            coefficients.push_back(Triplet(count, face_ID, coeff));
            count++;
        }
    }


    // Fill the sparse operator with the list of triplets
    node2faceAdj.setFromTriplets(coefficients.begin(), coefficients.end());
    node2faceAdj.makeCompressed();

    // 6 and 7 here to include the central node itself
    nNonZero = 12 * 6 + (NODE_NUM - 12) * 7;

    node2nodeAdj = SpMat(nNonZero, NODE_NUM);

    // Allocate memory for non-zero entries
    node2nodeAdj.reserve(nNonZero);

    std::vector<Triplet> coefficients2;
    coefficients2.reserve(nNonZero);

    // assign the non-zero coefficients and their column indexes in CSR format
    
    count = 0;
    for (unsigned i = 0; i < NODE_NUM; ++i)
    {
        friend_num = node_fnum(i);

        coefficients2.push_back(Triplet(count, i, coeff)); // Add itself

        count++;
        for (unsigned j = 0; j < friend_num; j++)
        {
            node_ID = node_friends(i, j);

            coefficients2.push_back(Triplet(count, node_ID, coeff));

            count++;
        }
    }
    

    // Fill the sparse operator with the list of triplets
    node2nodeAdj.setFromTriplets(coefficients2.begin(), coefficients2.end());
    node2nodeAdj.makeCompressed();

    return 1;
}

int Mesh::CalcCoriolisOperatorCoeffs(void)
{
    int fnum[2];
    int f_ID;
    double coeff;
    int n1, n2;

    //--------------- Define objects to construct matrix -----------------------

    // Here we assign the number of non-zero matrix elements for each operator.
    // There are 12 pentagons on the geodesic grid, each with 5 neighbours and
    // one central node (total of 6). The other cells are hexagons with 6
    // neighbours and one central node (total of 7). Additionally, there is a
    // row for each variable (total 3). Therefore, the total number
    // of non-zero coefficients in our sparse matrix is given as:
    int nNonZero = (12 * 5) * 9 + (FACE_NUM - 12 * 5) * 10;

    operatorCoriolis = SpMat(FACE_NUM, FACE_NUM);

    // Allocate memory for non-zero entries
    operatorCoriolis.reserve(nNonZero);

    std::vector<Triplet> coefficients;
    coefficients.reserve(nNonZero);

    for (unsigned i = 0; i < FACE_NUM; ++i)
    {
        n1 = face_nodes(i, 0);
        n2 = face_nodes(i, 1);

        fnum[0] = node_fnum(n1);
        fnum[1] = node_fnum(n2);

        // Loop through the faces adjacent to face i and use
        // interpolation weights to find the tangential velocity component
        for (unsigned k=0; k<2; k++) {
            for (unsigned j = 0; j < fnum[k]-1; j++)
            {
                f_ID = face_friends(i, k, j);
                coeff = -2.0 * globals->angVel.Value() * sin(face_centre_pos_sph(i, 0)) * face_interp_weights(i, k, j) * face_len(f_ID) * face_node_dist_r(i);
            
                coefficients.push_back(Triplet(i, f_ID, coeff));
            }
        }
    }

    // Fill the sparse operator with the list of triplets
    operatorCoriolis.setFromTriplets(coefficients.begin(), coefficients.end());
    operatorCoriolis.makeCompressed();

    return 1;
}

int Mesh::CalcLinearDragOperatorCoeffs(void)
{
    //--------------- Define objects to construct matrix -----------------------

    // Here we assign the number of non-zero matrix elements for each operator.
    // There are 12 pentagons on the geodesic grid, each with 5 neighbours and
    // one central node (total of 6). The other cells are hexagons with 6
    // neighbours and one central node (total of 7). Additionally, there is a
    // row for each variable (total 3). Therefore, the total number
    // of non-zero coefficients in our sparse matrix is given as:
    int nNonZero = FACE_NUM;

    // Define shape of operator
    operatorLinearDrag = SpMat(FACE_NUM, FACE_NUM);

    // Allocate memory for non-zero entries
    operatorLinearDrag.reserve(nNonZero);

    for (unsigned i = 0; i < FACE_NUM; i++)
    {
        operatorLinearDrag.insert(i, i) = -globals->alpha.Value();
    }

    operatorLinearDrag.makeCompressed();

    return 1;
}

int Mesh::CalcGradOperatorCoeffs(void)
{
    // Here we assign the number of non-zero matrix elements for each operator.
    // There are 12 pentagons on the geodesic grid, each with 5 neighbours and
    // one central node (total of 6). The other cells are hexagons with 6
    // neighbours and one central node (total of 7). Additionally, there is a
    // row for each variable (total 3). Therefore, the total number
    // of non-zero coefficients in our sparse matrix is given as:
    int nNonZero = 2 * FACE_NUM; //(NODE_NUM-12)*7 + 12*6;

    operatorGradient = SpMat(FACE_NUM, NODE_NUM);

    // Allocate memory for non-zero entries
    operatorGradient.reserve(nNonZero);

    std::vector<Triplet> coefficients;
    coefficients.reserve(nNonZero);

    double coeff;
    double dist_r;
    int inner_node_ID, outer_node_ID;
    for (int i = 0; i < FACE_NUM; i++)
    {
        inner_node_ID = face_nodes(i, 0);
        outer_node_ID = face_nodes(i, 1);

        dist_r = face_node_dist_r(i);

        // coeff = (-face_centre_m(i, 0)) / dist;
        // coefficients.push_back(Triplet(i, inner_node_ID, coeff));

        // coeff = (face_centre_m(i, 1)) / dist;
        // coefficients.push_back(Triplet(i, outer_node_ID, coeff));

        // Removed face_centre_m from grad calculation as now 
        // face_node_dist is the actual arc length (in metres)
        coeff = -dist_r;
        coefficients.push_back(Triplet(i, inner_node_ID, coeff));

        coeff = dist_r;
        coefficients.push_back(Triplet(i, outer_node_ID, coeff));
    }

    operatorGradient.setFromTriplets(coefficients.begin(), coefficients.end());
    operatorGradient.makeCompressed();

    return 1;
}

int Mesh::CalcCurlOperatorCoeffs(void)
{
    double area_r;
    int face_id;
    double node_dist;
    int t_ev;
    double coeff;

    //--------------- Define objects to construct matrix -----------------------

    // Here we assign the number of non-zero matrix elements for each operator.
    // The curl operator is computed at each vertex using the velocity normal to
    // each adjoining face. Each vertex has three faces, so the total number of
    // non-zero coefficients is:
    int nNonZeroCoeffs = VERTEX_NUM * 3;

    operatorCurl = SpMat(VERTEX_NUM, FACE_NUM);

    // Allocate memory for non-zero entries
    operatorCurl.reserve(nNonZeroCoeffs);

    std::vector<Triplet> coefficients;
    coefficients.reserve(nNonZeroCoeffs);

    // assign the non-zero coefficients
    for (unsigned i = 0; i < VERTEX_NUM; i++)
    {
        area_r = vertex_area_r(i);

        for (unsigned j = 0; j < 3; j++)
        {
            face_id = vertex_faces(i, j);

            // get distance between nodes crossing face
            // face_id
            node_dist = face_node_dist(face_id);

            // get indicator function (does a face velocity
            // contribute positively or negatively to the
            // curl at a vertex?)
            t_ev = vertex_face_dir(i, j);

            coeff = node_dist * (double)t_ev * area_r;

            coefficients.push_back(Triplet(i, face_id, coeff));
        }
    }

    operatorCurl.setFromTriplets(coefficients.begin(), coefficients.end());
    operatorCurl.makeCompressed();

    return 1;
}

int Mesh::CalcDivOperatorCoeffs(void)
{
    int f_num, face_id, dir;
    double area_r, edge_len, coeff;

    //--------------- Define objects to construct matrix -----------------------

    // Here we assign the number of non-zero matrix elements for each operator.
    // There are 12 pentagons on the geodesic grid, each with 5 faces.
    // The other cells are hexagons with 6 faces. Therefore, the total
    // number of non-zero coefficients in our sparse matrix is given as:
    int nNonZero = ((NODE_NUM - 12) * 6 + 12 * 5);

    operatorDivergence = SpMat(NODE_NUM, FACE_NUM);

    // Allocate memory for non-zero entries
    operatorDivergence.reserve(nNonZero);

    std::vector<Triplet> coefficients;
    coefficients.reserve(nNonZero);

    for (unsigned i = 0; i < NODE_NUM; i++)
    {
        f_num = node_fnum(i);

        for (unsigned j = 0; j < f_num; j++)
        {
            face_id = faces(i, j);

            // get edge length of current edge
            edge_len = face_len(face_id);
            dir = node_face_dir(i, j);
            area_r = cv_area_sph_r(i);

            coeff = -dir * edge_len * area_r;
            coefficients.push_back(Triplet(i, face_id, coeff));
        }
    }

    operatorDivergence.setFromTriplets(coefficients.begin(), coefficients.end());
    operatorDivergence.makeCompressed();

    return 1;
}

int Mesh::CalcFreeSurfaceSolver(void)
{
    double g, h, dt;
    g = globals->g.Value();
    h = globals->h.Value();
    dt = globals->timeStep.Value();

    // DenseMat I = DenseMat::Identity(NODE_NUM);
    auto I = DenseMat::Identity(NODE_NUM, NODE_NUM);

    SpMat Laplacian = -0.5 * 0.5 * g * h * dt * dt * (-operatorDivergence * operatorGradient);

    Laplacian += I.sparseView();

    cg.compute(Laplacian);
    return 1;
}

int Mesh::ReadGridFile()
{
    double r = globals->radius.Value();
    double rr = r*r;

    int count;

    std::string file_str = globals->path + SEP + "input_files" + SEP + "grid_l" + std::to_string(GRID_LVL) + ".h5";

    hid_t file_id = H5Fopen(file_str.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    //-----------------------------------------------------------------------------------------
    //------------------------------ LOAD NODE INFO -------------------------------------------
    //-----------------------------------------------------------------------------------------

    // LATITUDE OF NODE
    unsigned n_lat_size=0;
    double * n_lat_h5 = ReadHDF5Dataset<double>(&file_id, "/NODES/LAT/", H5T_NATIVE_DOUBLE, n_lat_size);

    // LONGITUDE OF NODE
    unsigned n_lon_size=0;
    double * n_lon_h5 = ReadHDF5Dataset<double>(&file_id, "/NODES/LON/", H5T_NATIVE_DOUBLE, n_lon_size);

    // AREA OF NODE
    unsigned n_area_size=0;
    double * n_area_h5 = ReadHDF5Dataset<double>(&file_id, "/NODES/AREA/", H5T_NATIVE_DOUBLE, n_area_size);

    // NUMBER OF NODES SURROUNDING NODE
    unsigned n_n_num_size=0;
    unsigned * n_n_num_h5 = ReadHDF5Dataset<unsigned>(&file_id, "/NODES/FRIENDS/NODES/FRIEND_NUM/", H5T_NATIVE_UINT, n_n_num_size);

    // ID'S OF NODES SURROUNDING NODE
    unsigned n_n_IDs_size=0;
    unsigned * n_n_IDs_h5 = ReadHDF5Dataset<unsigned>(&file_id, "/NODES/FRIENDS/NODES/ID/", H5T_NATIVE_UINT, n_n_IDs_size);

    // ARC DISTANCES BETWEEN NODE AND ITS SURROUNDING NODES
    unsigned n_n_arc_size=0;
    double * n_n_arc_h5 = ReadHDF5Dataset<double>(&file_id, "/NODES/FRIENDS/NODES/DISTANCE/", H5T_NATIVE_DOUBLE, n_n_arc_size);

    // ID'S OF FACES SURROUNDING NODE
    unsigned n_f_IDs_size=0;
    unsigned * n_f_IDs_h5 = ReadHDF5Dataset<unsigned>(&file_id, "/NODES/FRIENDS/FACES/ID/", H5T_NATIVE_UINT, n_f_IDs_size);

    // DIRECTION OF FACE RELATIVE TO CONTROL VOLUME
    unsigned n_f_dir_size=0;
    int * n_f_dir_h5 = ReadHDF5Dataset<int>(&file_id, "/NODES/FRIENDS/FACES/DIR/", H5T_NATIVE_INT, n_f_dir_size);

    // ID'S OF VERTICES SURROUNDING NODE
    unsigned n_v_IDs_size=0;
    unsigned * n_v_IDs_h5 = ReadHDF5Dataset<unsigned>(&file_id, "/NODES/FRIENDS/VERTICES/ID/", H5T_NATIVE_UINT, n_v_IDs_size);

    // AREAS SUBTENDED BY NODE AND VERTEX (USED IN CORIOLIS OPERATOR)
    unsigned n_v_area_size=0;
    double * n_v_area_h5 = ReadHDF5Dataset<double>(&file_id, "/NODES/FRIENDS/VERTICES/AREA/", H5T_NATIVE_DOUBLE, n_v_area_size);


    // Write node properties from the format in the hdf5 file to the appropriate arrays
    count = 0;
    for (unsigned i=0; i<n_lat_size; i++) {
        node_pos_sph(i, 0)  = n_lat_h5[i]*radConv;
        node_pos_sph(i, 1)  = n_lon_h5[i]*radConv;   
        node_fnum(i)        = n_n_num_h5[i];
        cv_area_sph(i)      = n_area_h5[i]*rr;          // CONVERT TO METRES^2
        cv_area_sph_r(i)    = 1.0/cv_area_sph(i);

        for (unsigned j=0; j<n_n_num_h5[i]; j++) {
            node_friends(i, j)      = n_n_IDs_h5[count];
            node_dists(i, j)        = n_n_arc_h5[count]*r;      // CONVERT TO METRES
            faces(i, j)             = n_f_IDs_h5[count];
            vertexes(i, j)          = n_v_IDs_h5[count];
            node_face_dir(i, j)     = n_f_dir_h5[count];
            node_vertex_area(i, j)  = n_v_area_h5[count]*rr;    // CONVERT OT METRES^2

            count++; 
        }
    }

    delete[] n_n_IDs_h5;
    delete[] n_n_arc_h5;
    delete[] n_lat_h5;
    delete[] n_lon_h5;
    delete[] n_n_num_h5;
    delete[] n_f_IDs_h5;
    delete[] n_area_h5;
    delete[] n_f_dir_h5;
    delete[] n_v_area_h5;
    delete[] n_v_IDs_h5;


    //-----------------------------------------------------------------------------------------
    //------------------------------ LOAD FACE INFO -------------------------------------------
    //-----------------------------------------------------------------------------------------

    // LATITUDE OF FACE CENTER
    unsigned f_lat_size=0;
    double * f_lat_h5 = ReadHDF5Dataset<double>(&file_id, "/FACES/LAT/", H5T_NATIVE_DOUBLE, f_lat_size);

    // LONGITUDE OF FACE CENTER
    unsigned f_lon_size=0;
    double * f_lon_h5 = ReadHDF5Dataset<double>(&file_id, "/FACES/LON/", H5T_NATIVE_DOUBLE, f_lon_size);

    // ARC LENGTH OF FACE
    unsigned f_arc_size=0;
    double * f_arc_h5 = ReadHDF5Dataset<double>(&file_id, "/FACES/LENGTH/", H5T_NATIVE_DOUBLE, f_arc_size);

    // FACE NORMAL VEC: LONGITUDE COMPONENT
    unsigned f_nlon_size=0;
    double * f_nlon_h5 = ReadHDF5Dataset<double>(&file_id, "/FACES/NORMAL_VEC_LON/", H5T_NATIVE_DOUBLE, f_nlon_size);

    // FACE NORMAL VEC: LATITUDE COMPONENT
    unsigned f_nlat_size=0;
    double * f_nlat_h5 = ReadHDF5Dataset<double>(&file_id, "/FACES/NORMAL_VEC_LAT/", H5T_NATIVE_DOUBLE, f_nlat_size);

    // INTERSECT LON WITH NODE-NODE ARC AND THE FACE
    unsigned f_lon_intersect_size=0;
    double * f_lon_intersect_h5 = ReadHDF5Dataset<double>(&file_id, "/FACES/INTERSECT_LON/", H5T_NATIVE_DOUBLE, f_lon_intersect_size);

    // INTERSECT LAT WITH NODE-NODE ARC AND THE FACE
    unsigned f_lat_intersect_size=0;
    double * f_lat_intersect_h5 = ReadHDF5Dataset<double>(&file_id, "/FACES/INTERSECT_LAT/", H5T_NATIVE_DOUBLE, f_lat_intersect_size);

    // LENGTH OF NODE-NODE ARC
    unsigned f_intersect_arc_size=0;
    double * f_intersect_arc_h5 = ReadHDF5Dataset<double>(&file_id, "/FACES/INTERSECT_LENGTH/", H5T_NATIVE_DOUBLE, f_intersect_arc_size);

    // IDS OF NODES EITHER SIDE OF FACE
    unsigned f_n_IDs_size=0;
    unsigned * f_n_IDs_h5 = ReadHDF5Dataset<unsigned>(&file_id, "/FACES/FRIENDS/NODES/ID/", H5T_NATIVE_UINT, f_n_IDs_size);

    // IDS OF VERTICES EITHER END OF FACE
    unsigned f_v_IDs_size=0;
    unsigned * f_v_IDs_h5 = ReadHDF5Dataset<unsigned>(&file_id, "/FACES/FRIENDS/VERTICES/ID/", H5T_NATIVE_UINT, f_v_IDs_size);
    
    // NUMBER OF FACE FRIENDS UPWIND
    unsigned f_f_num1_size=0;
    unsigned * f_f_num1_h5 = ReadHDF5Dataset<unsigned>(&file_id, "/FACES/FRIENDS/FACES1/FRIEND_NUM/", H5T_NATIVE_UINT, f_f_num1_size);

    // NUMBER OF FACE FRIENDS DOWNWIND
    unsigned f_f_num2_size=0;
    unsigned * f_f_num2_h5 = ReadHDF5Dataset<unsigned>(&file_id, "/FACES/FRIENDS/FACES2/FRIEND_NUM/", H5T_NATIVE_UINT, f_f_num2_size);

    // IDS OF UPWIND FACES
    unsigned f_f_ID1_size=0;
    unsigned * f_f_ID1_h5 = ReadHDF5Dataset<unsigned>(&file_id, "/FACES/FRIENDS/FACES1/ID/", H5T_NATIVE_UINT, f_f_ID1_size);

    // IDS OF DOWNWIND FACES
    unsigned f_f_ID2_size=0;
    unsigned * f_f_ID2_h5 = ReadHDF5Dataset<unsigned>(&file_id, "/FACES/FRIENDS/FACES2/ID/", H5T_NATIVE_UINT, f_f_ID2_size);

    // AREA OF FACE
    unsigned f_area_size=0;
    double * f_area_h5 = ReadHDF5Dataset<double>(&file_id, "/FACES/AREA/", H5T_NATIVE_DOUBLE, f_area_size);

    
    // Write node properties from the format in the hdf5 file to the appropriate arrays
    count = 0;
    unsigned count2 = 0;
    for (unsigned i=0; i<f_lat_size; i++) {
        face_centre_pos_sph(i, 0)       = f_lat_h5[i]*radConv;
        face_centre_pos_sph(i, 1)       = f_lon_h5[i]*radConv;   
        face_len(i)                     = f_arc_h5[i]*r;

        face_nodes(i, 0)                = f_n_IDs_h5[2*i];      // Upwind node    
        face_nodes(i, 1)                = f_n_IDs_h5[2*i+1];    // Downwind node

        face_vertexes(i, 0)             = f_v_IDs_h5[2*i];      // Should these be sorted into a "left" and "right" vertex?
        face_vertexes(i, 1)             = f_v_IDs_h5[2*i+1]; 

        face_normal_vec_map(i, 0)       = f_nlon_h5[i];
        face_normal_vec_map(i, 1)       = f_nlat_h5[i];

        face_intercept_pos_sph(i, 0)    = f_lat_intersect_h5[i]*radConv;
        face_intercept_pos_sph(i, 1)    = f_lon_intersect_h5[i]*radConv;

        face_node_dist(i)               = f_intersect_arc_h5[i]*r;
        face_node_dist_r(i)             = 1.0/face_node_dist(i);

        face_area(i)                    = f_area_h5[i]*rr;      // Convert to METRE^2

        // Need two seperate counters as a face can 
        // neighbour both a hexagon and a pentagon
        for (unsigned j=0; j<f_f_num1_h5[i]; j++) {
            face_friends(i, 0, j) = f_f_ID1_h5[count];
            count++;
        }
        for (unsigned j=0; j<f_f_num2_h5[i]; j++) {
            face_friends(i, 1, j) = f_f_ID2_h5[count2];
            count2++; 
        }

    }

    delete[] f_n_IDs_h5;
    delete[] f_v_IDs_h5;
    delete[] f_arc_h5;
    delete[] f_lat_h5;
    delete[] f_lon_h5;
    delete[] f_nlon_h5;
    delete[] f_nlat_h5;
    delete[] f_lat_intersect_h5;
    delete[] f_lon_intersect_h5;
    delete[] f_f_ID1_h5;
    delete[] f_f_ID2_h5;
    delete[] f_f_num1_h5;
    delete[] f_f_num2_h5;
    delete[] f_intersect_arc_h5;
    delete[] f_area_h5;

    //-----------------------------------------------------------------------------------------
    //------------------------------ LOAD VERTEX INFO -----------------------------------------
    //-----------------------------------------------------------------------------------------

    // VERTEX LAT
    unsigned v_lat_size=0;
    double * v_lat_h5 = ReadHDF5Dataset<double>(&file_id, "/VERTICES/LAT/", H5T_NATIVE_DOUBLE, v_lat_size);

    // VERTEX LON
    unsigned v_lon_size=0;
    double * v_lon_h5 = ReadHDF5Dataset<double>(&file_id, "/VERTICES/LON/", H5T_NATIVE_DOUBLE, v_lon_size);

    // VERTEX AREA
    unsigned v_area_size=0;
    double * v_area_h5 = ReadHDF5Dataset<double>(&file_id, "/VERTICES/AREA/", H5T_NATIVE_DOUBLE, v_area_size);

    // IDS OF FACES ATTACHED TO VERTEX
    unsigned v_f_IDs_size=0;
    unsigned * v_f_IDs_h5 = ReadHDF5Dataset<unsigned>(&file_id, "/VERTICES/FRIENDS/FACES/ID/", H5T_NATIVE_UINT, v_f_IDs_size);

    // DIRECTION OF FACES ATTACHED TO VERTEX
    unsigned v_f_dir_size=0;
    int * v_f_dir_h5 = ReadHDF5Dataset<int>(&file_id, "/VERTICES/FRIENDS/FACES/DIR/", H5T_NATIVE_INT, v_f_dir_size);

    // IDS OF NODES ATTACHED TO VERTEX
    unsigned v_n_IDs_size=0;
    unsigned * v_n_IDs_h5 = ReadHDF5Dataset<unsigned>(&file_id, "/VERTICES/FRIENDS/NODES/ID/", H5T_NATIVE_UINT, v_n_IDs_size);

    // SUBAREA THAT NODES SHARE WITH VERTEX
    unsigned v_n_subarea_size=0;
    double * v_n_subarea_h5 = ReadHDF5Dataset<double>(&file_id, "/VERTICES/FRIENDS/NODES/SUBAREA/", H5T_NATIVE_DOUBLE, v_n_subarea_size);

    // Write vertex properties from the format in the hdf5 file to the appropriate arrays
    count = 0;
    for (unsigned i=0; i<v_area_size; i++) {
        vertex_area(i)             = v_area_h5[i]*rr;               // CONVERT TO METRE^2
        vertex_area_r(i)           = 1.0/vertex_area(i);
        vertex_pos_sph(i, 0)       = v_lat_h5[i]*radConv;
        vertex_pos_sph(i, 1)       = v_lon_h5[i]*radConv;   

        for (unsigned j=0; j<3; j++) {
            vertex_faces(i, j)      = v_f_IDs_h5[count];
            vertex_nodes(i, j)      = v_n_IDs_h5[count];
            vertex_face_dir(i, j)   = v_f_dir_h5[count];

            // Area that vertex shares with each node, normalised by the area of the control volume
            vertex_R(i, j)          = v_n_subarea_h5[count]* rr * cv_area_sph_r( vertex_nodes(i, j) ); 

            count++; 
        }
    }

    count =0;
    for (unsigned i=0; i<NODE_NUM; i++) {
        int fnum = node_fnum(i);
    
        for (unsigned j=0; j<fnum; j++) {
            int v_ID = vertexes(i, j);
            centroid_pos_sph(i, j, 0) = vertex_pos_sph( v_ID, 0);
            centroid_pos_sph(i, j, 1) = vertex_pos_sph( v_ID, 1);
                
            count++; 
        }    
    }

    delete[] v_area_h5;
    delete[] v_lat_h5;
    delete[] v_lon_h5;
    delete[] v_f_IDs_h5;
    delete[] v_n_IDs_h5;
    delete[] v_f_dir_h5;

    H5Fclose(file_id);
   
    return 1;
};

template<typename T>
T * Mesh::ReadHDF5Dataset(hid_t * file_id, char const *dataset_name, hid_t H5TYPE, unsigned &data_size)
{
    hsize_t dims[1];    // Size 1 as we are assuming that every dataset is 1D

    hid_t dset_id = H5Dopen(*file_id, dataset_name, H5P_DEFAULT);
    hid_t dspace_id = H5Dget_space(dset_id);
    H5Sget_simple_extent_dims(dspace_id, dims, NULL);
    T * data = new T[ dims[0] ];
    H5Dread(dset_id, H5TYPE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(dset_id);
    H5Sclose(dspace_id);

    data_size = dims[0];
    return data;
};

// LEGACY FUNCTION
// Function to read in text file containing mesh information
// The four pieces of information read and stored are:
//      1. Node ID numbers
//      2. Node positions in spherical coords
//      3. Neighbouring node ID numbers
//      4. All centroid positions in spherical coords
int Mesh::ReadMeshFile(void)
{
    std::string line, val; // strings for column and individual number
    std::string file_str;  // string with path to mesh file.
    int i, node_id;

    // file_str = globals->path + SEP + "input_files" + SEP + "grid_l" + std::to_string(globals->geodesic_l.Value()) + ".txt";
    file_str = globals->path + SEP + "input_files" + SEP + "grid_l" + std::to_string(globals->geodesic_l.Value()) + ".txt";

    // in stream for input.in file
    std::ifstream gridFile(file_str, std::ifstream::in);

    if (gridFile.is_open())
    {
        outstring << std::endl
                  << "Found mesh file: " + file_str << std::endl;
        globals->Output->Write(OUT_MESSAGE, &outstring);

        std::getline(gridFile, line); // READ HEADER
        while (std::getline(gridFile, line))
        {
            // std::cout<<line<<std::endl;
            std::istringstream line_ss(line);
            std::getline(line_ss >> std::ws, val, ' '); // COL 0: Node ID number

            node_id = std::stoi(val);

            std::getline(line_ss >> std::ws, val, ' '); // COL 1: Node Latitude
            node_pos_sph(node_id, 0) = std::atof(val.c_str()) * radConv;

            std::getline(line_ss >> std::ws, val, ' '); // COL 2: Node Longitude
            node_pos_sph(node_id, 1) = std::atof(val.c_str()) * radConv;

            std::getline(line_ss >> std::ws, val, '{'); // Read up to friends open bracket
            for (i = 0; i < 5; i++)
            {
                std::getline(line_ss >> std::ws, val, ','); // Read each friend ID
                node_friends(node_id, i) = std::stoi(val);
            }
            std::getline(line_ss >> std::ws, val, '}'); // Read last friend ID
            node_friends(node_id, 5) = std::stoi(val);

            std::getline(line_ss >> std::ws, val, ',');
            std::getline(line_ss >> std::ws, val, '{'); // Read up to centroid list
            for (i = 0; i < 5; i++)
            {
                std::getline(line_ss >> std::ws, val, '('); // Read up to coord open bracket
                std::getline(line_ss >> std::ws, val, ','); // Read coord lat
                centroid_pos_sph(node_id, i, 0) = std::atof(val.c_str()) * radConv;

                std::getline(line_ss >> std::ws, val, ')'); // Read coord lon
                centroid_pos_sph(node_id, i, 1) = std::atof(val.c_str()) * radConv;
                std::getline(line_ss >> std::ws, val, ','); // Read end of first coord
            }
            std::getline(line_ss >> std::ws, val, '('); // Read up to last coord open bracket
            std::getline(line_ss >> std::ws, val, ','); // Read last coord lat
            centroid_pos_sph(node_id, 5, 0) = std::atof(val.c_str()) * radConv;

            std::getline(line_ss >> std::ws, val, ')'); // Read last coord lon
            centroid_pos_sph(node_id, 5, 1) = std::atof(val.c_str()) * radConv;

            std::getline(line_ss >> std::ws, val, '}'); // Finish line

            // UNCOMMENT THIS BLOCK TO VERIFY THAT MESH FILE IS BEING READ CORRECTLY
            // std::cout<<node_id<<'\t';
            // std::cout << node_pos_sph(node_id, 0) <<'\t'<< node_pos_sph(node_id, 1)<<'\t';
            // for (i=0; i<6; i++) {
            //   std::cout<<node_friends(node_id, i)<<'\t';
            // }
            // for (i=0; i<6; i++) {
            //   std::cout<<'('<<centroid_pos_sph(node_id, i, 0)<<", "<<centroid_pos_sph(node_id, i, 1)<<")\t";
            // }
            // std::cout<<std::endl;
        }

        gridFile.close();
    }
    else
    {
        outstring << "ERROR: GRID FILE NOT FOUND AT " + file_str << std::endl;
        globals->Output->Write(ERR_MESSAGE, &outstring);
        globals->Output->TerminateODIS();
    }

    return 1;
};

// int Mesh::ReadLatLonFile(void)
// {
//     int i, j, k, N_ll;
//
//     N_ll = (int)globals->dLat.Value();
//
//     const H5std_string DSET_CellID("cell_ID");
//     const H5std_string DSET_VInv("vandermonde_inv");
//     const H5std_string DSET_Rot("rotation");
//
//     std::string file_str;
//     file_str = globals->path + SEP
//                + "input_files" + SEP
//                + "grid_l" + std::to_string(globals->geodesic_l.Value())
//                + '_' + std::to_string(N_ll)
//                + 'x' + std::to_string(N_ll)
//                + "_test.h5";
//
//     // Define file name and open file
//     H5std_string FILE_NAME(file_str);
//     H5File file(FILE_NAME, H5F_ACC_RDONLY);
//
//     // Access dataspaces in latlon file
//     DataSet dset_cellID = file.openDataSet(DSET_CellID);
//     DataSet dset_vInv = file.openDataSet(DSET_VInv);
//     DataSet dset_rot = file.openDataSet(DSET_Rot);
//
//     // Create filespaces for the correct rank and dimensions
//     DataSpace fspace_cellID = dset_cellID.getSpace();
//     DataSpace fspace_vInv = dset_vInv.getSpace();
//     DataSpace fspace_rot = dset_rot.getSpace();
//
//     // Get number of dimensions in the files dataspace
//     int rank_cellID = fspace_cellID.getSimpleExtentNdims();
//     int rank_vInv = fspace_vInv.getSimpleExtentNdims();
//     int rank_rot = fspace_rot.getSimpleExtentNdims();
//
//     hsize_t dims_cellID[2];
//     hsize_t dims_vInv[3];
//     hsize_t dims_rot[1];
//
//     // Get size of each dimension
//     rank_cellID = fspace_cellID.getSimpleExtentDims( dims_cellID );
//     rank_vInv = fspace_vInv.getSimpleExtentDims( dims_vInv );
//     rank_rot = fspace_rot.getSimpleExtentDims( dims_rot );
//
//     // Create memoryspace to read the datasets
//     DataSpace mspace_cellID(2, dims_cellID);
//     DataSpace mspace_vInv(3, dims_vInv);
//     DataSpace mspace_rot(1, dims_rot);
//
//     // Create 1D arrays to store file data
//     int * cellID_1D;
//     double * vInv_1D;
//     double * rot_1D;
//
//     cellID_1D = new int[180/N_ll * 360/N_ll];
//     vInv_1D = new double[NODE_NUM * 6 * 6];
//     rot_1D = new double[NODE_NUM];
//
//     // Read in the data
//     dset_cellID.read( cellID_1D, PredType::NATIVE_INT, mspace_cellID, fspace_cellID );
//     dset_vInv.read( vInv_1D, PredType::NATIVE_DOUBLE, mspace_vInv, fspace_vInv );
//
//     // Load Array classes with 1D dynamic arrays
//     int count = 0;
//     for (i=0; i<NODE_NUM; i++)
//     {
//         for (j=0; j<6; j++)
//         {
//             for (k=0; k<6; k++)
//             {
//                 V_inv(i, j, k) = vInv_1D[count];
//                 count++;
//             }
//         }
//     }
//
//     count = 0;
//
//     int ID;
//     double lat1, lat2, lon1, lon2;
//     double *m, *x, *y;
//     double r;
//
//     m = new double;
//
//     double test_solution_gg[NODE_NUM];
//     double test_solution_ll[180/N_ll][360/N_ll];
//     r = 1.0;//globals->radius.Value();
//
//     for (i=0; i<180/N_ll; i++)
//     {
//         for (j=0; j<360/N_ll; j++)
//         {
//             cell_ID(i, j) = cellID_1D[count];
//             count++;
//
//             // get cell ID which contains current lat-lon grid point
//             ID = cell_ID(i, j);
//
//             // get sph coords of cell with current lat-lon node
//             lat1 = node_pos_sph(ID,0);
//             lon1 = node_pos_sph(ID,1);
//
//             // calculate regular lat-lon position in radians
//             lat2 = (90.0 - (double)(i*N_ll))*radConv;
//             lon2 = (double)(j*N_ll)*radConv;
//
//             // std::cout<<lat2/radConv<<' '<<lat1/radConv<<std::endl;
//             // std::cout<<lat1/radConv<<' '<<lon2/radConv<<' '<<lon1/radConv<<std::endl;
//
//             test_solution_ll[i][j] = cos(3. * lat2) * sin(5 * lon2);
//
//             // std::cout<<lon2<<' '<<lat2<<' '<<test_solution_ll[i][j]<<std::endl;
//
//             // set pointers to mapped coords of lat-lon node
//             x = &ll_map_coords(i, j, 0);
//             y = &ll_map_coords(i, j, 1);
//             // *m = 0.0;
//
//             // call mapping function to find x y of current lat-lon node
//             mapAtPoint(*m, *x, *y, lat1, lat2, lon1, lon2, r);
//         }
//     }
//
//     delete m;
//
//
//     delete[] cellID_1D;
//     delete[] vInv_1D;
//     delete[] rot_1D;
//
// };

// int Mesh::ReadWeightingFile(void)
// {
//     int i, N_ll;

//     N_ll = (int)globals->dLat.Value();

//     const H5std_string DSET_cols("column index");
//     const H5std_string DSET_rows("row index");
//     const H5std_string DSET_data("weights");

//     std::string file_str;
//     file_str = globals->path + SEP
//                + "input_files" + SEP
//                + "grid_l" + std::to_string(globals->geodesic_l.Value())
//                + '_' + std::to_string(N_ll)
//                + 'x' + std::to_string(N_ll)
//                + "_weights.h5";

//     // Define file name and open file
//     H5std_string FILE_NAME(file_str);
//     H5File file(FILE_NAME, H5F_ACC_RDONLY);

//     // Access dataspaces in latlon file
//     DataSet dset_cols = file.openDataSet(DSET_cols);
//     DataSet dset_rows = file.openDataSet(DSET_rows);
//     DataSet dset_data = file.openDataSet(DSET_data);

//     // Create filespaces for the correct rank and dimensions
//     DataSpace fspace_cols = dset_cols.getSpace();
//     DataSpace fspace_rows = dset_rows.getSpace();
//     DataSpace fspace_data = dset_data.getSpace();

//     // Get number of dimensions in the files dataspace
//     int rank_cols = fspace_cols.getSimpleExtentNdims();
//     int rank_rows = fspace_rows.getSimpleExtentNdims();
//     int rank_data = fspace_data.getSimpleExtentNdims();

//     hsize_t dims_cols[1];    // length no of non-zero elements
//     hsize_t dims_rows[1];
//     hsize_t dims_data[1];    // length no of non-zero elements

//     // Get size of each dimension
//     rank_cols = fspace_cols.getSimpleExtentDims( dims_cols );
//     rank_rows = fspace_rows.getSimpleExtentDims( dims_rows );
//     rank_data = fspace_data.getSimpleExtentDims( dims_data );

//     // Create memoryspace to read the datasets
//     DataSpace mspace_cols(1, dims_cols);
//     DataSpace mspace_rows(1, dims_rows);
//     DataSpace mspace_data(1, dims_data);

//     // Create 1D arrays to store file data
//     interpCols =    new    int[ dims_cols[0] ];
//     interpRows =    new    int[ dims_rows[0] ];
//     interpWeights = new double[ dims_data[0] ];

//     // Read in the data
//     dset_cols.read( interpCols, PredType::NATIVE_INT, mspace_cols, fspace_cols );
//     dset_rows.read( interpRows, PredType::NATIVE_INT, mspace_rows, fspace_rows );
//     dset_data.read( interpWeights, PredType::NATIVE_DOUBLE, mspace_data, fspace_data );

//     // Create  matrix handle using loaded data
//     // sparse_index_base_t index_type = SPARSE_INDEX_BASE_ZERO;     // we employ 0-based indexing.
//     // sparse_status_t err;

//     // sm1(rows,cols,nnz,outerIndexPtr, // read-write
//     //                            innerIndices,values);

//     int nrows = (360/N_ll)*(180/N_ll);
//     int ncols = 3*NODE_NUM;

//     interpMatrix = Eigen::Map<SpMat>(nrows, ncols, dims_data[0], interpRows, interpCols, interpWeights );

//     // std::cout<<interpMatrix.cols()<<' '<<dims_data[0]<<std::endl;

//     // interpMatrix.makeCompressed();

//     // interpMatrix = SpMat(nrows, ncols);

//     // // Allocate memory for non-zero entries
//     // interpMatrix.reserve(dims_data[0]);

//     // std::vector<Triplet> coefficients;
//     // coefficients.reserve(dims_data[0]);

//     // interpMatrix = new sparse_matrix_t;
//     // err = mkl_sparse_d_create_csr(interpMatrix, index_type, nrows, ncols, interpRows, interpRows+1, interpCols, interpWeights);

//     // // Sparse interpolation matrix successfully created. Now we must provide
//     // // additional information to the matrix handle for optimization purposes

//     // // Expected number of calls...?
//     // // sparse_operation_t operation = SPARSE_OPERATION_NON_TRANSPOSE;
//     // matrix_descr descrp;
//     // descrp.type = SPARSE_MATRIX_TYPE_GENERAL;
//     // descrp.mode = SPARSE_FILL_MODE_LOWER;
//     // descrp.diag = SPARSE_DIAG_NON_UNIT;

//     // sparse_operation_t operation = SPARSE_OPERATION_NON_TRANSPOSE;

//     // int numberOfExpectedCalls = 100000000;  // Expect a big number!
//     // err = mkl_sparse_set_dotmv_hint(*interpMatrix, operation, descrp, numberOfExpectedCalls);
//     // err = mkl_sparse_set_memory_hint (*interpMatrix, SPARSE_MEMORY_AGGRESSIVE);
//     // if (err>0) {
//     //     std::cout<<"Intel matrix optimization error!"<<std::endl;
//     // }

//     // // Now optimize matrix
//     // err = mkl_sparse_optimize(*interpMatrix);

//     return 1;

// };
