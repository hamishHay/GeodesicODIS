#include "mesh.h"
#include "globals.h"
#include "math.h"
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include "mathRoutines.h"
// #include "sphericalHarmonics.h"

// #include <mkl.h>
// #include <mkl_spblas.h>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cstring>

#include "H5Cpp.h"

#include <Eigen/Sparse>

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

//Mesh::Mesh():Mesh(2., 2.) {}; //default constructor

Mesh::Mesh(Globals * Globals, int N, int face_N, int vertex_N, int N_ll, int l_max)
   :node_pos_sph(N,2),
    node_pos_map(N,7,2),        // len 7 to include central node coords (even though it's zero)
    node_friends(N,6),
    node_dists(N,6),
    centroid_pos_sph(N,6,2),
    centroid_pos_map(N,6,2),
    node_m(N,7),                // len 7 to include central (parent) node
    centroid_m(N,6),            // len 6 as there is are only 5 or 6 centroids
    node_vel_trans(N,7,2),      // In the last dimension, element 1 is cos_a, element 2 is sin_a
    control_vol_edge_len(N,6),
    control_vol_edge_centre_pos_map(N,6,2),
    control_vol_edge_centre_m(N,6),
    control_vol_edge_normal_map(N,6,2),
    control_volume_surf_area_map(N),
    control_volume_mass(N),
    node_friend_element_areas_map(N,6,3),
    centroid_node_dists_map(N,6),

    trigLat(N,2),
    trigLon(N,2),
    trig2Lat(N,2),
    trig2Lon(N,2),
    trigSqLat(N,2),
    trigSqLon(N,2),

    Pbar_20(N),
    Pbar_22(N),

    // Pbar_lm(N),

    faces(N, 6),
    face_dir(face_N, 2),
    face_normal_vec_map(face_N, 2),
    face_normal_vec_dir(face_N),
    face_nodes(face_N, 2), //---> ID of nodes that belong to each face
    face_node_vel_trans(face_N, 2, 2),
    face_vertexes(face_N, 2), //---> ID of vertexes that belong to each face
    // face_friends(face_N, 2),
    face_interp_friends(face_N, 10),
    face_interp_weights(face_N, 10),
    face_node_pos_map(face_N, 2, 2),
    face_node_dist(face_N),
    face_centre_m(face_N, 2),
    face_centre_pos_sph(face_N, 2),
    face_intercept_pos_sph(face_N, 2),
    face_area(face_N),
    face_len(face_N),
    node_face_dir(N, 6),
    node_face_vel_trans(N, 6, 2),
    face_is_boundary(face_N),

    vertexes(N, 6),
    vertex_pos_sph(vertex_N, 2),
    vertex_R(vertex_N, 3),
    vertex_nodes(vertex_N, 3),
    vertex_faces(vertex_N, 3),
    node_R(N, 6),


    sh_matrix(N, (l_max+1)*(l_max+2) - 6), //-6 is to ignore degrees 0 and 1
    sh_matrix_fort(N*((l_max+1)*(l_max+2) - 6)),

    cell_is_boundary(N),


    trigMLon(N, l_max+1, 2)
{

    globals = Globals;                   // define reference to all constants

    node_num = globals->node_num;
    face_num = globals->face_num;
    vertex_num = globals->vertex_num;

    // Read in grid file
    ReadMeshFile();



    // Calculate mapping coordinates
    CalcMappingCoords();

    CalcNodeDists();

    

    // Find control volume edge lengths in mapping coords
    CalcControlVolumeEdgeLengths();

    // Find control volume edge midpoints in mapping coords
    CalcControlVolumeEdgeCentres();


    // Find control volume edge outward normal unit vectors
    CalcControlVolumeEdgeNormals();

    // Calculate velocity transform factors
    CalcVelocityTransformFactors();

    // Find control volume area for each node
    CalcControlVolumeArea();

    CalcControlVolumeMass();

    // Find the area of each sub triangle within each element
    CalcElementAreas();

    CalcCentNodeDists();

    // Evaluate trig functions at every node
    CalcTrigFunctions();

    CalcMaxTimeStep();

    CalcLegendreFuncs();

    CalcControlVolumeVertexR();

    AssignFaces();

    CalcVelocityTransformFactors();

    DefineBoundaryCells();

    CalcGradOperatorCoeffs();

    CalcDivOperatorCoeffs();

    // // CalcLaplaceOperatorCoeffs();

    CalcCoriolisOperatorCoeffs();

    CalcLinearDragOperatorCoeffs();

    // // GeneratePressureSolver();

    // GenerateMomentumOperator();

    if (globals->surface_type != FREE) {
        ReadWeightingFile();
    }

    CalcControlVolumeInterpMatrix();
    

    globals->Output->DumpGridData(this);
};

int Mesh::CalcControlVolumeInterpMatrix(void)
{

    interpSolvers = new Eigen::ColPivHouseholderQR<Eigen::MatrixXd>[node_num];

    int row_count = 0;
    int row_target = 6;
    for (int i=0; i<node_num; i++) {
        
        int f_num = 6;
        if ( node_friends(i, 5) == -1 ) f_num--;

        Vandermonde V(f_num+1, 6);

        int j = 0;
        for (j=0; j<f_num+1; j++) {
            double x, y; 
           
            x = node_pos_map(i, j, 0);
            y = node_pos_map(i, j, 1);

            V(j, 0) = 1.0;
            V(j, 1) = x;
            V(j, 2) = y;
            V(j, 3) = pow(x, 2.0);
            V(j, 4) = x*y;
            V(j, 5) = pow(y, 2.0); 
            
        }
        row_count += j;

        interpSolvers[i].compute(V);
        
    }

};

int Mesh::CalcControlVolumeVertexR(void)
{
    int vert_count = 0;
    int friend_num = 6;

    int node_count[node_num];

    for (int i=0; i<node_num; i++)
    {
        node_count[i] = 0;
        for (int j=0; j<friend_num; j++) {
            vertexes(i, j) = -1;
        }
    }

    for (int i=0; i<node_num; i++) {
        friend_num = 6;
        if (node_friends(i,5) < 0) friend_num--;

        for (int j=0; j<friend_num; j++) {
            if (vertexes(i, j) < 0) {
                vertexes(i, j) = vert_count;
                vertex_pos_sph(vert_count, 0) = centroid_pos_sph(i, j, 0);
                vertex_pos_sph(vert_count, 1) = centroid_pos_sph(i, j, 1);

                int f1 = node_friends(i, (j)%friend_num);
                int f2 = node_friends(i, (j+1)%friend_num);

                int f_num1 = 6;
                int f_num2 = 6;
                if (node_friends(f1,5) < 0) f_num1--;
                if (node_friends(f2,5) < 0) f_num2--;

                vertex_nodes(vert_count, 0) = i;
                vertex_nodes(vert_count, 1) = f1;
                vertex_nodes(vert_count, 2) = f2;

                node_count[i]++;
                for (int k=0; k<f_num1; k++) {
                    if (centroid_pos_sph(f1, k, 0) < vertex_pos_sph(vert_count, 0) +1e-9 && centroid_pos_sph(f1, k, 0) > vertex_pos_sph(vert_count, 0) -1e-9 &&
                        centroid_pos_sph(f1, k, 1) < vertex_pos_sph(vert_count, 1) +1e-9 && centroid_pos_sph(f1, k, 1) > vertex_pos_sph(vert_count, 1) -1e-9)

                    {
                        vertexes(f1, k) = vert_count;
                    }
                    //         {
                }

                for (int k=0; k<f_num2; k++) {
                    if (centroid_pos_sph(f2, k, 0) < vertex_pos_sph(vert_count, 0) +1e-9 && centroid_pos_sph(f2, k, 0) > vertex_pos_sph(vert_count, 0) -1e-9 &&
                        centroid_pos_sph(f2, k, 1) < vertex_pos_sph(vert_count, 1) +1e-9 && centroid_pos_sph(f2, k, 1) > vertex_pos_sph(vert_count, 1) -1e-9)

                    {
                        vertexes(f2, k) = vert_count;
                    }
                    //         {
                }

                vert_count++;
            }
        }
    }

    std::cout<<"VERTEX NUM: "<<vert_count<<' '<<vertex_num<<std::endl;
};

int Mesh::AssignFaces(void)
{
    struct face_friends
    {
        int ID;         // Global ID of face
        int rel_ID;     // ID of face relative to parent face
        double angle;   // Angle between face and parent face
    };


    int friend_num=6;
    int face_count;

    // Mark all faces as unassigned (<0)
    for (int i=0; i<node_num; i++) {
        for (int j=0; j<friend_num; j++) {
            faces(i, j) = -1;
        }
    }

    // Loop through all faces and assign:
    //  A unique ID (faces)
    //  The direction of the face normal (in or out)
    //  Distance between the nodes either side of the face (face_node_dist)
    face_count = 0;
    for (int i=0; i<node_num; i++) {
        friend_num = 6;
        if (node_friends(i,5) < 0) friend_num--;

        for (int j=0; j<friend_num; j++) {
            if (faces(i, j) < 0) {
                faces(i, j) = face_count;
                face_nodes(face_count, 0) = i;
                face_is_boundary(face_count) = 0;
                // face_friends(face_count, 0) = j;
                face_dir(face_count, 0) = 1;
                face_node_dist(face_count) = node_dists(i, j);
                node_face_dir(i, j) = 1;

                int f = node_friends(i, j);
                int friend_num2 = 6;
                if (node_friends(f,5) < 0) friend_num2--;

                double sph1[2], sph2[2], sph3[2], sph4[2], sphc[2];
                for (int k=0; k<2; k++) {
                    sph1[k] = centroid_pos_sph(i, (j+friend_num-1)%friend_num, k);
                    sph2[k] = centroid_pos_sph(i, (j)%friend_num, k);
                }

                double dist=0.0;
                double r = globals->radius.Value();
                distanceBetweenSph(dist, sph1[0], sph2[0], sph1[1], sph2[1], r);

                face_len(face_count) = fabs(dist);

                // get position of face centre
                midpointBetweenSph(sph1, sph2, sphc);

                if (sphc[1] < 0.0) sphc[1] += 2*pi;

                face_centre_pos_sph(face_count, 0) = sphc[0];
                face_centre_pos_sph(face_count, 1) = sphc[1];

                double lat1, lon1, lat2, lon2;
                double m = 0.0;
                double x1 = 0.0;
                double y1 = 0.0;
                double x2 = 0.0;
                double y2 = 0.0;

                lat1 = sphc[0];           lon1 = sphc[1];
                lat2 = node_pos_sph(i,0); lon2 = node_pos_sph(i, 1);

                mapAtPoint(m, x1, y1, lat1, lat2, lon1, lon2, r);
                face_centre_m(face_count, 0) = m;

                lat2 = node_pos_sph(f,0); lon2 = node_pos_sph(f, 1);

                mapAtPoint(m, x2, y2, lat1, lat2, lon1, lon2, r);
                face_centre_m(face_count, 1) = m;

                face_node_dist(face_count) = sqrt( pow(x2-x1, 2.0) + pow(y2-y1, 2.0));

                double normal_vec[2];
                normalVectorBetweenMap(normal_vec, sph1, sph2);
                face_normal_vec_map(face_count, 0) = normal_vec[0];
                face_normal_vec_map(face_count, 1) = normal_vec[1];

                if (normal_vec[0] > 0 && control_vol_edge_normal_map(i, (j+friend_num-1)%friend_num, 0) > 0) face_dir(face_count, 0) = 1;
                else face_dir(face_count, 0) = -1;

                // Calculate intersection point between the face length vector
                // and the vector between the shared nodes
                for (int k=0; k<2; k++) {
                    sph3[k] = node_pos_sph(i, k);
                    sph4[k] = node_pos_sph(f, k);
                }

                intersectPointSph(sph1, sph2, sph3, sph4, sphc);

                face_intercept_pos_sph(face_count, 0) = sphc[0];
                face_intercept_pos_sph(face_count, 1) = sphc[1];

                // Find friend on opposite side of the face
                // and update the values accordingly
                int k=0;
                for (int j2=0; j2<friend_num2; j2++) {
                    if (node_friends(f, j2) == i) {
                        faces(f, j2) = face_count;
                        face_dir(face_count, 1) = -face_dir(face_count, 0);
                        face_nodes(face_count, 1) = f;
                        node_face_dir(f, j2) = -1;
                        // face_dir(f, j2) = -1;

                    }
                }

                face_count++;
            }
        }
    }

    // Now each face must be looped over in order to calculate interpolation weights
    // This means finding the vertices and faces that are joined to the parent face, 
    // in a anticlockwise order
    for (int i=0; i<face_num; i++)
    {
        int j_add = 0;
        for (int k=0; k<2; k++)
        {
            int node_ID = face_nodes(i, k);

            friend_num = 6;
            if (node_friends(node_ID, 5) < 0) friend_num--;

            // Get list of faces attached to node
            std::vector<face_friends> friends_list(friend_num);
            // friends_list[0].ID = i;
            // friends_list[0].angle = 360.0;

            // Loop through the friend faces, calculating the angle between
            int j = 0,  added = 0;
            while (added < friend_num) {
                int ID = faces(node_ID, j);
                friends_list[added].ID = ID;

                double angle = 0.0;
                double sph1[2], sph2[2], sphc[2];

                for (int ii=0; ii<2; ii++) {
                    sph1[ii] = face_intercept_pos_sph(i,ii);
                    sph2[ii] = face_intercept_pos_sph(ID,ii);
                    sphc[ii] = node_pos_sph(node_ID, ii);
                }
                angleBetweenMap(sph1, sph2, sphc, angle);

                friends_list[added].angle = angle*180./pi;
                if (friends_list[added].angle > 0.0+1e-8) friends_list[added].angle -= 360.0;

                added++;
                j++;
            }

            // Get list of vertexes attached to node
            // Loop through the friend vertexes, calculating the angle between
            std::vector<face_friends> face_vertex_list(friend_num);
            for (j=0; j<friend_num; j++) {
                face_vertex_list[j].ID = vertexes(node_ID, j);
                face_vertex_list[j].rel_ID = j;

                double angle = 0.0;
                double sph1[2], sph2[2], sphc[2];

                for (int ii=0; ii<2; ii++) {
                    sph1[ii] = face_intercept_pos_sph(i, ii);
                    sph2[ii] = vertex_pos_sph(vertexes(node_ID, j), ii);
                    sphc[ii] = node_pos_sph(node_ID, ii);
                }
                angleBetweenMap(sph1, sph2, sphc, angle);

                face_vertex_list[j].angle = angle*180./pi;
                if (face_vertex_list[j].angle > 0.0+1e-8) face_vertex_list[j].angle -= 360.0;
            }

            // order faces relative to face i in an anti-clockwise direction
            std::sort( friends_list.begin( ), friends_list.end( ), [ ]( const face_friends &lhs, const face_friends &rhs )
            { return lhs.angle > rhs.angle; });

            // order vertexes relative to face i in an anti-clockwise direction
            std::sort( face_vertex_list.begin( ), face_vertex_list.end( ), [ ]( const face_friends &lhs, const face_friends &rhs )
            { return lhs.angle > rhs.angle; });

            // Now we know the right order of faces and vertexes relative to face
            // i, so we can now compute the R weights
            double area_cv=0.0;
            for (j=0; j<friend_num; j++) {
                double a1=0.0, a2=0.0, area=0.0;
                double r = globals->radius.Value();

                // get intersect position
                int v_ID = face_vertex_list[j].ID;
                int v_ID2 = face_vertex_list[(j+1)%friend_num].ID;

                double sph_int[2], sph1[2], sph2[2];

                sph1[0] = node_pos_sph(node_ID, 0);
                sph1[1] = node_pos_sph(node_ID, 1);

                sph2[0] = vertex_pos_sph(v_ID, 0);
                sph2[1] = vertex_pos_sph(v_ID, 1);

                // calc triangle area for face_intersect1-node-vertex
                sph_int[0] = vertex_pos_sph(v_ID2, 0);
                sph_int[1] = vertex_pos_sph(v_ID2, 1);


                triangularAreaSph(a1, sph1[0], sph2[0], sph_int[0], sph1[1], sph2[1], sph_int[1], r);

                // sum areas and divide by cv area to get R_iv
                area_cv += a1;

            }

            double area_t=0.0;
            std::vector<double> R_weights(friend_num);
            for (j=0; j<friend_num; j++) {
                double a1=0.0, a2=0.0, area=0.0;
                double r = globals->radius.Value();

                // get intersect position
                int v_ID = face_vertex_list[j].ID;
                int f_ID = friends_list[j].ID;
                double sph_int[2], sph1[2], sph2[2];

                if (j==1) {
                    face_vertexes(f_ID, 0) = v_ID;
                    face_vertexes(f_ID, 1) = face_vertex_list[j-1].ID;
                }

                sph1[0] = node_pos_sph(node_ID, 0);
                sph1[1] = node_pos_sph(node_ID, 1);

                sph2[0] = vertex_pos_sph(v_ID, 0);
                sph2[1] = vertex_pos_sph(v_ID, 1);

                // calc triangle area for face_intersect1-node-vertex
                sph_int[0] = face_intercept_pos_sph(f_ID, 0);
                sph_int[1] = face_intercept_pos_sph(f_ID, 1);


                triangularAreaSph(a1, sph1[0], sph2[0], sph_int[0], sph1[1], sph2[1], sph_int[1], r);

                // calc triangle area for face_intersect2-node-vertex
                f_ID = friends_list[(j+1)%friend_num].ID;
                sph_int[0] = face_intercept_pos_sph(f_ID, 0);
                sph_int[1] = face_intercept_pos_sph(f_ID, 1);


                triangularAreaSph(a2, sph1[0], sph2[0], sph_int[0], sph1[1], sph2[1], sph_int[1], r);

                area = (a1+a2)/area_cv;

                R_weights[j] = area;

                area_t += area;

                bool added = false;
                int count = 0;
                while (!added) {
                    if (vertex_nodes(v_ID, count) == node_ID) {
                        vertex_R(v_ID, count) = R_weights[j];
                        added = true;
                    }
                    count++;
                }

                // if vertex_nodes(v_ID, 0 or 1 or 2) == node_ID, vertex_R(v_ID, 0 or 1 or 2) = R_weights[j]

            }


            for (j=0; j<friend_num; j++) {

                //Calc each weight w_ee' for e
                if (i != friends_list[j].ID) {
                    face_interp_friends(i, j_add) = friends_list[j].ID;
                    face_interp_weights(i, j_add) = 0.0;
                    for (int j2=0; j2<j; j2++)
                    {
                        face_interp_weights(i, j_add) += R_weights[j2];
                    }
                    face_interp_weights(i, j_add) -= 0.5;


                    for (int j2=0; j2<friend_num; j2++)
                    {
                        if (friends_list[j].ID == faces(node_ID, j2)) {
                            face_interp_weights(i, j_add) *= node_face_dir(node_ID, j2);
                        }
                    }

                    // THIS STEP IS VERY IMPORTANT
                    // The edge under interpolation may point in the positive
                    // or negative vorticity direction; this ensures that this
                    // happens under the asusmption that the weights are
                    // calculated in an anti-clockwise sense
                    int t_ev = 1;
                    if (face_nodes(i, 0) == node_ID) t_ev = 1;
                    else t_ev = -1;

                    face_interp_weights(i, j_add) *= t_ev;

                    j_add++;
                }
            }
        }
    }

    for (int i=0; i<vertex_num; i++) {
        vertex_faces(i,0) = -1;
        vertex_faces(i,1) = -1;
        vertex_faces(i,2) = -1;

    }

    for (int i=0; i<face_num; i++)
    {
        int v1 = face_vertexes(i,0);
        int v2 = face_vertexes(i,1);

        for (int j=0; j<3; j++)
        {
            if (vertex_faces(v1, j) < 0)
            {
                vertex_faces(v1, j) = i;
                // vertex_face_dir(v1, j) = ?

                double sphf[2], sphv[2];
                sphf[0] = face_centre_pos_sph(i, 0);
                sphf[1] = face_centre_pos_sph(i, 1);
                sphv[0] = vertex_pos_sph(v1, 0);
                sphv[1] = vertex_pos_sph(v1, 1);

                double lat1, lon1, lat2, lon2;
                double m = 0.0;
                double xv = 0.0;
                double yv = 0.0;
                double fnx = 0.0;
                double fny = 0.0;

                lat1 = sphc[0];           lon1 = sphc[1];
                lat2 = sphv[0];           lon2 = sphv[1]);

                mapAtPoint(m, xv, yv, lat1, lat2, lon1, lon2, r);

                fnx = face_normal_vec_map(i, 0);
                fny = face_normal_vec_map(i, 1);

                double cross = xv*fny - yv*fnx;

                // take cross product between face normal 
                // and vertex vector. The sign indicates if 
                // the face normal points in the clockwise or
                // anticlockwise direction!



                if (cross > 0) vertex_face_dir(v1, j) = 1; // clockwise
                else vertex_face_dir(v1, j) = -1;          // anti-clockwise

                break;
            } 
        }

        for (int j=0; j<3; j++)
        {
           int face_ID = vertex_faces(v1, j);
           int node_ID = face_nodes(i, 0);

           //check if vertex has inner or outer 
        }
        
    }

    for (int i=0; i<face_num; i++)
    {
        int v1 = face_vertexes(i,0);
        int v2 = face_vertexes(i,1);

        double sph1[2], sph2[2], sph3[2];

        sph1[0] = vertex_pos_sph(v1, 0);
        sph1[1] = vertex_pos_sph(v1, 1);

        sph2[0] = vertex_pos_sph(v2, 0);
        sph2[1] = vertex_pos_sph(v2, 1);

        double a=0.0, area=0.0;
        double r = globals->radius.Value();

        for (int k=0; k<2; k++)
        {
            int n1 = face_nodes(i, k);

            sph3[0] = node_pos_sph(n1, 0);
            sph3[1] = node_pos_sph(n1, 1);

            triangularAreaSph(a, sph1[0], sph2[0], sph3[0], sph1[1], sph2[1], sph3[1], r);

            area += a;
        }

        face_area(i) = area;

        face_area(i) = face_node_dist(i)*face_len(i)/2.0;
    }

    std::cout<<"FACE NUM: "<<face_count<<' '<<face_num<<std::endl;
    std::cout<<"NODE NUM: "<<node_num<<std::endl;
}

int Mesh::DefineBoundaryCells(void)
{
    // simply define if a control volume is either interior to a boundary (0),
    // shares a face with a boundary (1), or is completely outside of the
    // boundary (2).
    int i, j, f, friend_num;
    double lat1, lat2, lon1, lon2, dist;
    double r = 1.0;

    // r = globals->radius.Value();

    for (i=0; i<node_num; i++)
    {
        f = node_friends(i,5);

        friend_num = 6;                                     // Assume hexagon (6 centroids)
        if (f == -1) {
            friend_num = 5;                                 // Check if pentagon (5 centroids)
            node_dists(i,5) = -1.0;
        }

        lat1 = 0.0;                     // set map coords for first
        lon1 = pi;

        lat2 = node_pos_sph(i, 0);                     // set map coords for second centroid
        lon2 = node_pos_sph(i, 1);


        distanceBetweenSph(dist, lat1, lat2, lon1, lon2, r);
        // dist /= r;

        // std::cout<<dist*180./pi<<std::endl;
        if (dist*180./pi <= 90) cell_is_boundary(i) = 0;
        else cell_is_boundary(i) = 1;

        // if (cell_is_boundary(i) == 0)
        // {
        //   for (j=0; j<friend_num; j++)                        // Loop through all centroids in the control volume
        //   {
        //     int f_ID = faces(i,j);
        //     lat2 = face_intercept_pos_sph(f_ID, 0);                     // set map coords for second centroid
        //     lon2 = face_intercept_pos_sph(f_ID, 1);
        //
        //     distanceBetweenSph(dist, lat1, lat2, lon1, lon2, r);   //
        //     // dist /= r;
        //     // std::cout<<i<<'\t'<<dist*180./pi<<std::endl;
        //     if (dist*180./pi < 20) {
        //       // cell_is_boundary(i) = 1;
        //       face_is_boundary(f_ID) = 1;
        //       // break;
        //     }
        //   }
        // }
        // if (cell_is_boundary(i) == 2) std::cout<<i<<'\t'<<cell_is_boundary(i)<<std::endl;

        cell_is_boundary(i) = 0;
    }



}

// Function to calculate the cosine(alpha) and sine(alpha) velocity tranform
// factors from Lee and Macdonald (2009). These factors are constant in time,
// and so are pre-calculated here to be easily accessed during the spatial
// integration.
int Mesh::CalcVelocityTransformFactors(void)
{
    int i, j, f;
    double * cos_a, * sin_a;
    double lat1, lat2, lon1, lon2;

    for (i=0; i<node_num; i++)
    {
        lat1 = node_pos_sph(i,0);
        lon1 = node_pos_sph(i,1);

        // Set pointers to address of variables we want to change
        cos_a = &node_vel_trans(i,0,0);
        sin_a = &node_vel_trans(i,0,1);

        // Pass pointers by reference
        velTransform(*cos_a, *sin_a, lat1, lat1, lon1, lon1);

        // Loop through all friends (len is 7 as first element is parent node)
        // Assign node mapping factors and coordinates
        for (j = 1; j<7; j++)
        {
            cos_a = &node_vel_trans(i,j,0);
            sin_a = &node_vel_trans(i,j,1);

            f = node_friends(i,j-1);            // j-1 because len node_friends[i] = 6;
            switch (f) {
            case -1:                          // if node is pentagon center
                *cos_a = -1.0;
                *sin_a = -1.0;
                break;
            default:                          // if node is a hexagon center
                lat2 = node_pos_sph(f, 0);
                lon2 = node_pos_sph(f, 1);

                // assign transform factors values to arrays
                velTransform(*cos_a, *sin_a, lat1, lat2, lon1, lon2);

                break;
            }
        }
    }

    for (i = 0; i < face_num; i++)
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
    double * x, * y, * m, r;
    double lat1, lat2, lon1, lon2;

    r = globals->radius.Value();

    for (i=0; i<node_num; i++)
    {
        lat1 = node_pos_sph(i,0);
        lon1 = node_pos_sph(i,1);

        // Set pointers to address of variables we want to change
        m = &node_m(i,0);
        x = &node_pos_map(i,0,0);
        y = &node_pos_map(i,0,1);

        // Pass pointers by reference
        mapAtPoint(*m, *x, *y, lat1, lat1, lon1, lon1, r);

        // Loop through all friends (len is 7 as first element is parent node)
        // Assign node mapping factors and coordinates
        for (j = 1; j<7; j++)
        {
            m = &node_m(i,j);
            x = &node_pos_map(i,j,0);
            y = &node_pos_map(i,j,1);

            f = node_friends(i,j-1);            // j-1 because len node_friends[i] = 6;
            switch (f) {
                case -1:                          // if node is pentagon center
                    *m = -1.0;
                    *x = -1.0;
                    *y = -1.0;
                    break;
                default:                          // if node is a hexagon center
                    lat2 = node_pos_sph(f, 0);
                    lon2 = node_pos_sph(f, 1);

                    // assign mapped values to arrays
                    mapAtPoint(*m, *x, *y, lat1, lat2, lon1, lon2, r);

                    // if (i==137 && (f==1427 || f==1453 || f ==1457))
                    // {
                    //   double dist = 0.0;
                    //   std::cout<<i<<' '<<f<<' '<<*x<<' '<<*y<<' '<<lat2*180/pi<<std::endl;
                    //   // distanceBetween(dist, *xc, *x1, *yc, *y1);
                    //   // std::cout<<i<<' '<<f<<' '<<dist<<std::endl;
                    //   // std::cout<<i<<' '<<nx<<' '<<ny<<std::endl;
                    //   // std::cout<<i<<' '<<m<<std::endl;
                    //   // std::cout<<i<<' '<<mesh->node_dists(i, (j+1)%friend_num)<<std::endl;
                    //   // for (int k=0; k<3; k++) std::cout<<' '<<f1<<' '<<(*element_areas)(i,(j)%friend_num,k);
                    //   // for (int k=0; k<3; k++) std::cout<<' '<<f2<<' '<<(*element_areas)(i,(j+1)%friend_num,k);
                    //   // for (int k=0; k<3; k++) std::cout<<' '<<f2<<' '<<(*element_areas)(i,(j+2)%friend_num,k);
                    // }
                    // if (i==1453 && (f==1457 || f==137 || f ==1427))
                    // {
                    //   double dist = 0.0;
                    //   // distanceBetween(dist, *xc, *x1, *yc, *y1);
                    //   // std::cout<<i<<' '<<f<<' '<<dist<<std::endl;
                    //   // std::cout<<i<<' '<<j<<' '<<f<<std::endl;
                    //   std::cout<<i<<' '<<f<<' '<<*x<<' '<<*y<<' '<<lat2*180/pi<<std::endl;
                    // }
                    break;
            }
        }
    }

    for (i=0; i<node_num; i++)
    {
        lat1 = node_pos_sph(i,0);
        lon1 = node_pos_sph(i,1);

        // Assign centroid mapping factors and coordinates
        for (j = 0; j<6; j++)
        {
            m = &centroid_m(i,j);
            x = &centroid_pos_map(i,j,0);
            y = &centroid_pos_map(i,j,1);

            f = node_friends(i,j);
            switch (f) {
            case -1:                          // if node is pentagon center
                *m = -1.0;
                *x = -1.0;
                *y = -1.0;
                break;
            default:                          // if node is a hexagon center
                lat2 = centroid_pos_sph(i, j, 0);
                lon2 = centroid_pos_sph(i, j, 1);

                // assign mapped values to arrays
                mapAtPoint(*m, *x, *y, lat1, lat2, lon1, lon2, r);
                break;
            }
        }
    }

    return 1;
};

// Function to find the length of control volume edges for each node, in mapping
// coordinates. Values are stored in control_vol_edge_len 2D array.
int Mesh::CalcControlVolumeEdgeLengths(void)
{
    int i, j, f, friend_num;
    double * x1, * y1, * x2, * y2, * edge_len;

    for (i=0; i<node_num; i++)
    {

        f = node_friends(i,5);

        friend_num = 6;                                     // Assume hexagon (6 centroids)
        if (f == -1) {
            friend_num = 5;                                   // Check if pentagon (5 centroids)
            control_vol_edge_len(i,5) = -1.0;
        }

        for (j=0; j<friend_num; j++)                        // Loop through all centroids in the control volume
        {
            edge_len = &control_vol_edge_len(i,j);            // set pointer to edge length array

            x1 = &centroid_pos_map(i, j, 0);                  // set map coords for first centroid
            y1 = &centroid_pos_map(i, j, 1);

            x2 = &centroid_pos_map(i, (j+1)%friend_num, 0);   // set map coords for second centroid
            y2 = &centroid_pos_map(i, (j+1)%friend_num, 1);   // automatically loops around using %

            distanceBetween(*edge_len, *x1, *x2, *y1, *y2);   // calculate distance between the two centroids.
            // Edge_len automatically assigned the length
        }
    }

    return 1;
};

// Function to find the spherical distance between a node and each of
// its neighbours.
int Mesh::CalcNodeDists(void)
{
    int i, j, f, friend_num;
    double * lat1, * lat2, * lon1, * lon2, * edge_len;
    double r;

    r = globals->radius.Value();

    for (i=0; i<node_num; i++)
    {

        f = node_friends(i,5);

        friend_num = 6;                                     // Assume hexagon (6 centroids)
        if (f == -1) {
            friend_num = 5;                                 // Check if pentagon (5 centroids)
            node_dists(i,5) = -1.0;
        }

        lat1 = &node_pos_sph(i, 0);                     // set map coords for first
        lon1 = &node_pos_sph(i, 1);
        for (j=0; j<friend_num; j++)                        // Loop through all centroids in the control volume
        {
            edge_len = &node_dists(i,j);                    // set pointer to edge length array

            f = node_friends(i,j);

            lat2 = &node_pos_sph(f, 0);                     // set map coords for second centroid
            lon2 = &node_pos_sph(f, 1);

            distanceBetweenSph(*edge_len, *lat1, *lat2, *lon1, *lon2, r);   // calculate distance between the two centroids.
            // Edge_len automatically assigned the length

            // std::cout<<*edge_len<<std::endl;
        }
    }


    return 1;
};


int Mesh::CalcMaxTimeStep(void)
{
  int i, j, f, friend_num;
  double dist;
  double g, h_max, dt;

  g = globals->g.Value();
  h_max = globals->h.Value();


  dt = globals->timeStep.Value();

  if (globals->surface_type == FREE ||
      globals->surface_type == FREE_LOADING ||
      globals->surface_type == LID_LOVE)
  {
      if (globals->surface_type == LID_LOVE) g *= -(globals->shell_factor_beta[globals->l_max.Value()]-1.0);
      // std::cout<<g<<std::endl;
      for (i=0; i<node_num; i++)
      {
          f = node_friends(i,5);
          friend_num = 6;                                     // Assume hexagon (6 centroids)
          if (f == -1) {
              friend_num = 5;                                   // Check if pentagon (5 centroids)
          }
          for (j=0; j<friend_num; j++)                      // Loop through all centroids in the control volume
          {
              dist = node_dists(i,j)*0.25;                 // consider one quater node-node distance
              dt = std::min(dt, dist/sqrt(g*h_max));
          }
      }

      dt *= 0.85;         // take some caution
  }

  std::cout<<"DT: "<<dt<<std::endl;

  globals->timeStep.SetValue(dt);



  return 1;
}


// Function to find the length of control volume edges for each node, in mapping
// coordinates. Values are stored in control_vol_edge_len 2D array.
int Mesh::CalcCentNodeDists(void)
{
    int i, j, f, friend_num;
    double * x1, * y1, * x2, * y2, * edge_len;

    for (i=0; i<node_num; i++)
    {

        f = node_friends(i,5);

        friend_num = 6;                                     // Assume hexagon (6 centroids)
        if (f == -1) {
            friend_num = 5;                                   // Check if pentagon (5 centroids)
            centroid_node_dists_map(i,5) = -1.0;
        }

        for (j=0; j<friend_num; j++)                        // Loop through all centroids in the control volume
        {
            f = node_friends(i,j);

            edge_len = &centroid_node_dists_map(i,j);            // set pointer to edge length array

            x1 = &node_pos_map(i, 0, 0);   // set map coords for second centroid
            y1 = &node_pos_map(i, 0, 1);   // automatically loops around using %

            x2 = &centroid_pos_map(i, j%friend_num, 0);   // set map coords for second centroid
            y2 = &centroid_pos_map(i, j%friend_num, 1);   // automatically loops around using %

            distanceBetween(*edge_len, *x1, *x2, *y1, *y2);   // calculate distance between the two centroids.
            // Edge_len automatically assigned the length
        }
    }

    return 1;
};

// Function to find the midpoint of control volume edges for each node, in
// mapping coordinates. Values are stored in control_vol_edge_centre_pos_map 2D
// array.
int Mesh::CalcControlVolumeEdgeCentres(void)
{
    int i, j, f, friend_num;
    double * x1, * y1, * x2, * y2, * xc, * yc, * m;
    double lat1, lat2, lat3, lon1, lon2, lon3;
    double r = globals->radius.Value();

    for (i=0; i<node_num; i++)
    {

        f = node_friends(i,5);

        friend_num = 6;                                     // Assume hexagon (6 centroids)
        if (f == -1) {
            friend_num = 5;                                   // Check if pentagon (5 centroids)
            control_vol_edge_len(i,5) = -1.0;
            control_vol_edge_centre_m(i,5) = -1.0;
        }

        for (j=0; j<friend_num; j++)                        // Loop through all centroids in the control volume
        {
            f = node_friends(i,j);

            xc = &control_vol_edge_centre_pos_map(i,j,0);     // set pointer to edge length array
            yc = &control_vol_edge_centre_pos_map(i,j,1);

            x1 = &centroid_pos_map(i, j, 0);                  // set map coords for first centroid
            y1 = &centroid_pos_map(i, j, 1);

            x2 = &centroid_pos_map(i, (j+1)%friend_num, 0);   // set map coords for second centroid
            y2 = &centroid_pos_map(i, (j+1)%friend_num, 1);   // automatically loops around using %

            midpointBetween(*xc, *yc, *x1, *x2, *y1, *y2);    // calculate center coords between the two centroids.
            // xc and yc automatically assigned the coords

            m = &control_vol_edge_centre_m(i,j);              // assign pointer to midpoint map factor

            // lat1 = centroid_pos_sph(i, j, 0);                 // get first centroid coords
            // lon1 = centroid_pos_sph(i, j, 1);

            lat1 = node_pos_sph(i, 0);                 // get first centroid coords
            lon1 = node_pos_sph(i, 1);

            lat2 = centroid_pos_sph(i, j, 0);  // get second centroid coords
            lon2 = centroid_pos_sph(i, j, 1);

            // lat3 = centroid_pos_sph(i, (j+1)%friend_num, 0);  // get second centroid coords
            // lon3 = centroid_pos_sph(i, (j+1)%friend_num, 1);
            //
            // lat2 = 0.5*(lat2 + lat3);
            // lon2 = 0.5*(lon2 + lon3);

            mapFactorAtPoint(*m, lat1, lat2, lon1, lon2);     // calculate map factor and store


            // Calculate here as the midpoint to the cv edge is only known in
            // mapped coordinates
            *m = (4.*pow(r,2.0) + pow(*xc, 2.0) + pow(*yc, 2.0))/(4.*pow(r,2.0));

            // std::cout<<*m<<'\t'<<(4.*r*r + (*xc)*(*xc)+ (*yc)*(*yc))/(4.*r*r)<<std::endl;
        }
    }

    return 1;
};

// Function to find the unit normal vector to the control volume edge, in
// mapping coordinates. Values are stored in control_vol_edge_normal_map 2Dd
// array.
int Mesh::CalcControlVolumeEdgeNormals(void)
{
    int i, j, f, friend_num;
    double * x1, * y1, * x2, * y2, * xn, * yn;

    for (i=0; i<node_num; i++)
    {

        f = node_friends(i,5);

        friend_num = 6;                                     // Assume hexagon (6 centroids)
        if (f == -1) {
            friend_num = 5;                                   // Check if pentagon (5 centroids)
            control_vol_edge_len(i,5) = -1.0;
            control_vol_edge_normal_map(i,5,0) = -1.0;        // set pointer to edge length array
            control_vol_edge_normal_map(i,5,1) = -1.0;
        }



        for (j=0; j<friend_num; j++)                        // Loop through all centroids in the control volume
        {

            xn = &control_vol_edge_normal_map(i,j,0);         // set pointer to edge length array
            yn = &control_vol_edge_normal_map(i,j,1);

            x1 = &centroid_pos_map(i, j, 0);                  // set map coords for first centroid
            y1 = &centroid_pos_map(i, j, 1);

            x2 = &centroid_pos_map(i, (j+1)%friend_num, 0);   // set map coords for second centroid
            y2 = &centroid_pos_map(i, (j+1)%friend_num, 1);   // automatically loops around using %

            // double g2 = (*y2-*y1)/(*x2-*x1);
            // double c2 = *y2 - g2*(*x2);

            normalVectorBetween(*xn, *yn, *x1, *x2, *y1, *y2);    // calculate center coords between the two centroids.

            // x1 = &node_pos_map(i, 0, 0);                  // set map coords for first node_pos_map
            // y1 = &node_pos_map(i, 0, 1);
            //
            // x2 = &node_pos_map(i, (j+2)%friend_num, 0);   // set map coords for second node_pos_map
            // y2 = &node_pos_map(i, (j+2)%friend_num, 1);   // automatically loops around using %
            // double xn2 = *x2/std::sqrt((*x2)*(*x2)+(*y2)*(*y2));
            // double yn2 = *y2/std::sqrt((*x2)*(*x2)+(*y2)*(*y2));
            //
            // double g1 = (*y2)/(*x2);
            // double c1 = *y2 - g2*(*x2);
            //
            // double xc = (c2-c1)/(g1-g2);
            // double yc = g1*(xc) + c1;
            //
            // double cosTheta = acos(xn2*(*xn) + yn2*(*yn));
            // // double dist = sqrt(xc*xc + yc*yc);
            //
            // // normalVectorBetween(xn2, yn2, *x1, *x2, *y1, *y2);
            //
            // std::cout<<i<<' '<<cosTheta*180./pi<<std::endl;
            // xc and yc automatically assigned the coords

            // For no-normal flow at the edges, can we just set the normal vector components to zero?

            // if boundary face, *xn = 0.0, *yn=0.0;

            // if (cell_is_boundary(i) &&
            //     cell_is_boundary( node_friends(i, (j+1)%friend_num) ))
            // {
            //     *xn = 0.0;
            //     *yn = 0.0;
            // }

        }
    }

    return 1;
};

// Function to find the area of each node's control volume. THe function loops
// over each pair of centroids belonging to a node, while summing the area of
// each triangle to find the total area of the hexagon or pentagon.
int Mesh::CalcControlVolumeArea(void)
{
    int i, j, f, friend_num;
    double * x1, * y1, * x2, * y2, * xc, * yc, *t_area, area;

    for (i=0; i<node_num; i++)
    {

        f = node_friends(i,5);

        friend_num = 6;                                     // Assume hexagon (6 centroids)
        if (f == -1) {
            friend_num = 5;                                   // Check if pentagon (5 centroids)
        }

        xc = &node_pos_map(i,0,0);                          // set pointer to edge length array
        yc = &node_pos_map(i,0,1);

        t_area = &control_volume_surf_area_map(i);

        area = 0.0;
        for (j=0; j<friend_num; j++)                        // Loop through all centroids in the control volume
        {
            f = node_friends(i,j);

            x1 = &centroid_pos_map(i, j, 0);                  // set map coords for first centroid
            y1 = &centroid_pos_map(i, j, 1);

            x2 = &centroid_pos_map(i, (j+1)%friend_num, 0);   // set map coords for second centroid
            y2 = &centroid_pos_map(i, (j+1)%friend_num, 1);   // automatically loops around using %

            triangularArea(*t_area, *xc, *yc, *x1, *x2, *y1, *y2);     // calculate center coords between the two centroids.
            // xc and yc automatically assigned the coords
            area += *t_area;
        }
        control_volume_surf_area_map(i) = area;
    }

    return 1;
};


int Mesh::CalcControlVolumeMass(void)
{
  int i, j, k, as, ae, f, friend_num;
  double * lat1, * lon1, * lat2, * lon2, * lat3, * lon3, *t_area, r;
  double area, avg_area;
  double ax, ay, az;
  double bx, by, bz;
  double cx, cy, cz;
  double * vol;
  double vol1, vol2;

  // r = globals->radius.Value();//# + globals->shell_thickness.Value();
  // r_core = r - globals->h.Value();
  // r_ocean = r;
  // double mass_sum = 0.0;
  //
  // avg_area = 0.0;
  // for (i=0; i<node_num; i++)
  // {
  //
  //     f = node_friends(i,5);
  //
  //     friend_num = 6;                    // Assume hexagon (6 centroids)
  //     if (f == -1) {
  //         friend_num = 5;                // Check if pentagon (5 centroids)
  //     }
  //
  //     lat1 = &node_pos_sph(i,0);
  //     lon1 = &node_pos_sph(i,1);
  //
  //     vol = &control_volume_mass(i);     // set pointer to area of sub element
  //     *vol = 0.0;
  //     for (j=0; j<friend_num; j++)                       // Loop through all centroids in the control volume
  //     {
  //
  //         lat2 = &centroid_pos_sph(i,j%friend_num, 0);
  //         lon2 = &centroid_pos_sph(i,j%friend_num, 1);
  //
  //         lat3 = &centroid_pos_sph(i,(j+1)%friend_num, 0);
  //         lon3 = &centroid_pos_sph(i,(j+1)%friend_num, 1);
  //
  //         sph2cart(ax, ay, az, r_ocean, radConv*0.5 - *lat1, *lon1);
  //         sph2cart(bx, by, bz, r_ocean, radConv*0.5 - *lat2, *lon2);
  //         sph2cart(cx, cy, cz, r_ocean, radConv*0.5 - *lat3, *lon3);
  //
  //         volumeSphericalTriangle(vol1, ax, bx, cx, ay, by, cy, az, bz, cz, r_ocean);
  //
  //         sph2cart(ax, ay, az, r_core, radConv*0.5 - *lat1, *lon1);
  //         sph2cart(bx, by, bz, r_core, radConv*0.5 - *lat2, *lon2);
  //         sph2cart(cx, cy, cz, r_core, radConv*0.5 - *lat3, *lon3);
  //
  //         volumeSphericalTriangle(vol2, ax, bx, cx, ay, by, cy, az, bz, cz, r_core);
  //
  //         *vol += fabs(vol1)-fabs(vol2);
  //     }
  //     *vol *= 1000.0;
  // }
  //
  // for (i=0; i<node_num; i++)
  // {
  //   vol = &control_volume_mass(i);
  //   mass_sum += *vol;
  // }
  // std::cout<<"MASS: "<<mass_sum<<std::endl;
  // return 1;

  r = globals->radius.Value();//# + globals->shell_thickness.Value();

  double mass_sum = 0.0;

  avg_area = 0.0;
  for (i=0; i<node_num; i++)
  {

      f = node_friends(i,5);

      friend_num = 6;                    // Assume hexagon (6 centroids)
      if (f == -1) {
          friend_num = 5;                // Check if pentagon (5 centroids)
      }

      lat1 = &node_pos_sph(i,0);
      lon1 = &node_pos_sph(i,1);

      t_area = &control_volume_mass(i);
      *t_area = 0.0;
      for (j=0; j<friend_num; j++)                       // Loop through all centroids in the control volume
      {

          lat2 = &centroid_pos_sph(i,j%friend_num, 0);
          lon2 = &centroid_pos_sph(i,j%friend_num, 1);

          lat3 = &centroid_pos_sph(i,(j+1)%friend_num, 0);
          lon3 = &centroid_pos_sph(i,(j+1)%friend_num, 1);


          triangularAreaSph(area, *lat1, *lat2, *lat3, *lon1, *lon2, *lon3, r);     // calculate subelement area
          *t_area += area;

      }
      avg_area += *t_area;
  }
  avg_area /= node_num;
  for (i=0; i<node_num; i++)
  {
     t_area = &control_volume_mass(i);
    *t_area *= avg_area/(*t_area) * 1000.0 * globals->h.Value();

    mass_sum += *t_area;
  }

  return 1;
}

// Function to fine the areas of each subelement. Each subelement is a triangle
// subtended by the central node, one friend, and the centroid belonging to the
// friend. Values are stored in node_friend_element_areas_map
int Mesh::CalcElementAreas(void)
{
    int i, j, k, as, ae, f, f1, f2, friend_num;
    double * x1, * y1, * x2, * y2, * xc, * yc, *t_area;

    for (i=0; i<node_num; i++)
    {

        f = node_friends(i,5);

        friend_num = 6;                    // Assume hexagon (6 centroids)
        if (f == -1) {
            friend_num = 5;                // Check if pentagon (5 centroids)
        }
        for (j=0; j<friend_num; j++)                       // Loop through all centroids in the control volume
        {
            xc = &centroid_pos_map(i,j,0);                 // set pointer element centroid
            yc = &centroid_pos_map(i,j,1);

            double lat0 = centroid_pos_sph(i,j,0);
            double lon0 = centroid_pos_sph(i,j,1);
            f = node_friends(i,j);
            f1 = node_friends(i,j);
            f2 = node_friends(i,(j+1)%friend_num);

            double lat1 = node_pos_sph(i, 0);
            double lon1 = node_pos_sph(i, 1);

            double lat2 = node_pos_sph(f1, 0);
            double lon2 = node_pos_sph(f1, 1);
            double area= 0.0;
            double rr = globals->radius.Value();

            triangularAreaSph(area, lat0, lat1, lat2, lon0, lon1, lon2, rr);
            node_friend_element_areas_map(i, j, 0) = area;

            lat1 = node_pos_sph(f1, 0);
            lon1 = node_pos_sph(f1, 1);

            lat2 = node_pos_sph(f2, 0);
            lon2 = node_pos_sph(f2, 1);

            triangularAreaSph(area, lat0, lat1, lat2, lon0, lon1, lon2, rr);
            node_friend_element_areas_map(i, j, 1) = area;

            lat1 = node_pos_sph(f2, 0);
            lon1 = node_pos_sph(f2, 1);

            lat2 = node_pos_sph(i, 0);
            lon2 = node_pos_sph(i, 1);

            triangularAreaSph(area, lat0, lat1, lat2, lon0, lon1, lon2, rr);
            node_friend_element_areas_map(i, j, 2) = area;

            // if (i==137 && (f==1427))// || f==1453 || f ==1457))
            // {
            //   // double dist = 0.0;
            //   // std::cout<<i<<' '<<f<<' '<<*x1<<' '<<*y1<<' '<<*x2<<' '<<*y2<<std::endl;
            //   // distanceBetween(dist, *xc, *x1, *yc, *y1);
            //   std::cout<<i<<' '<<f<<' '<<area<<std::endl;
            //   // std::cout<<i<<' '<<nx<<' '<<ny<<std::endl;
            //   // std::cout<<i<<' '<<m<<std::endl;
            //   // std::cout<<i<<' '<<mesh->node_dists(i, (j+1)%friend_num)<<std::endl;
            //   // for (int k=0; k<3; k++) std::cout<<' '<<f1<<' '<<(*element_areas)(i,(j)%friend_num,k);
            //   // for (int k=0; k<3; k++) std::cout<<' '<<f2<<' '<<(*element_areas)(i,(j+1)%friend_num,k);
            //   // for (int k=0; k<3; k++) std::cout<<' '<<f2<<' '<<(*element_areas)(i,(j+2)%friend_num,k);
            // }

            for (k=0; k < 3; k++)
            {
                t_area = &node_friend_element_areas_map(i, j, k);     // set pointer to area of sub element


                // messy code TODO - try and tidy
                as = j+k;
                ae = (j+(k+1)%3)%(friend_num+2);
                if (k==0) as = 0;
                else if (k==2) {
                    ae = 0;
                }

                if (j == friend_num-1) {
                    if (k==1) ae = 1;
                    else if (k==2) as = 1;
                }

                x1 = &node_pos_map(i, as, 0);   // set map coords for first vertex of the element
                y1 = &node_pos_map(i, as, 1);

                x2 = &node_pos_map(i, ae, 0);   // set map coords for second vertex on the element
                y2 = &node_pos_map(i, ae, 1);

                // triangularArea(*t_area, *xc, *yc, *x1, *x2, *y1, *y2);     // calculate subelement area

                // double lat1 = node_pos_sph(i, 0);
                // double lon1 = node_pos_sph(i, 1);
                //
                // double lat2 = node_pos_sph(f, 0);
                // double lon2 = node_pos_sph(f, 1);
                // double area= 0.0;
                // double rr = globals->radius.Value();
                // triangularAreaSph(area, lat0, lat1, lat2, lon0, lon1, lon2, rr);

                // if (i==137 && (f==1427))// || f==1453 || f ==1457))
                // {
                //   double dist = 0.0;
                //   // std::cout<<i<<' '<<f<<' '<<*x1<<' '<<*y1<<' '<<*x2<<' '<<*y2<<std::endl;
                //   distanceBetween(dist, *xc, *x1, *yc, *y1);
                //   // std::cout<<i<<' '<<f<<' '<<*t_area<<std::endl;
                //   // std::cout<<i<<' '<<nx<<' '<<ny<<std::endl;
                //   // std::cout<<i<<' '<<m<<std::endl;
                //   // std::cout<<i<<' '<<mesh->node_dists(i, (j+1)%friend_num)<<std::endl;
                //   // for (int k=0; k<3; k++) std::cout<<' '<<f1<<' '<<(*element_areas)(i,(j)%friend_num,k);
                //   // for (int k=0; k<3; k++) std::cout<<' '<<f2<<' '<<(*element_areas)(i,(j+1)%friend_num,k);
                //   // for (int k=0; k<3; k++) std::cout<<' '<<f2<<' '<<(*element_areas)(i,(j+2)%friend_num,k);
                // }
                // if (i==1453 && (f==137))// || f==1457 || f ==1427))
                // {
                //   double dist = 0.0;
                //   distanceBetween(dist, *xc, *x1, *yc, *y1);
                //   // std::cout<<i<<' '<<f<<' '<<dist<<std::endl;
                //   std::cout<<i<<' '<<f<<' '<<*t_area<<std::endl;
                //   // std::cout<<i<<' '<<f<<' '<<as<<' '<<ae<<std::endl;
                //   // std::cout<<i<<' '<<f<<' '<<*x1<<' '<<*y1<<' '<<*x2<<' '<<*y2<<std::endl;
                // }

            }
        }
    }

    return 1;
};

// Function to evaluate trigonometric functions over latitude and longitude for
// every node.
int Mesh::CalcTrigFunctions(void)
{
    int i;
    double lat, lon;

    for (i=0; i<node_num; i++)
    {
        lat = node_pos_sph(i, 0);
        lon = node_pos_sph(i, 1);

        trigLat(i,0) = cos(lat);
        trigLat(i,1) = sin(lat);

        trigLon(i,0) = cos(lon);
        trigLon(i,1) = sin(lon);

        trig2Lat(i,0) = cos(2.0 * lat);
        trig2Lat(i,1) = sin(2.0 * lat);

        trig2Lon(i,0) = cos(2.0 * lon);
        trig2Lon(i,1) = sin(2.0 * lon);

        trigSqLat(i,0) = pow(cos(lat), 2.0);
        trigSqLat(i,1) = pow(sin(lat), 2.0);

        trigSqLon(i,0) = pow(cos(lon), 2.0);
        trigSqLon(i,1) = pow(sin(lon), 2.0);
    }

    return 1;
};

int Mesh::CalcLegendreFuncs(void)
{
    int i, l, m, l_max;
    double cosCoLat;

    Array2D<double> * temp_legendre;    // temp array for legendre polynomials
    Array2D<double> * temp_dlegendre;   // temp array for legendre polynomial derivs

    //-----------------------------

    l_max = globals->l_max.Value();

    temp_legendre = new Array2D<double>(2 * (l_max+1), 2 * (l_max+1));
    temp_dlegendre = new Array2D<double>(2 * (l_max+1), 2 * (l_max+1));

    Array3D<double> Pbar_lm(node_num, l_max+1, l_max+1);
    Array3D<double> Pbar_cosMLon(node_num, l_max+1, l_max+1);
    Array3D<double> Pbar_sinMLon(node_num, l_max+1, l_max+1);

    for (i=0; i<node_num; ++i)
    {
        cosCoLat = cos(pi*0.5 - node_pos_sph(i,0));     // cos(colatitude) of node i

        getLegendreFuncs(cosCoLat, *temp_legendre, l_max);
        getLegendreFuncsDeriv(cosCoLat, *temp_dlegendre, l_max);



        for (l = 0; l < l_max+1; l++)
        {
            for (m = 0; m <= l; m++)
            {

                // assign legendre pol at node i for degree l and order m
                Pbar_lm(i, l, m) = (*temp_legendre)(l, m);

                // assign legendre pol derivative at node i for degree l and order m
                // Pbar_lm_deriv(i, l, m) = (*temp_dlegendre)(l, m)*sin(pi*0.5 - node_pos_sph(i, 0));

                if ((l==2) && (m==0)) Pbar_20(i) = (*temp_legendre)(l, m);
                if ((l==2) && (m==2)) Pbar_22(i) = (*temp_legendre)(l, m);

            }
        }

        // calculate cos(m*longitude) and sin(m*longitude) for node i
        for (m=0; m < l_max+1; m++)
        {
            trigMLon(i, m, 0) = cos( (double)m * node_pos_sph( i, 1 ) );
            trigMLon(i, m, 1) = sin( (double)m * node_pos_sph( i, 1 ) );
        }

        // int count = 0;
        for (l = 0; l < l_max+1; l++)
        {
            for (m = 0; m <= l; m++)
            {
                Pbar_cosMLon(i, l, m) = Pbar_lm(i, l, m)*trigMLon(i, m, 0);
                Pbar_sinMLon(i, l, m) = Pbar_lm(i, l, m)*trigMLon(i, m, 1);

                // Pbar_deriv_cosMLon(i, l, m) = Pbar_lm_deriv(i, l, m)*trigMLon(i, m, 0);
                // Pbar_deriv_sinMLon(i, l, m) = Pbar_lm_deriv(i, l, m)*trigMLon(i, m, 1);

            }
        }

        int count = 0;
        for (l = 2; l < l_max+1; l++)
        {
            for (m = 0; m <= l; m++)
            {
                sh_matrix(i, count) = Pbar_cosMLon(i, l, m);
                if (globals->surface_type == FREE_LOADING) sh_matrix(i, count) *= globals->loading_factor[l];
                else if (globals->surface_type == LID_MEMBR ||
                         globals->surface_type == LID_LOVE) sh_matrix(i, count) *= globals->shell_factor_beta[l];


                count++;

                sh_matrix(i, count) = Pbar_sinMLon(i, l, m);
                if (globals->surface_type == FREE_LOADING) sh_matrix(i, count) *= globals->loading_factor[l];
                else if (globals->surface_type == LID_MEMBR ||
                         globals->surface_type == LID_LOVE) sh_matrix(i, count) *= globals->shell_factor_beta[l];

                count++;
            }
        }
    }

    int count = 0;
    for (m = 0; m < (l_max+1)*(l_max+2) - 6; m++)
    {
        for (i=0; i<node_num; i++)
        {
            sh_matrix_fort(count) = sh_matrix(i, m);
            // std::cout<<sh_matrix_fort(count)<<std::endl;
            count++;
        }
    }

    // free up memory!
    delete temp_legendre;
    delete temp_dlegendre;

    return 1;
}

int Mesh::CalcCoriolisOperatorCoeffs(void)
{
    int i, j, j2, j3, f, f_num;

    //--------------- Construct sparse matrix in CSR format --------------------

    //--------------- Define objects to construct matrix -----------------------

    // Here we assign the number of non-zero matrix elements for each operator.
    // There are 12 pentagons on the geodesic grid, each with 5 neighbours and
    // one central node (total of 6). The other cells are hexagons with 6
    // neighbours and one central node (total of 7). Additionally, there is a
    // row for each variable (total 3). Therefore, the total number
    // of non-zero coefficients in our sparse matrix is given as:
    int nNonZeroGrad = (12*5)*9 + (face_num - 12*5)*10;

    operatorCoriolis = SpMat(face_num, face_num);

    // Allocate memory for non-zero entries
    operatorCoriolis.reserve(nNonZeroGrad);

    std::vector<Triplet> coefficients;
    coefficients.reserve(nNonZeroGrad);


    int     * colIndx;    // column index for each non-zero element
    int     * rowIndxX;    // row indices for CSR matrix format (x component)
    // int     * rowIndxY;    // row indices for CSR matrix format (y component)
    double * nzCoeffsX;   // non-zero grad operator coeffis in x direction

    int * rowStartX;      // start row indices (x component)
    int * rowEndX;        // start row indices (x component)


    int         error;    // error message from MKL

    // // Define c (0) based indexing
    // sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;

    // // Define all operations the original (non-transpose) matrix
    // sparse_operation_t operation = SPARSE_OPERATION_NON_TRANSPOSE;

    //---------------- Initialise objects into memory --------------------------

    rowStartX = new int[face_num];
    rowEndX   = new int[face_num];
    rowIndxX  = new int[face_num + 1];


    nzCoeffsX = new double[nNonZeroGrad];
    colIndx   = new int[nNonZeroGrad];

    // assign the non-zero coefficients and their column indexes in CSR format
    int count=0;

    for (int i=0; i<face_num; ++i)
    {
        double v_tang = 0.0;
        int friend_num = 10;

        int n1, n2;
        n1 = face_nodes(i,0);
        n2 = face_nodes(i,1);
        if (node_friends(n1,5) < 0) friend_num--;
        if (node_friends(n2,5) < 0) friend_num--;

        for (int j=0; j<friend_num; j++) {
            int f_ID = face_interp_friends(i, j);
            colIndx[count] = f_ID;
            nzCoeffsX[count] = -2.0 * globals->angVel.Value() * sin(face_intercept_pos_sph(i, 0))*face_interp_weights(i,j)*face_len(f_ID)/face_node_dist(i);

            count++;

            double coeff = -2.0 * globals->angVel.Value() * sin(face_intercept_pos_sph(i, 0))*face_interp_weights(i,j)*face_len(f_ID)/face_node_dist(i);
            // double coeff = -2.0 * globals->angVel.Value() * sin(face_intercept_pos_sph(f_ID, 0)) * face_interp_weights(i, j) * face_len(f_ID) / face_node_dist(i);
            coefficients.push_back( Triplet(i, f_ID, coeff) );

            // v_tang += v_tm1(f_ID)*grid->face_interp_weights(i,j)*grid->face_len(f_ID);
        }

        // if (i<face_num-1) rowIndxX[i+1] = count;
        // v_tang /= grid->face_node_dist(i);

        // dvdt(i) -= 2.0 * omega * sin(grid->face_intercept_pos_sph(i, 0)) * v_tang;
        //
        // dvdt(i) -= v_tm1(i)*globals->alpha.Value();
    }

    // Fill the sparse operator with the list of triplets
    operatorCoriolis.setFromTriplets(coefficients.begin(), coefficients.end());
    operatorCoriolis.makeCompressed();

    return 1;
}

int Mesh::CalcLinearDragOperatorCoeffs(void)
{
    int i, j, j2, j3, f, f_num;

    //--------------- Construct sparse matrix in CSR format --------------------

    //--------------- Define objects to construct matrix -----------------------

    // Here we assign the number of non-zero matrix elements for each operator.
    // There are 12 pentagons on the geodesic grid, each with 5 neighbours and
    // one central node (total of 6). The other cells are hexagons with 6
    // neighbours and one central node (total of 7). Additionally, there is a
    // row for each variable (total 3). Therefore, the total number
    // of non-zero coefficients in our sparse matrix is given as:
    int nNonZeroGrad = face_num;

    // Define shape of operator
    operatorLinearDrag = SpMat(face_num, face_num);

    // Allocate memory for non-zero entries
    operatorLinearDrag.reserve(face_num);

    int     * colIndx;    // column index for each non-zero element
    int     * rowIndxX;    // row indices for CSR matrix format (x component)
    int     * rowIndxY;    // row indices for CSR matrix format (y component)
    double * nzCoeffsX;   // non-zero grad operator coeffis in x direction

    int * rowStartX;      // start row indices (x component)
    int * rowEndX;        // start row indices (x component)


    int         error;    // error message from MKL

    // Define c (0) based indexing
    // sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;

    // // Define all operations the original (non-transpose) matrix
    // sparse_operation_t operation = SPARSE_OPERATION_NON_TRANSPOSE;

    //---------------- Initialise objects into memory --------------------------

    rowStartX = new int[face_num];
    rowEndX   = new int[face_num];
    rowIndxX  = new int[face_num + 1];


    nzCoeffsX = new double[nNonZeroGrad];
    colIndx   = new int[nNonZeroGrad];

    // assign the non-zero coefficients and their column indexes in CSR format
    int count=0;


    Triplet T;
    for (i=0; i<face_num; i++)
    {

        //   nzCoeffsX[count] = -globals->alpha.Value();
        //   colIndx[count] = i;

          operatorLinearDrag.insert(i,i) = -globals->alpha.Value();

        //   count++;

          // nzCoeffsX[count] = -globals->alpha.Value();
          // colIndx[count] = i*3+1;
          // count++;
    }

    // assign the row indexes for the sparse matrix in CSR format
    // rowIndxX[0] = 0;         // first element is always zero

    // for (i=1; i<face_num; i++) {
    //   rowIndxX[i] = i;
    // }

    // // last element is always the number of non-zero coefficients
    // rowIndxX[face_num] = nNonZeroGrad;

    // for (i=0; i<face_num; i++) {
    //   rowStartX[i] = rowIndxX[i];
    //   rowEndX[i] = rowIndxX[i+1]; }

    // Create the operator for the x and y components
    // operatorLinearDrag = new sparse_matrix_t;
    // // sparse_matrix_t * operatorLinearDragTemp = new sparse_matrix_t;
    // error = mkl_sparse_d_create_csr(operatorLinearDrag, indexing, face_num, face_num, rowStartX, rowEndX, colIndx, nzCoeffsX);

    // matrix_descr descrp;
    // descrp.type = SPARSE_MATRIX_TYPE_GENERAL;
    // error = mkl_sparse_copy (*operatorLinearDragTemp, descrp, operatorLinearDrag);
    //
    // delete[] rowStartX;
    // delete[] rowEndX;
    // delete[] rowIndxX;
    // delete[] colIndx;
    // delete[] nzCoeffsX;
    // delete operatorLinearDragTemp;

    return 1;
}

int Mesh::CalcGradOperatorCoeffs(void)
{
    int i, j, j2, j3, f, f_num;

    //--------------- Construct sparse matrix in CSR format --------------------

    //--------------- Define objects to construct matrix -----------------------

    // Here we assign the number of non-zero matrix elements for each operator.
    // There are 12 pentagons on the geodesic grid, each with 5 neighbours and
    // one central node (total of 6). The other cells are hexagons with 6
    // neighbours and one central node (total of 7). Additionally, there is a
    // row for each variable (total 3). Therefore, the total number
    // of non-zero coefficients in our sparse matrix is given as:
    int nNonZeroGrad = 2*face_num;//(node_num-12)*7 + 12*6;

    operatorGradient = SpMat(face_num, node_num);

    // Allocate memory for non-zero entries
    operatorGradient.reserve(nNonZeroGrad);

    std::vector<Triplet> coefficients;
    coefficients.reserve(nNonZeroGrad);

    int     * colIndx;    // column index for each non-zero element
    int     * rowIndxX;    // row indices for CSR matrix format (x component) component)
    double * nzCoeffsX;   // non-zero grad operator coeffis in x direction

    int * rowStartX;      // start row indices (x component)
    int * rowEndX;        // start row indices (x component)


    int         error;    // error message from MKL

    // Define c (0) based indexing
    // sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;

    // // Define all operations the original (non-transpose) matrix
    // sparse_operation_t operation = SPARSE_OPERATION_NON_TRANSPOSE;

    //---------------- Initialise objects into memory --------------------------

    rowStartX = new int[face_num];
    rowEndX   = new int[face_num];
    rowIndxX  = new int[face_num + 1];


    nzCoeffsX = new double[nNonZeroGrad];
    colIndx   = new int[nNonZeroGrad];

    // assign the non-zero coefficients and their column indexes in CSR format
    int count=0;

    double coeff;
    for (i=0; i<face_num; i++)
    {
        int inner_node_ID, outer_node_ID;
        inner_node_ID = face_nodes(i, 0);
        outer_node_ID = face_nodes(i, 1);

        double dist;
        dist = face_node_dist(i);

        colIndx[count] = inner_node_ID;
        nzCoeffsX[count] = -globals->g.Value()*(-face_centre_m(i,0))/dist;
        count++;

        colIndx[count] = outer_node_ID;
        nzCoeffsX[count] = -globals->g.Value()*(face_centre_m(i,1))/dist;
        count++;

        coeff = (-face_centre_m(i, 0)) / dist;
        coefficients.push_back(Triplet(i, inner_node_ID, coeff));

        coeff = (face_centre_m(i, 1)) / dist;
        coefficients.push_back(Triplet(i, outer_node_ID, coeff));
    }

    operatorGradient.setFromTriplets(coefficients.begin(), coefficients.end());
    operatorGradient.makeCompressed();

    rowIndxX[0] = 0;         // first element is always zero
    for (i=1; i<face_num; i++) {
      rowIndxX[i] = 2*i;
    }

    // last element is always the number of non-zero coefficients
    rowIndxX[face_num] = nNonZeroGrad;

    for (i=0; i<face_num; i++) {
      rowStartX[i] = rowIndxX[i];
      rowEndX[i] = rowIndxX[i+1]; }

    // Create the operator for the x and y components
    // operatorGradient = new sparse_matrix_t;
    // error = mkl_sparse_d_create_csr(operatorGradient, indexing, face_num, node_num, rowStartX, rowEndX, colIndx, nzCoeffsX);

    // matrix_descr descrp;
    // descrp.type = SPARSE_MATRIX_TYPE_GENERAL;


    // // error = mkl_sparse_set_dotmv_hint (*operatorMomentum, operation, descrp, 1000000000);
    // error = mkl_sparse_set_mv_hint (*operatorGradient, operation, descrp, 100000000);
    // error = mkl_sparse_set_memory_hint (*operatorGradient, SPARSE_MEMORY_AGGRESSIVE);
    // error = mkl_sparse_optimize (*operatorGradient);

    // delete operatorGradientX;

    // delete[] rowStartX;
    // delete[] rowEndX;
    // delete[] rowIndxX;
    // // delete[] nzCoeffsX;
    // delete[] colIndx;

    return 1;
}

int Mesh::CalcDivOperatorCoeffs(void)
{
    int i, j, j2, j3, f, f_num;

    // why is this line here? Where would it be better placed?

    //--------------- Construct sparse matrix in CSR format --------------------

    //--------------- Define objects to construct matrix -----------------------

    // Here we assign the number of non-zero matrix elements for each operator.
    // There are 12 pentagons on the geodesic grid, each with 5 faces.
    // The other cells are hexagons with 6 faces. Therefore, the total
    // number of non-zero coefficients in our sparse matrix is given as:
    int nNonZeroGrad = ((node_num-12)*6 + 12*5);

    operatorDivergence = SpMat(node_num, face_num);

    // Allocate memory for non-zero entries
    operatorDivergence.reserve(nNonZeroGrad);

    std::vector<Triplet> coefficients;
    coefficients.reserve(nNonZeroGrad);

    int     * colIndxX;    // column index for each non-zero element
    int     * rowIndxX;    // row indices for CSR matrix format (x component) component)
    double * nzCoeffsX;   // non-zero grad operator coeffis in x direction

    int * rowStartX;      // start row indices (x component)
    int * rowEndX;        // start row indices (x component)

    int         error;    // error message from MKL

    // Define c (0) based indexing
    // sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;

    // // Define all operations the original (non-transpose) matrix
    // sparse_operation_t operation = SPARSE_OPERATION_NON_TRANSPOSE;

    //---------------- Initialise objects into memory --------------------------

    rowStartX = new int[node_num];
    rowEndX   = new int[node_num];
    rowIndxX  = new int[node_num + 1];

    nzCoeffsX = new double[nNonZeroGrad];
    colIndxX   = new int[nNonZeroGrad];

    // assign the non-zero coefficients and their column indexes in CSR format
    int count=0;

    for (i=0; i<node_num; i++)
    {
        f_num = 6;
        if (node_friends(i, 5) == -1) f_num = 5;

        for (j=0; j<f_num; j++)
        {
            int face_id = faces(i,j);

            // get edge length of current edge
            double edge_len = face_len(face_id);
            int dir = node_face_dir(i, j);
            double area = control_volume_surf_area_map(i);

            nzCoeffsX[count]  =  -dir*edge_len/area;
            colIndxX[count] = face_id;
            count++;

            double coeff = -dir * edge_len / area;
            coefficients.push_back(Triplet(i, face_id, coeff));
        }

        if (i<node_num-1) rowIndxX[i+1] = count;
    }

    operatorDivergence.setFromTriplets(coefficients.begin(), coefficients.end());
    operatorDivergence.makeCompressed();

    // assign the row indexes for the sparse matrix in CSR format
    rowIndxX[0] = 0;         // first element is always zero
    rowIndxX[node_num] = nNonZeroGrad;

    // last element is always the number of non-zero coefficients
    // rowIndxY[node_num] = nNonZeroGrad;

    for (i=0; i<node_num; i++) {
      rowStartX[i] = rowIndxX[i];
      rowEndX[i] = rowIndxX[i+1]; }

    // Create the operator for the x and y components
    // operatorDivergenceX = new sparse_matrix_t;
    // error = mkl_sparse_d_create_csr(operatorDivergenceX, indexing, node_num, face_num, rowStartX, rowEndX, colIndxX, nzCoeffsX);

    // matrix_descr descrp;
    // descrp.type = SPARSE_MATRIX_TYPE_GENERAL;

    // error = mkl_sparse_set_mv_hint (*operatorDivergenceX, operation, descrp, 1000000000);
    // error = mkl_sparse_set_memory_hint (*operatorDivergenceX, SPARSE_MEMORY_AGGRESSIVE);
    // error = mkl_sparse_optimize (*operatorDivergenceX);

    // delete operatorDivergenceX;
    //
    // delete[] rowStartX;
    // delete[] rowEndX;
    // delete[] rowIndxX;
    // delete[] nzCoeffsX;
    // delete[] colIndxX;

    return 1;
}

// // Function to generate a sparse matrix A for the Laplacian operator, convert it
// // to CSR format, and create a solver object that is capable of solving the
// // linear system A*p = d to find the pressure field, p, give a rhs vector d.
// int Mesh::CalcLaplaceOperatorCoeffs(void)
// {
//     int i, j, j2, j3, f, f_num;
//     // Now create the gradient operators using the coefficients above

//     struct node {
//       int ID;
//       int F_NUM;
//     };


//     //--------------- Construct sparse matrix in CSR format --------------------

//     //--------------- Define objects to construct matrix -----------------------

//     double * nzCoeffs_grad_x;   // non-zero grad operator coeffis in x direction
//     double * nzCoeffs_grad_y;   // non-zero grad operator coeffis in y direction
//     double * nzCoeffs_div_x;    // non-zero div operator coeffis in x direction
//     double * nzCoeffs_div_y;    // non-zero div operator coeffis in y direction

//     // Here we assign the number of non-zero matrix elements for each operator.
//     // There are 12 pentagons on the geodesic grid, each with 5 neighbours and
//     // one central node (total of 6). The other cells are hexagons with 6
//     // neighbours and one central node (total of 7). Therefore, the total number
//     // of non-zero coefficients in our sparse matrix is given as:
//     int nNonZero = (node_num-12)*7 + 12*6;

//     int     * colIndx;    // column index for each non-zero element
//     int     * rowIndx;    // row indices for CSR matrix format
//     double * nzCoeffs;    // array for each non-zero coefficient

//     int         error;    // error message from MKL

//     int * rowStart, * rowEnd; // row indexes for CSR sparse matrix format


//     //---------------- Initialise objects into memory --------------------------

//     rowStart = new int[node_num];
//     rowEnd   = new int[node_num];

//     colIndx  = new int[nNonZero];
//     nzCoeffs = new double[nNonZero];
//     rowIndx  = new int[node_num + 1];

//     nzCoeffs_grad_x = new double[nNonZero];
//     nzCoeffs_grad_y = new double[nNonZero];
//     nzCoeffs_div_x  = new double[nNonZero];
//     nzCoeffs_div_y  = new double[nNonZero];

//     // assign the non-zero coefficients and their column indexes in CSR format
//     int count=0;
//     for (i=0; i<node_num; i++)
//     {
//         f_num = 6;               // Assume hexagon (6 neighbours)
//         f = node_friends(i,5);
//         if (f == -1) f_num = 5;  // Check if pentagon (5 neighbours)


//         std::vector<node> cols(f_num+1);
//         cols[0].ID = i;
//         cols[0].F_NUM = 0;

//         // Add each friend to the node list
//         for (j=0; j<f_num; j++) {
//             cols[j+1].ID = node_friends(i, j);
//             cols[j+1].F_NUM = j+1; }

//         // Order the node list so that lowest node ID's come first
//         // This ensures each friend corresponds to its correct column
//         std::sort(cols.begin(), cols.end(),
//           [](auto const &a, auto const &b) { return a.ID < b.ID; });

//         // Using the ordered list to construct row i coefficients in the operator
//         double cos_a, sin_a;
//         for (j=0; j<f_num+1; j++) {
//           cos_a = node_vel_trans(i, cols[j].F_NUM, 0);
//           sin_a = node_vel_trans(i, cols[j].F_NUM, 1);

//           nzCoeffs_grad_x[count] = grad_coeffs(i, cols[j].F_NUM, 0);
//           nzCoeffs_grad_y[count] = grad_coeffs(i, cols[j].F_NUM, 1);

//           // Because the gradient operator returns a vector, we have to transform
//           // the vector into mapped coordinates. We do this here by premultiplying
//           // the components of the divergence operator by the transpose of the
//           // the transform matrix defined in node_vel_trans.
//           nzCoeffs_div_x[count]  =  1e6*(cos_a*div_coeffs(i, cols[j].F_NUM, 0) - sin_a*div_coeffs(i,  cols[j].F_NUM, 1));
//           nzCoeffs_div_y[count]  =  1e6*(sin_a*div_coeffs(i, cols[j].F_NUM, 0) + cos_a*div_coeffs(i,  cols[j].F_NUM, 1));

//           // Assign the column index. These operators act on each
//           // velocity component, which occupies every 3rd index, starting at
//           // column index 2
//           colIndx[count] = cols[j].ID;
//           // colIndxY[count] = cols[j].ID*3+1;
//           count++; }



//     }



//     // assign the row indexes for the sparse matrix in CSR format

//     rowIndx[0] = 0;         // first element is always zero

//     for (i=0; i<node_num; i++) {
//       if (node_friends(i,5) == -1) {
//         rowIndx[i+1] = rowIndx[i]+6;
//       }
//       else {
//         rowIndx[i+1] = rowIndx[i]+7;
//       }
//     }

//     // last element is always the number of non-zero coefficients
//     rowIndx[node_num] = nNonZero;

//     for (i=0; i<node_num; i++) {
//       rowStart[i] = rowIndx[i];
//       rowEnd[i] = rowIndx[i+1]; }


//     //----------------- Set up sparse matrix for div and grad ------------------


//     // define handles to the sparse matrices for each component of the grad and
//     // div operators. The Laplacian matrix will be constructed from these.
//     sparse_matrix_t * DIV_X, * DIV_Y, * GRAD_X, * GRAD_Y;
//     sparse_matrix_t * LAP_X, * LAP_Y;

//     DIV_X   = new sparse_matrix_t;
//     DIV_Y   = new sparse_matrix_t;
//     GRAD_X  = new sparse_matrix_t;
//     GRAD_Y  = new sparse_matrix_t;
//     LAP_X   = new sparse_matrix_t;
//     LAP_Y   = new sparse_matrix_t;



//     // Define c (0) based indexing
//     sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;

//     // Define all operations the original (non-transpose) matrix
//     sparse_operation_t operation = SPARSE_OPERATION_NON_TRANSPOSE;

//     // Assign non-zero values to matrices for the components of div and grad operators
//     error = mkl_sparse_d_create_csr(DIV_X,  indexing, node_num, node_num, rowStart, rowEnd, colIndx, nzCoeffs_div_x);
//     error = mkl_sparse_d_create_csr(DIV_Y,  indexing, node_num, node_num, rowStart, rowEnd, colIndx, nzCoeffs_div_y);

//     error = mkl_sparse_d_create_csr(GRAD_X, indexing, node_num, node_num, rowStart, rowEnd, colIndx, nzCoeffs_grad_x);
//     error = mkl_sparse_d_create_csr(GRAD_Y, indexing, node_num, node_num, rowStart, rowEnd, colIndx, nzCoeffs_grad_y);

//     // Compute the multiplication of each component of the div dot grad
//     error = mkl_sparse_spmm(operation, *DIV_X, *GRAD_X, LAP_X); // div_x * grad_x
//     error = mkl_sparse_spmm(operation, *DIV_Y, *GRAD_Y, LAP_Y); // div_y * grad_y

//     // Sum each component to give the total Grad^2 operator (Laplacian)
//     // This is equivalent to (div_x * grad_x) + (div_y + grad_y)
//     sparse_matrix_t * LAP = new sparse_matrix_t;
//     error = mkl_sparse_d_add(operation, *LAP_X, 1.0, *LAP_Y, LAP);



//     delete DIV_X;
//     delete DIV_Y;
//     delete GRAD_X;
//     delete GRAD_Y;
//     delete LAP_X;
//     delete LAP_Y;

//     delete[] rowStart;
//     delete[] rowEnd;
//     delete[] colIndx;
//     delete[] nzCoeffs_div_x;
//     delete[] nzCoeffs_div_y;
//     delete[] nzCoeffs_grad_x;
//     delete[] nzCoeffs_grad_y;


//     int    * rowStartLap;
//     int      * rowEndLap;
//     int     * colIndxLap;
//     double * nzCoeffsLap;

//     int nc, nr; //nrows, ncols
//     int nNonZeroLap;

//     rowStartLap = NULL;
//     rowEndLap   = NULL;
//     colIndxLap  = NULL;
//     nzCoeffsLap = NULL;

//     // Now extract the the row and column index arrays, and nonzero coefficients
//     // from the Laplacian operator. These will be used to set up the solver for
//     // the problem L*p = d below, where L is the Laplacian operator and is to be
//     // inverted to find the solution p.
//     error = mkl_sparse_d_export_csr (*LAP, &indexing, &nc, &nr, &rowStartLap, &rowEndLap, &colIndxLap, &nzCoeffsLap);

//     nNonZeroLap = rowEndLap[nr - 1];

//     int * colIndxU = new int[nNonZeroLap];
//     int * colIndxV = new int[nNonZeroLap];
//     int * rowStartU = new int[2*node_num];
//     int * rowStartV = new int[2*node_num];
//     int * rowEndU = new int[2*node_num];
//     int * rowEndV = new int[2*node_num];

//     double * nzCoeffsLap2 = new double[nNonZeroLap];

//     for (i=0; i<nNonZeroLap; i++)
//     {
//       colIndxU[i] = 3*colIndxLap[i];
//       colIndxV[i] = 3*colIndxLap[i]+1;

//       nzCoeffsLap2[i] = nzCoeffsLap[i];
//     }

//     rowStartU[0] = 0.0;
//     rowStartV[0] = 0.0;
//     for (i=1; i<node_num; i++)
//     {
//       rowStartU[2*i-1] = rowStartLap[i];
//       rowStartU[2*i] = rowStartLap[i];}

//     for (i=0; i<node_num; i++)
//     {
//       rowStartV[2*i] = rowStartLap[i];
//       rowStartV[2*i+1] = rowStartLap[i];

//     }
//     rowStartU[2*node_num-1] = nNonZeroLap;
//     for (i=0; i<2*node_num-1; i++)
//     {
//       rowEndU[i] = rowStartU[i+1];
//       rowEndV[i] = rowStartV[i+1];
//     }


//     rowEndU[2*node_num-1] = nNonZeroLap;
//     rowEndV[2*node_num-1] = nNonZeroLap;

//     // assign the row indexes for the sparse matrix in CSR format
//     // for (i=0; i<2*node_num; i++)
//     // {
//     //   std::cout<<i<<' '<<rowStartV[i]<<' '<<rowEndV[i]<<' '<<nNonZeroLap<<std::endl;
//     // }


//     sparse_matrix_t * LAP2_X = new sparse_matrix_t;
//     sparse_matrix_t * LAP2_Y = new sparse_matrix_t;
//     error = mkl_sparse_d_create_csr(LAP2_X,  indexing, 2*node_num, 3*node_num, rowStartU, rowEndU, colIndxU, nzCoeffsLap2);
//     error = mkl_sparse_d_create_csr(LAP2_Y,  indexing, 2*node_num, 3*node_num, rowStartV, rowEndV, colIndxV, nzCoeffsLap2);

//     // Assing memory to the handle for the Laplacian operator matrix
//     operatorLaplacian = new sparse_matrix_t;  // defined in mesh.h
//     error = mkl_sparse_d_add(operation, *LAP2_X, 1.0, *LAP2_Y, operatorLaplacian);


//     delete[] rowStartU;
//     delete[] rowEndU;
//     delete[] rowStartV;
//     delete[] rowEndV;
//     delete[] colIndxU;
//     delete[] colIndxV;
//     // delete[] nzCoeffsLap;
//     delete[] nzCoeffsLap2;

//     delete LAP2_X;
//     delete LAP2_Y;

//     return 1;
// }



// int Mesh::GenerateMomentumOperator(void)
// {
//     int error;
//     // Define c (0) based indexing
//     sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;

//     // Define all operations the original (non-transpose) matrix
//     sparse_operation_t operation = SPARSE_OPERATION_NON_TRANSPOSE;

//     operatorMomentum = new sparse_matrix_t;
//     // int error = mkl_sparse_d_add(operation, *operatorLaplacian, 1.0, *operatorGradient, operatorMomentum);
//     // int error = mkl_sparse_d_add(operation, *operatorCoriolis, 1.0, *operatorGradient, operatorMomentum);
//     // error = mkl_sparse_d_add(operation, *operatorMomentum, 1.0, *operatorCoriolis, operatorMomentum);

//     if (globals->fric_type == LINEAR) mkl_sparse_d_add(operation, *operatorCoriolis, 1.0, *operatorLinearDrag, operatorMomentum);
//     else mkl_sparse_d_add(operation, *operatorLinearDrag, 0.0, *operatorCoriolis, operatorMomentum);

//     mkl_sparse_order (*operatorMomentum);


//     matrix_descr descrp;
//     descrp.type = SPARSE_MATRIX_TYPE_GENERAL;

//     int iter_num = 9999999;
//     // error = mkl_sparse_set_dotmv_hint (*operatorMomentum, operation, descrp, 1000000000);
//     error = mkl_sparse_set_mv_hint (*operatorMomentum, operation, descrp, iter_num);
//     // error = mkl_sparse_set_mm_hint (*operatorMomentum, operation, descrp, SPARSE_LAYOUT_COLUMN_MAJOR, 1, iter_num);
//     error = mkl_sparse_set_memory_hint (*operatorMomentum, SPARSE_MEMORY_AGGRESSIVE);
//     error = mkl_sparse_optimize (*operatorMomentum);
//     error = mkl_sparse_optimize (*operatorMomentum);



//     // op_mom = new SpMat(2*node_num,3*node_num);
//     //
//     // int    * rowStartLap;
//     // int      * rowEndLap;
//     // int     * colIndxLap;
//     // double * nzCoeffsLap;
//     //
//     // int nc, nr; //nrows, ncols
//     // int nNonZeroLap;
//     //
//     // rowStartLap = NULL;
//     // rowEndLap   = NULL;
//     // colIndxLap  = NULL;
//     // nzCoeffsLap = NULL;
//     //
//     // // Now extract the the row and column index arrays, and nonzero coefficients
//     // // from the Laplacian operator. These will be used to set up the solver for
//     // // the problem L*p = d below, where L is the Laplacian operator and is to be
//     // // inverted to find the solution p.
//     // error = mkl_sparse_d_export_csr (*operatorMomentum, &indexing, &nc, &nr, &rowStartLap, &rowEndLap, &colIndxLap, &nzCoeffsLap);
//     //
//     // nNonZeroLap = rowEndLap[nr - 1];
//     //
//     // // for (i=0; i<2*node_num; i++)
//     // // {
//     // //   rowIndx[i] = rowStartLap[i];
//     // // }
//     // // rowIndx[node_num] = nNonZeroLap;
//     //
//     // int count = 0;
//     // int num_row;
//     // for (int i=0; i<2*node_num; ++i)
//     // {
//     //   num_row = rowEndLap[i]-rowStartLap[i];
//     //
//     //   for (int j=0; j<num_row; ++j)
//     //   {
//     //
//     //     (*op_mom).insert(i, colIndxLap[rowEndLap[i]+j]) = nzCoeffsLap[count++];
//     //
//     //   }
//     // }
//     //
//     // (*op_mom).makeCompressed();

//     // Free space!
//     // delete operatorLaplacian;
//     // delete operatorGradient;
//     delete operatorCoriolis;
//     delete operatorLinearDrag;

//     return 1;
// }
// Function to generate a sparse matrix A for the Laplacian operator, convert it
// to CSR format, and create a solver object that is capable of solving the
// linear system A*p = d to find the pressure field, p, give a rhs vector d.
// int Mesh::GeneratePressureSolver(void)
// {
//     int i, j, j2, j3, f, f_num;
//     double areaCV, *coeff;
//     double D;


//     double * hVar = new double[node_num];
//     double a_ratio = 0.0;
//     for (i=0; i<node_num; i++) hVar[i] = 1.0;//globals->h.Value()*(1. + a_ratio*trigLat(i,0)*pow(trigMLon(i,4,0),4.0));


//     // }

//     Array2D<double> * pressure_matrix;
//     pressure_matrix = new Array2D<double>(node_num, 7);

//     // why is this line here? Where would it be better placed?

//     //-------------- Calculate the laplacian matrix coefficients ---------------

//     for (i=0; i<node_num; i++)
//     {

//         f_num = 6;                  // Assume hexagon (6 centroids)
//         f = node_friends(i,5);
//         if (f == -1) f_num = 5;     // Check if pentagon (5 centroids)

//         // Control volume area
//         areaCV = control_volume_surf_area_map(i);

//         // Calculate the central node coefficient first
//         coeff = &(*pressure_matrix)(i, 0);
//         *coeff = 0.0;
//         for (j=0; j<f_num; j++)
//         {
//             j2 = (j-1)%f_num;
//             if (j2 < 0) j2 += f_num;
//             f = node_friends(i, j);

//             //          this part is the spatially
//             //          varying scalar, TODO change
//             //          to ocean thickness!!!!!
//             //        -----------------------------
//             D =  0.5*(hVar[i] + hVar[f]) * control_vol_edge_len(i, j2) /  ( control_vol_edge_centre_m(i, j2) * node_dists(i, j) );
//             // 0.5*(trigLon(i, 1) + trigLon(f, 1)) *
//             *coeff += D;
//         }
//         *coeff *= -1./areaCV; // Need the negative sign here!

//         // Now calculate coeffs for neighbouring nodes (matrix diagonals)
//         for (j=0; j<f_num; j++)
//         {
//             coeff = &(*pressure_matrix)(i, j+1);

//             j2 = (j-1)%f_num;
//             if (j2 < 0) j2 += f_num;
//             f = node_friends(i,j);


//             D = 0.5*(hVar[i] + hVar[f]) * control_vol_edge_len(i, j2) /  ( control_vol_edge_centre_m(i, j2) * node_dists(i, j) );
//             // 0.5*(trigLon(i, 1) + trigLon(f, 1)) *

//             *coeff = D;
//             *coeff *= 1./areaCV; // No negative sign here!
//         }
//     }


//     //--------------- Construct sparse matrix in CSR format --------------------

//     //--------------- Define objects to construct matrix -----------------------

//     double * nzCoeffs_grad_x;   // non-zero grad operator coeffis in x direction
//     double * nzCoeffs_grad_y;   // non-zero grad operator coeffis in y direction
//     double * nzCoeffs_div_x;    // non-zero div operator coeffis in x direction
//     double * nzCoeffs_div_y;    // non-zero div operator coeffis in y direction

//     // Here we assign the number of non-zero matrix elements for each operator.
//     // There are 12 pentagons on the geodesic grid, each with 5 neighbours and
//     // one central node (total of 6). The other cells are hexagons with 6
//     // neighbours and one central node (total of 7). Therefore, the total number
//     // of non-zero coefficients in our sparse matrix is given as:
//     int nNonZero = (node_num-12)*7 + 12*6;

//     int     * colIndx;    // column index for each non-zero element
//     int     * rowIndx;    // row indices for CSR matrix format
//     double * nzCoeffs;    // array for each non-zero coefficient

//     int         error;    // error message from MKL

//     int * rowStart, * rowEnd; // row indexes for CSR sparse matrix format


//     //---------------- Initialise objects into memory --------------------------

//     rowStart = new int[node_num];
//     rowEnd   = new int[node_num];

//     colIndx  = new int[nNonZero];
//     nzCoeffs = new double[nNonZero];
//     rowIndx  = new int[node_num + 1];

//     nzCoeffs_grad_x = new double[nNonZero];
//     nzCoeffs_grad_y = new double[nNonZero];
//     nzCoeffs_div_x  = new double[nNonZero];
//     nzCoeffs_div_y  = new double[nNonZero];

//     // assign the non-zero coefficients and their column indexes in CSR format
//     int count=0;
//     for (i=0; i<node_num; i++)
//     {
//         f_num = 6;               // Assume hexagon (6 neighbours)
//         f = node_friends(i,5);
//         if (f == -1) f_num = 5;  // Check if pentagon (5 neighbours)

//         nzCoeffs[count] = (*pressure_matrix)(i, 0);  // assign diagonal element

//         nzCoeffs_grad_x[count] = grad_coeffs(i, 0, 0);
//         nzCoeffs_grad_y[count] = grad_coeffs(i, 0, 1);

//         nzCoeffs_div_x[count]  = div_coeffs(i, 0, 0);
//         nzCoeffs_div_y[count]  = div_coeffs(i, 0, 1);
//         // nzCoeffs_div_x[count]  = hVar[i] * div_coeffs(i, 0, 0);
//         // nzCoeffs_div_y[count]  = hVar[i] * div_coeffs(i, 0, 1);

//         colIndx[count] = i;

//         count++;
//         double cos_a, sin_a;
//         for (j=0; j<f_num; j++)                       // assign off-diagonals
//         {
//             nzCoeffs[count] = (*pressure_matrix)(i, j+1);
//             colIndx[count]  = node_friends(i, j);

//             cos_a = node_vel_trans(i, j+1, 0);
//             sin_a = node_vel_trans(i, j+1, 1);

//             nzCoeffs_grad_x[count] = grad_coeffs(i, j+1, 0);
//             nzCoeffs_grad_y[count] = grad_coeffs(i, j+1, 1);

//             f = node_friends(i,j);

//             // Because the gradient operator returns a vector, we have to transform
//             // the vector into mapped coordinates. We do this here by premultiplying
//             // the components of the divergence operator by the transpose of the
//             // the transform matrix defined in node_vel_trans.
//             nzCoeffs_div_x[count]  =  cos_a*div_coeffs(i, j+1, 0) - sin_a*div_coeffs(i,  j+1, 1);
//             nzCoeffs_div_y[count]  =  sin_a*div_coeffs(i, j+1, 0) + cos_a*div_coeffs(i,  j+1, 1);

//             // nzCoeffs_div_x[count]  =  hVar[f] * (cos_a*div_coeffs(i, j+1, 0) - sin_a*div_coeffs(i,  j+1, 1));
//             // nzCoeffs_div_y[count]  =  hVar[f] * (sin_a*div_coeffs(i, j+1, 0) + cos_a*div_coeffs(i,  j+1, 1));

//             count++;
//         }


//     }

//     // assign the row indexes for the sparse matrix in CSR format

//     rowIndx[0] = 0;         // first element is always zero

//     // number of non-zero elements for pentagons
//     for (i=0; i<12; i++) rowIndx[i+1] = 6*(i+1);

//     // number of non-zero elements for hexagons
//     for (i=12; i<node_num; i++) rowIndx[i+1] = (12*6) + 7*(i-11);

//     // last element is always the number of non-zero coefficients
//     rowIndx[node_num] = nNonZero;

//     for (i=0; i<node_num; i++)
//     {
//       rowStart[i] = rowIndx[i];
//       rowEnd[i] = rowIndx[i+1];
//     }

//     //----------------- Set up sparse matrix for div and grad ------------------


//     // define handles to the sparse matrices for each component of the grad and
//     // div operators. The Laplacian matrix will be constructed from these.
//     sparse_matrix_t * DIV_X, * DIV_Y, * GRAD_X, * GRAD_Y;
//     sparse_matrix_t * LAP_X, * LAP_Y;

//     DIV_X   = new sparse_matrix_t;
//     DIV_Y   = new sparse_matrix_t;
//     GRAD_X  = new sparse_matrix_t;
//     GRAD_Y  = new sparse_matrix_t;
//     LAP_X   = new sparse_matrix_t;
//     LAP_Y   = new sparse_matrix_t;

//     // Assing memory to the handle for the Laplacian operator matrix
//     operatorLaplacian = new sparse_matrix_t;  // defined in mesh.h

//     // Define c (0) based indexing
//     sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;

//     // Define all operations the original (non-transpose) matrix
//     sparse_operation_t operation = SPARSE_OPERATION_NON_TRANSPOSE;

//     // Assign non-zero values to matrices for the components of div and grad operators
//     error = mkl_sparse_d_create_csr(DIV_X,  indexing, node_num, node_num, rowStart, rowEnd, colIndx, nzCoeffs_div_x);
//     error = mkl_sparse_d_create_csr(DIV_Y,  indexing, node_num, node_num, rowStart, rowEnd, colIndx, nzCoeffs_div_y);

//     error = mkl_sparse_d_create_csr(GRAD_X, indexing, node_num, node_num, rowStart, rowEnd, colIndx, nzCoeffs_grad_x);
//     error = mkl_sparse_d_create_csr(GRAD_Y, indexing, node_num, node_num, rowStart, rowEnd, colIndx, nzCoeffs_grad_y);

//     // Compute the multiplication of each component of the div dot grad
//     error = mkl_sparse_spmm(operation, *DIV_X, *GRAD_X, LAP_X); // div_x * grad_x
//     error = mkl_sparse_spmm(operation, *DIV_Y, *GRAD_Y, LAP_Y); // div_y * grad_y

//     // Sum each component to give the total Grad^2 operator (Laplacian)
//     // This is equivalent to (div_x * grad_x) + (div_y + grad_y)
//     error = mkl_sparse_d_add(operation, *LAP_X, 1.0, *LAP_Y, operatorLaplacian);

//     // matrix_descr descrp;
//     // descrp.type = SPARSE_MATRIX_TYPE_GENERAL;
//     // error = mkl_sparse_set_mv_hint (*operatorLaplacian, operation, descrp, 100000000);
//     // error = mkl_sparse_set_memory_hint (*operatorLaplacian, SPARSE_MEMORY_AGGRESSIVE);
//     // error = mkl_sparse_optimize (*operatorLaplacian);

//     // error = mkl_sparse_d_create_csr(operatorLaplacian, indexing, node_num, node_num, rowStart, rowEnd, colIndx, nzCoeffs);

//     // error = mkl_sparse_d_create_csr(operatorLaplacian,  indexing, node_num, node_num, rowStart, rowEnd, colIndx, nzCoeffs_div_x);
//     // error = mkl_sparse_d_create_csr(operatorLaplacian2,  indexing, node_num, node_num, rowStart, rowEnd, colIndx, nzCoeffs_div_y);


//     // Remove uneeded matrices from memory
//     delete DIV_X;
//     delete DIV_Y;
//     delete LAP_Y;
//     delete LAP_X;
//     delete GRAD_X;
//     delete GRAD_Y;


//     // Create objects for the row and column indices, and non-zero coefficients
//     // for the Laplace operator. Note that, although the Laplace operator was
//     // just create and assigned above, we do not yet know its structure or the
//     // number and values of the coefficients in its matrix.

//     int    * rowStartLap;
//     int      * rowEndLap;
//     int     * colIndxLap;
//     double * nzCoeffsLap;

//     int nc, nr; //nrows, ncols
//     int nNonZeroLap;

//     rowStartLap = NULL;
//     rowEndLap   = NULL;
//     colIndxLap  = NULL;
//     nzCoeffsLap = NULL;

//     // Now extract the the row and column index arrays, and nonzero coefficients
//     // from the Laplacian operator. These will be used to set up the solver for
//     // the problem L*p = d below, where L is the Laplacian operator and is to be
//     // inverted to find the solution p.
//     error = mkl_sparse_d_export_csr (*operatorLaplacian, &indexing, &nc, &nr, &rowStartLap, &rowEndLap, &colIndxLap, &nzCoeffsLap);

//     nNonZeroLap = rowEndLap[nr - 1];

//     for (i=0; i<node_num; i++)
//     {
//       rowIndx[i] = rowStartLap[i];
//     }
//     rowIndx[node_num] = nNonZeroLap;



//     // std::cout<<rowStartLap[node_num]<<std::endl;

//     //----------------- Set up the intel MKL solver for L*p = d ----------------

//     // create the solver object and assign it to a handle
//     MKL_INT opt = MKL_DSS_ZERO_BASED_INDEXING;
//     error = dss_create(pressureSolverHandle, opt);

//     if (error != MKL_DSS_SUCCESS)
//     {
//         outstring << "ERROR: MKL couldn't create the solver object and exited with error code "<<error<< std::endl;
//         globals->Output->Write(ERR_MESSAGE, &outstring);
//         globals->Output->TerminateODIS();
//     }

//     // define the structure of the sparse matrix of size N*N with nNonZero elements
//     opt =  MKL_DSS_SYMMETRIC_STRUCTURE;//MKL_DSS_SYMMETRIC_STRUCTURE;
//     // error = dss_define_structure(pressureSolverHandle, opt, rowIndx, node_num, node_num, colIndx, nNonZero);
//     error = dss_define_structure(pressureSolverHandle, opt, rowIndx, node_num, node_num, colIndxLap, nNonZeroLap);

//     if (error != MKL_DSS_SUCCESS)
//     {
//         outstring << "ERROR: MKL couldn't create the laplacian matrix structure and exited with error code "<<error<< std::endl;
//         globals->Output->Write(ERR_MESSAGE, &outstring);
//         globals->Output->TerminateODIS();
//     }

//     // reorder the matrix to optimize the permutation order
//     opt = MKL_DSS_AUTO_ORDER;
//     error = dss_reorder(pressureSolverHandle, opt, 0);

//     if (error != MKL_DSS_SUCCESS)
//     {
//         outstring << "ERROR: MKL couldn't reorder the laplacian matrix and exited with error code "<<error<< std::endl;
//         globals->Output->Write(ERR_MESSAGE, &outstring);
//         globals->Output->TerminateODIS();
//     }

//     // LU factorize the matrix using the nzCoeffs for direct solution
//     opt = MKL_DSS_POSITIVE_DEFINITE;
//     // error = dss_factor_real(pressureSolverHandle, opt, nzCoeffs);
//     error = dss_factor_real(pressureSolverHandle, opt, nzCoeffsLap);

//     if (error != MKL_DSS_SUCCESS)
//     {
//         outstring << "ERROR: MKL couldn't factorize the laplacian matrix and exited with error code "<<error<< std::endl;
//         globals->Output->Write(ERR_MESSAGE, &outstring);
//         globals->Output->TerminateODIS();
//     }

//     delete pressure_matrix;

//     delete nzCoeffs_grad_x;
//     delete nzCoeffs_grad_y;
//     delete nzCoeffs_div_x;
//     delete nzCoeffs_div_y;

//     return 1;
// }


// Function to read in text file containing mesh information
// The four pieces of information read and stored are:
//      1. Node ID numbers
//      2. Node positions in spherical coords
//      3. Neighbouring node ID numbers
//      4. All centroid positions in spherical coords
int Mesh::ReadMeshFile(void)
{
    std::string line, val;                 // strings for column and individual number
    std::string file_str;                  // string with path to mesh file.
    int i, node_id;

    // file_str = globals->path + SEP + "input_files" + SEP + "grid_l" + std::to_string(globals->geodesic_l.Value()) + ".txt";
    file_str = globals->path + SEP + "input_files" + SEP + "grid_l" + std::to_string(globals->geodesic_l.Value()) + ".txt";

    // in stream for input.in file
    std::ifstream gridFile(file_str, std::ifstream::in);

    if (gridFile.is_open())
    {
        outstring << std::endl << "Found mesh file: " + file_str << std::endl;
        globals->Output->Write(OUT_MESSAGE, &outstring);

        std::getline(gridFile, line);                               // READ HEADER
        while (std::getline(gridFile, line))
        {
            // std::cout<<line<<std::endl;
            std::istringstream line_ss(line);
            std::getline(line_ss >> std::ws,val,' ');                 // COL 0: Node ID number

            node_id = std::stoi(val);

            std::getline(line_ss >> std::ws,val,' ');                 // COL 1: Node Latitude
            node_pos_sph(node_id,0) = std::atof(val.c_str())*radConv;

            std::getline(line_ss >> std::ws,val,' ');                 // COL 2: Node Longitude
            node_pos_sph(node_id,1) = std::atof(val.c_str())*radConv;

            std::getline(line_ss >> std::ws,val,'{');                 // Read up to friends open bracket
            for (i = 0; i<5; i++)
            {
                std::getline(line_ss >> std::ws,val,',');               // Read each friend ID
                node_friends(node_id,i) = std::stoi(val);
            }
            std::getline(line_ss >> std::ws,val,'}');                 // Read last friend ID
            node_friends(node_id,5) = std::stoi(val);

            std::getline(line_ss >> std::ws,val,',');
            std::getline(line_ss >> std::ws,val,'{');                 // Read up to centroid list
            for (i = 0; i<5; i++)
            {
                std::getline(line_ss >> std::ws,val,'(');               // Read up to coord open bracket
                std::getline(line_ss >> std::ws,val,',');               // Read coord lat
                centroid_pos_sph(node_id,i,0) = std::atof(val.c_str())*radConv;

                std::getline(line_ss >> std::ws,val,')');               // Read coord lon
                centroid_pos_sph(node_id,i,1) = std::atof(val.c_str())*radConv;
                std::getline(line_ss >> std::ws,val,',');               // Read end of first coord
            }
            std::getline(line_ss >> std::ws,val,'(');                 // Read up to last coord open bracket
            std::getline(line_ss >> std::ws,val,',');                 // Read last coord lat
            centroid_pos_sph(node_id,5,0) = std::atof(val.c_str())*radConv;

            std::getline(line_ss >> std::ws,val,')');                 // Read last coord lon
            centroid_pos_sph(node_id,5,1) = std::atof(val.c_str())*radConv;

            std::getline(line_ss >> std::ws,val,'}');                 // Finish line


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
//     vInv_1D = new double[node_num * 6 * 6];
//     rot_1D = new double[node_num];
//
//     // Read in the data
//     dset_cellID.read( cellID_1D, PredType::NATIVE_INT, mspace_cellID, fspace_cellID );
//     dset_vInv.read( vInv_1D, PredType::NATIVE_DOUBLE, mspace_vInv, fspace_vInv );
//
//     // Load Array classes with 1D dynamic arrays
//     int count = 0;
//     for (i=0; i<node_num; i++)
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
//     double test_solution_gg[node_num];
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

int Mesh::ReadWeightingFile(void)
{
    int i, N_ll;


    N_ll = (int)globals->dLat.Value();

    const H5std_string DSET_cols("column index");
    const H5std_string DSET_rows("row index");
    const H5std_string DSET_data("weights");

    std::string file_str;
    file_str = globals->path + SEP
               + "input_files" + SEP
               + "grid_l" + std::to_string(globals->geodesic_l.Value())
               + '_' + std::to_string(N_ll)
               + 'x' + std::to_string(N_ll)
               + "_weights.h5";

    // Define file name and open file
    H5std_string FILE_NAME(file_str);
    H5File file(FILE_NAME, H5F_ACC_RDONLY);

    // Access dataspaces in latlon file
    DataSet dset_cols = file.openDataSet(DSET_cols);
    DataSet dset_rows = file.openDataSet(DSET_rows);
    DataSet dset_data = file.openDataSet(DSET_data);

    // Create filespaces for the correct rank and dimensions
    DataSpace fspace_cols = dset_cols.getSpace();
    DataSpace fspace_rows = dset_rows.getSpace();
    DataSpace fspace_data = dset_data.getSpace();

    // Get number of dimensions in the files dataspace
    int rank_cols = fspace_cols.getSimpleExtentNdims();
    int rank_rows = fspace_rows.getSimpleExtentNdims();
    int rank_data = fspace_data.getSimpleExtentNdims();

    hsize_t dims_cols[1];    // length no of non-zero elements
    hsize_t dims_rows[1];
    hsize_t dims_data[1];    // length no of non-zero elements

    // Get size of each dimension
    rank_cols = fspace_cols.getSimpleExtentDims( dims_cols );
    rank_rows = fspace_rows.getSimpleExtentDims( dims_rows );
    rank_data = fspace_data.getSimpleExtentDims( dims_data );

    // Create memoryspace to read the datasets
    DataSpace mspace_cols(1, dims_cols);
    DataSpace mspace_rows(1, dims_rows);
    DataSpace mspace_data(1, dims_data);

    // Create 1D arrays to store file data
    interpCols =    new    int[ dims_cols[0] ];
    interpRows =    new    int[ dims_rows[0] ];
    interpWeights = new double[ dims_data[0] ];

    // Read in the data
    dset_cols.read( interpCols, PredType::NATIVE_INT, mspace_cols, fspace_cols );
    dset_rows.read( interpRows, PredType::NATIVE_INT, mspace_rows, fspace_rows );
    dset_data.read( interpWeights, PredType::NATIVE_DOUBLE, mspace_data, fspace_data );


    // Create  matrix handle using loaded data
    // sparse_index_base_t index_type = SPARSE_INDEX_BASE_ZERO;     // we employ 0-based indexing.
    // sparse_status_t err;

    // sm1(rows,cols,nnz,outerIndexPtr, // read-write
    //                            innerIndices,values);
    

    int nrows = (360/N_ll)*(180/N_ll);
    int ncols = 3*node_num;

    interpMatrix = Eigen::Map<SpMat>(nrows, ncols, dims_data[0], interpRows, interpCols, interpWeights );


    // std::cout<<interpMatrix.cols()<<' '<<dims_data[0]<<std::endl;

    // interpMatrix.makeCompressed();

    // interpMatrix = SpMat(nrows, ncols);



    // // Allocate memory for non-zero entries
    // interpMatrix.reserve(dims_data[0]);

    // std::vector<Triplet> coefficients;
    // coefficients.reserve(dims_data[0]);


    // interpMatrix = new sparse_matrix_t;
    // err = mkl_sparse_d_create_csr(interpMatrix, index_type, nrows, ncols, interpRows, interpRows+1, interpCols, interpWeights);


    // // Sparse interpolation matrix successfully created. Now we must provide
    // // additional information to the matrix handle for optimization purposes

    // // Expected number of calls...?
    // // sparse_operation_t operation = SPARSE_OPERATION_NON_TRANSPOSE;
    // matrix_descr descrp;
    // descrp.type = SPARSE_MATRIX_TYPE_GENERAL;
    // descrp.mode = SPARSE_FILL_MODE_LOWER;
    // descrp.diag = SPARSE_DIAG_NON_UNIT;

    // sparse_operation_t operation = SPARSE_OPERATION_NON_TRANSPOSE;

    // int numberOfExpectedCalls = 100000000;  // Expect a big number!
    // err = mkl_sparse_set_dotmv_hint(*interpMatrix, operation, descrp, numberOfExpectedCalls);
    // err = mkl_sparse_set_memory_hint (*interpMatrix, SPARSE_MEMORY_AGGRESSIVE);
    // if (err>0) {
    //     std::cout<<"Intel matrix optimization error!"<<std::endl;
    // }

    // // Now optimize matrix
    // err = mkl_sparse_optimize(*interpMatrix);

    return 1;

};