#include "mesh.h"
#include "globals.h"
#include "array3d.h"
#include "array2d.h"
#include "array1d.h"
#include "spatialOperators.h"
#include <math.h>
#include <iostream>
#include <mkl.h>

void setBeta(Array1D<double> & beta, int m, int n, int node_num, Array2D<double> & sph)
{
    int i;
    double lon, lat;
    for (i = 0 ; i < node_num; i++)
    {
        lat = sph(i, 0);
        lon = sph(i, 1);
        // beta(i) = 1e6*trigMLon(i, m, 0) * trigNLat(i, n, 0) * trigNLat(i, n, 0)
        //           * trigNLat(i, n, 0) * trigNLat(i, n, 0);
      beta(i) = 1e6*cos(m*lon) * pow(cos(n*lat), 4.0);
    }
};

void setU(Array2D<double> & u, int m, int n, int node_num, Array2D<double> & sph)
{
    int i;
    double lon, lat;

    for (i = 0 ; i < node_num; i++)
    {
        lat = sph(i, 0);
        lon = sph(i, 1);

        u(i, 0) = -4. * 1e6* n * sin(lon) * cos(m*lon);
        u(i, 0) *= sin(n*lat);
        u(i, 0) *= pow(cos(n*lat), 3.0);


        u(i, 1) = -m * 1e6*sin(lon) * sin(m*lon);
        u(i, 1) *= pow(cos(n*lat),4.0)/cos(lat);

    }
};

void setDivU(Array1D<double> & divU, int m, int n, int node_num, double r, Array2D<double> & sph)
{
    int i;
    double lon, lat;

    for (i = 0 ; i < node_num; i++)
    {
        lat = sph(i, 0);
        lon = sph(i, 1);

        divU(i)  = 2.*m*sin(lon)*sin(m*lon);
        divU(i) -= cos(lon)*cos(m*lon);
        divU(i) *= sin(n*lat);
        divU(i) *= pow(cos(n*lat), 3.0);
        divU(i) *= 1e6*4.*n/(r * cos(lat));
    }
};

void setGradBeta(Array2D<double> & gradBeta, int m, int n, int node_num, double r,  Array2D<double> & sph)
{
    int i;
    double lon, lat;
    for (i = 0 ; i < node_num; i++)
    {
        lat = sph(i, 0);
        lon = sph(i, 1);

        gradBeta(i, 1) =  -1e6 * 4.0 * n * cos(m*lon);
        gradBeta(i, 1) *= sin(n*lat);
        gradBeta(i, 1) *= pow(cos(n*lat), 3.0);
        gradBeta(i, 1) /= r;

        gradBeta(i, 0) = -m * 1e6 * sin(m*lon);
        gradBeta(i, 0) *= pow(cos(n*lat), 4.0);
        gradBeta(i, 0) /= r*cos(lat);

        // gradBeta(i, 1) =  -1e6 * 4.0 * n * trigMLon(i, m, 0);
        // gradBeta(i, 1) *= trigNLat(i, n, 1);
        // gradBeta(i, 1) *= trigNLat(i, n, 0)*trigNLat(i, n, 0)*trigNLat(i, n, 0);
        // gradBeta(i, 1) /= r;
        //
        // gradBeta(i, 0) = -m * 1e6 * trigMLon(i, m, 1);
        // gradBeta(i, 0) *= trigNLat(i, n, 0)*trigNLat(i, n, 0)
        //            *trigNLat(i, n, 0)*trigNLat(i, n, 0);
        // gradBeta(i, 0) /= r*trigLat(i, 0);
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
    // Array2D<double> * u;

    Array1D<double> gradTest2(mesh->face_num);
    Array1D<double> beta_faces(mesh->face_num);
    Array2D<double> beta_grad(mesh->face_num, 2);
    Array1D<double> beta_nodes(mesh->node_num);
    Array1D<double> div_nodes(mesh->node_num);
    Array1D<double> div_nodes_test(mesh->node_num);
    Array2D<double> u_faces(mesh->face_num, 2);
    Array2D<double> u(mesh->node_num, 2);
    Array1D<double> u_faces_normal(mesh->face_num);

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
    // u           = new Array2D<double>(node_num, 2);

    gradBeta    = new Array2D<double>(node_num, 2);
    gradTest    = new Array2D<double>(node_num, 2);
    laplaceBeta = new Array1D<double>(node_num);
    divU        = new Array1D<double>(node_num);
    divTest     = new Array1D<double>(node_num);


    // trigLon     = &(mesh->trigLon);
    // trigLat     = &(mesh->trigLat);
    // trigSqLon   = &(mesh->trigSqLon);
    // trigSqLat   = &(mesh->trigSqLat);
    // trigMLon    = &(mesh->trigMLon);
    // trigNLat    = new Array3D<double>(node_num, nMax+1, 2);
    //
    // for (i=0; i < node_num; i++)
    // {
    //     for (n=0; n < nMax+1; n++)
    //     {
    //         (*trigNLat)(i,n,0) = cos(mesh->node_pos_sph(i,0));
    //         (*trigNLat)(i,n,1) = sin(mesh->node_pos_sph(i,0));
    //     }
    // }

    int a = 0;
    double sum = -1.0;
    for (n=nMin; n < nMax + 1; n++)
    {
        for (m = mMin; m < mMax + 1; m++)
        {
            setBeta(beta_faces, m, n, mesh->face_num, mesh->face_centre_pos_sph);
            setBeta(beta_nodes, m, n, mesh->node_num, mesh->node_pos_sph);

            setU(u_faces, m, n, mesh->face_num, mesh->face_centre_pos_sph);

            setU(u, m, n, mesh->node_num, mesh->node_pos_sph);

            // for (i=0; i<mesh->face_num; i++) {
            //     double u_temp = u_faces(i,0);
            //     double v_temp = u_faces(i,1);
            //     double u1, v1;
            //     double nx, ny;
            //     double cos_a, sin_a;
            //
            //     // CONVERT TO MAPPED VELOCITIES
            //     // cos_a = mesh->node_face_vel_trans(i, j, 0);
            //     // sin_a = mesh->node_face_vel_trans(i, j, 1);
            //     // u1 = u_temp * cos_a + v_temp * sin_a;
            //     // v1 = -u_temp * sin_a + v_temp * cos_a;
            //     u1 = u_temp;
            //     v1 = v_temp;
            //
            //     nx = mesh->face_normal_vec_map(i, 0);
            //     ny = mesh->face_normal_vec_map(i, 1);
            //
            //     u_faces_normal(i) = (nx*u1 + ny*v1);
            //     // u_faces(i,1) = (nx*u1 + ny*v1)*ny;
            //
            // }

            setDivU(div_nodes, m, n, node_num, r, mesh->node_pos_sph);


            setGradBeta(beta_grad, m, n, mesh->face_num, r, mesh->face_centre_pos_sph);

            // velocityDivergence(mesh, *divTest, u, sum, -1.0);

            for (i=0; i<mesh->node_num; i++) {
                int friend_num = 6;
                if (mesh->node_friends(i,5) < 0) friend_num--;

                double div = 0.0;
                double area = mesh->control_volume_surf_area_map(i);
                for (int j=0; j<friend_num; j++) {
                    int face_id = mesh->faces(i, j);

                    double u_temp = u_faces(face_id,0);
                    double v_temp = u_faces(face_id,1);
                    double u1, v1;
                    double nx, ny;
                    double cos_a, sin_a;

                    // CONVERT TO MAPPED VELOCITIES
                    // We do not have to convert to mapped velocities because
                    // the velocity vector is already normal to the control
                    // volume face!
                    // cos_a = mesh->node_face_vel_trans(i, j, 0);
                    // sin_a = mesh->node_face_vel_trans(i, j, 1);
                    // u1 = u_temp * cos_a + v_temp * sin_a;
                    // v1 = -u_temp * sin_a + v_temp * cos_a;

                    u1 = u_temp;
                    v1 = v_temp;
                    nx = mesh->face_normal_vec_map(face_id, 0);
                    ny = mesh->face_normal_vec_map(face_id, 1);

                    double len = mesh->face_len(face_id);

                    double vel = nx*u1 + ny*v1;

                    int dir = mesh->node_face_dir(i, j);
                    // if (i==753 && face_id ==4166) dir=1;
                    div += vel*len*dir;

                    if (i==973 || i==974) std::cout<<i<<' '<<face_id<<' '<<dir<<' '<<mesh->face_nodes(face_id,0)<<' '<<mesh->face_nodes(face_id,1)<<' '<<len<<' '<<nx<<' '<<ny<<std::endl;

                    // if (mesh->node_face_dir(i, j) < 0) {
                    //     std::cout<<"NODE "<<i<<std::endl;
                    //     break;
                    // }
                }
                div /= area;
                div_nodes_test(i) = div;
                // if (i==0) std::cout<<std::endl;

                // if (true) {
                //     std::cout<<i<<' ';
                //     std::cout<<div_nodes(i)<<'\t';
                //     // std::cout<<(*gradBeta)(i,1)<<'\t';
                //     // std::cout<<(*u)(i,a)/(r*(*trigLon)(i,1))<<'\t';
                //     std::cout<<div<<'\t';
                //     std::cout<<fabs(div_nodes(i) - div)/((div_nodes(i)))*100.0<<'\t';
                //     // std::cout<<(*gradTest)(i,1)<<'\t';
                //     // std::cout<<du_vel[2*i]<<'\t';
                //     // std::cout<<(*gradTest)(i,1)<<'\t';
                //     // std::cout<<fabs((*gradBeta)(i,a) - (*gradTest)(i,a))/(*gradBeta)(i,a)*100.0<<std::endl;
                //
                //     std::cout<<std::endl;
                //     // std::cout<<(*divU)(i)<<'\t';
                //     // std::cout<<(*divTest)(i)<<std::endl;
                //     // std::cout<<(*divU)(i) - (*divTest)(i)<<std::endl;
                //
                // }
            }

            for (i=0; i<mesh->face_num; i++) {
                double inner, outer;
                int node_in, node_out;
                int dir_in, dir_out;
                double dist;

                node_in = mesh->face_nodes(i, 0);
                node_out = mesh->face_nodes(i, 1);
                dist = mesh->face_node_dist(i);

                inner = beta_nodes(node_in);
                outer = beta_nodes(node_out);

                // if (i==30661) std::cout<<dist<<' '<<node_in<<' '<<node_out<<' '<<(inner-outer)/dist<<std::endl;

                gradTest2(i) = -(inner*mesh->face_centre_m(i, 0) - outer*mesh->face_centre_m(i, 1))/dist;

                double u_temp = beta_grad(i,0);
                double v_temp = beta_grad(i,1);
                double u1, v1;
                double nx, ny;
                double cos_a, sin_a;

                // CONVERT TO MAPPED VELOCITIES
                // cos_a = mesh->node_vel_trans(mesh->face_nodes(i, 0), mesh->face_friends(i,0), 0);
                // sin_a = mesh->node_vel_trans(mesh->face_nodes(i, 0), mesh->face_friends(i,0), 1);
                u1 = u_temp * cos_a + v_temp * sin_a;
                v1 = -u_temp * sin_a + v_temp * cos_a;

                nx = mesh->face_normal_vec_map(i, 0);
                ny = mesh->face_normal_vec_map(i, 1);

                if (true) {
                    std::cout<<i<<' ';
                    std::cout<<u1*nx + v1*ny<<'\t';
                    // std::cout<<(*gradBeta)(i,1)<<'\t';
                    // std::cout<<(*u)(i,a)/(r*(*trigLon)(i,1))<<'\t';
                    std::cout<<gradTest2(i)<<'\t';
                    std::cout<<fabs((u1*nx + v1*ny) - gradTest2(i))/((u1*nx + v1*ny))*100.0<<'\t';
                    // std::cout<<(*gradTest)(i,1)<<'\t';
                    // std::cout<<du_vel[2*i]<<'\t';
                    // std::cout<<(*gradTest)(i,1)<<'\t';
                    // std::cout<<fabs((*gradBeta)(i,a) - (*gradTest)(i,a))/(*gradBeta)(i,a)*100.0<<std::endl;

                    std::cout<<std::endl;
                    // std::cout<<(*divU)(i)<<'\t';
                    // std::cout<<(*divTest)(i)<<std::endl;
                    // std::cout<<(*divU)(i) - (*divTest)(i)<<std::endl;

                }
            }



            // pressureGradient(mesh, *gradTest, *beta, node_num, -1.0);

            // sparse_operation_t operation = SPARSE_OPERATION_NON_TRANSPOSE;
            // matrix_descr descript;
            // // struct matrix_descr {
            // //     sparse_matrix_type_t type;
            // // } descript;
            // descript.type = SPARSE_MATRIX_TYPE_GENERAL;
            // double alpham = 1.0;
            // double betam = 1.0;
            // int error;
            //
            // // double * u_vel, * v_vel;
            // // u_vel = new double[node_num];
            // // v_vel = new double[node_num];
            // double * du_vel, * dv_vel;
            // du_vel = new double[2*node_num];
            // // dv_vel = new double[node_num];
            // double * vec = new double[3*node_num];
            // for (i=0; i<node_num; i++)
            // {
            //     // dvdt(i,0) = 0.0;
            //     // dvdt(i,1) = 0.0;
            //     // pressure(i) = 0.001*i;
            //     // u_vel[i] = velocity(i,0);
            //     // v_vel[i] = velocity(i,1);
            //     du_vel[2*i] = 0.0;
            //     du_vel[2*i+1] = 0.0;
            //     // dv_vel[i] = dvdt(i,1);
            //     vec[3*i] = 0.0;
            //     vec[3*i+1] = 0.0;
            //     vec[3*i+2] = (*beta)(i);
            // }
            //
            // error = mkl_sparse_d_mv(operation, -1/1.308, *(mesh->operatorGradient), descript, vec, betam, du_vel);

            // for (i = 0; i < mesh->node_num; i++) std::cout<<i<<' '<<beta_nodes(i)<<std::endl;

            // for (i = 0; i < mesh->node_num; i++)
            // {
            //     double u_temp = beta_grad(i,0);
            //     double v_temp = beta_grad(i,1);
            //     double u1, v1;
            //     double nx, ny;
            //     double cos_a, sin_a;
            //
            //     // CONVERT TO MAPPED VELOCITIES
            //     cos_a = mesh->node_vel_trans(mesh->face_nodes(i, 0), mesh->face_friends(i,0), 0);
            //     sin_a = mesh->node_vel_trans(mesh->face_nodes(i, 0), mesh->face_friends(i,0), 1);
            //     u1 = u_temp * cos_a + v_temp * sin_a;
            //     v1 = -u_temp * sin_a + v_temp * cos_a;
            //
            //     nx = mesh->face_normal_vec_map(i, 0);
            //     ny = mesh->face_normal_vec_map(i, 1);
            //
            //     if (true) {//i==973 || i==974) {
            //         // std::cout<<i<<' '<<mesh->face_centre_pos_sph(i,0)<<' '<<mesh->face_centre_pos_sph(i,1)<<' ';
            //         std::cout<<i<<' '<<div_nodes(i)<<'\t';
            //         std::cout<<div_nodes_test(i)<<'\t';
            //         std::cout<<fabs(div_nodes(i) - div_nodes_test(i))/fabs(div_nodes(i))*100.0<<'\t';
            //         // std::cout<<(*gradBeta)(i,1)<<'\t';
            //         // std::cout<<(*u)(i,a)/(r*(*trigLon)(i,1))<<'\t';
            //         // std::cout<<gradTest2(i)<<'\t';
            //         // std::cout<<(*gradTest)(i,1)<<'\t';
            //         // std::cout<<du_vel[2*i]<<'\t';
            //         // std::cout<<(*gradTest)(i,1)<<'\t';
            //         // std::cout<<fabs((*gradBeta)(i,a) - (*gradTest)(i,a))/(*gradBeta)(i,a)*100.0<<std::endl;
            //
            //         std::cout<<std::endl;
            //         // std::cout<<(*divU)(i)<<'\t';
            //         // std::cout<<(*divTest)(i)<<std::endl;
            //         // std::cout<<(*divU)(i) - (*divTest)(i)<<std::endl;
            //
            //     }
            //
            //
            //     (*gradTest)(i,0) = 0.0;
            //     (*gradTest)(i,1) = 0.0;
            //     (*divTest)(i) = 0.0;
            // }

            // 10182, 10183

                // set u vector
                // set laplacian solution
                // set div solution
                // set grad solution
        }
    }


    // Define vector u

    // Define analytical solutions to grad beta, laplacian beta, and div u
};
