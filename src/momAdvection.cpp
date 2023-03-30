#include "mesh.h"
#include "globals.h"
#include "array1d.h"
#include "interpolation.h"
#include "gridConstants.h"


#include <Eigen/Dense>

void calculateMomentumAdvection(Globals * globals,
                                Mesh * mesh,
                                Array1D<double> & dvdt,
                                Array1D<double> & vel,
                                Array1D<double> & thickness_n,
                                Array1D<double> & Ekin)
{
    double rot_rate = globals->angVel.Value();

    // -------------------------------------------------------
    // Step 1: Calculate curl of the velocity field, given at 
    // vertex points
    // -------------------------------------------------------

    Array1D<double> vorticity_v(VERTEX_NUM);

    Eigen::Map<Eigen::VectorXd> vel_eig(&vel(0), FACE_NUM);
    Eigen::Map<Eigen::VectorXd> vorticity_eig(&vorticity_v(0), VERTEX_NUM);

    // Perform sparse matrix * vector operation
    vorticity_eig = mesh->operatorCurl * vel_eig;

    // -------------------------------------------------------
    // Step 2: Calculate absolute vorticity field, given at
    // vertex points, and
    // -------------------------------------------------------
    // Step 3: Calculate potential voritcity field using 
    // thickness field interpolated to vertexes
    // -------------------------------------------------------
    {
        double f;
        double thickness;
        double node_ID;
        for (int i=0; i<VERTEX_NUM; i++) {
            f = -2 * rot_rate * mesh->vertex_sinlat(i);

#ifdef TEST_SW2    
            // SWE TEST CASE 2. In this case the coriolis parameter has to be
            // redifined with Eq 96 in Williamson et al (1992). 
            double lat = mesh->vertex_pos_sph(i, 0);
            double lon = mesh->vertex_pos_sph(i, 1);
            double a = pi/4.0;
            f = -2 * rot_rate * (-cos(lon)*cos(lat)*sin(a) + sin(lat)*cos(a));
#endif
 
            thickness = 0.0;
            for (int j=0; j<3; j++)
            {
                node_ID = mesh->vertex_nodes(i, j);

                // Eq. 25, Ringler et al (2010)
                thickness += thickness_n(node_ID) * mesh->cv_area_sph(node_ID) * mesh->vertex_R(i, j);
            }

            thickness *= mesh->vertex_area_r(i);

            vorticity_v(i) = (vorticity_v(i) + f)/thickness;
        }
    }

    // -------------------------------------------------------
    // Step 4: Calculate vorticity from vertii to edges
    // -------------------------------------------------------
    // -------------------------------------------------------
    // Step 5: Calculate thickness from node to edges
    // -------------------------------------------------------


    Array1D<double> vorticity_e(FACE_NUM); //vorticity at faces
    Array1D<double> thickness_e(FACE_NUM);
    {
        int v1, v2;
        int n1, n2;
        for (int i=0; i<FACE_NUM; i++) {
            v1 = mesh->face_vertexes(i, 0);
            v2 = mesh->face_vertexes(i, 1);

            vorticity_e(i) = 0.5 * (vorticity_v(v1) + vorticity_v(v2) );
        
            n1 = mesh->face_nodes(i, 0);
            n2 = mesh->face_nodes(i, 1);

            thickness_e(i) = 0.5 * (thickness_n(n1) + thickness_n(n2));
        }
    }

    // -------------------------------------------------------
    // Step 6: Calculate vorticity flux - note this is very
    // similiar to the coriolis operator
    // -------------------------------------------------------
    {   
        int n1, n2, f_ID;
        int fnum[2];
        double q_e, F_tang_q, q_e2, F_e, coeff;
        for (int i = 0; i < FACE_NUM; ++i)
        {
            // Loops similar to coriolis and interpolation operators
            n1 = mesh->face_nodes(i, 0);
            n2 = mesh->face_nodes(i, 1);

            fnum[0] = mesh->node_fnum(n1);
            fnum[1] = mesh->node_fnum(n2);

            q_e = vorticity_e(i);

            // vorticity flux at face i
            F_tang_q = 0.0;
            for (unsigned k=0; k<2; k++) {
                for (unsigned j = 0; j < fnum[k]-1; j++)
                {
                    // friend face ID
                    f_ID = mesh->face_friends(i, k, j);

                    // vorticity at edge of friend
                    q_e2 = vorticity_e(f_ID);

                    // thickness flux at edge of friend
                    F_e = thickness_e(f_ID) * vel(f_ID);

                    // Eq. 49 from Ringler et al (2010)
                    //                      w_ee                             l_e                      
                    coeff = mesh->face_interp_weights(i, k, j) * mesh->face_len(f_ID);
                    

                    F_tang_q += coeff * F_e * (q_e + q_e2);

                }
            }
            
            //                                 1/d_e
            dvdt(i) -= -F_tang_q * 0.5 * mesh->face_node_dist_r(i);

        }
    }


    // -------------------------------------------------------
    // Step 7: Calculate kinetic energy of each control 
    // volume
    // -------------------------------------------------------

    

    Array2D<double> vel_xyz(NODE_NUM, 3);
    interpolateVelocityCartRBF(globals, mesh, vel_xyz, vel);

    for (int i=0; i<NODE_NUM; i++) {
        Ekin(i) = 0.5 * (vel_xyz(i,0)*vel_xyz(i,0) + vel_xyz(i,1)*vel_xyz(i,1)+ vel_xyz(i,2)*vel_xyz(i,2));
    }

    // -------------------------------------------------------
    // Step 8: Remove the gradient of the kinetic energy 
    // from the advection term to account for the non-rotational
    // part of momentum advection. (as it appears on the RHS)
    // -------------------------------------------------------

    Eigen::Map<Eigen::VectorXd> dvdt_eig(&dvdt(0), FACE_NUM);
    Eigen::Map<Eigen::VectorXd> Ek_eig(&Ekin(0), NODE_NUM);

    // Perform sparse matrix * vector operation
    dvdt_eig -= mesh->operatorGradient * Ek_eig;

    // Old GRAD(KE) calculation
    // Array1D<double> Ekin(NODE_NUM);
    // Array1D<double> vt(face_num);

    // calculateKE(globals, mesh, Ekin, vel);

    // for (int i=0; i<FACE_NUM; ++i)
    // {
    // //     double v_tang = 0.0;
    // //     double v_norm = 0.0;
    // //     int friend_num = 10;

    // //     int n1, n2;
    // //     n1 = mesh->face_nodes(i,0);
    // //     n2 = mesh->face_nodes(i,1);
    // //     if (mesh->node_friends(n1,5) < 0) friend_num--;
    // //     if (mesh->node_friends(n2,5) < 0) friend_num--;

    // //     for (int j=0; j<friend_num; j++) {
    // //         int f_ID = mesh->face_interp_friends(i, j);
    // //         v_tang += vel(f_ID)*mesh->face_interp_weights(i,j)*mesh->face_len(f_ID);
    // //     }
    // //     v_tang /= mesh->face_node_dist(i);

    // //     vt(i) = v_tang;

    //     // advection(i) = 0.0;
    // }

    // for (int i=0; i<NODE_NUM; i++) {
    //     int friend_num = 6;
    //     if (mesh->node_friends(i, 5) < 0)
    //         friend_num--;
        
    //     double E = 0.0;
    //     double Ai = 0.0;
    //     double Ae;
    //     for (int j = 0; j < friend_num; j++)
    //     {
    //         int f_ID = mesh->faces(i, j);

    //         double ue = vel(f_ID);
    //         // double v, u, ut,vt, nx, ny;
    //         // double cosa, sina;
    //         Ae = mesh->face_area(f_ID);
    //         // nx = mesh->control_vol_edge_normal_map(i, j, 0);
    //         // ny = mesh->control_vol_edge_normal_map(i, j, 1);
    //         // cosa = 1.0;//mesh->node_face_vel_trans(i,j,0);
    //         // sina = 0.0;//mesh->node_face_vel_trans(i,j,1);

    //         // ut = ue*nx; 
    //         // vt = ue*ny;

    //         // u = cosa*ut + sina*vt;
    //         // v = -sina*ut + cosa*vt; 

    //         // std::cout<

    //         // ny = mesh->face_node_*ue;
    //         // Ae = mesh->node_friend_element_areas_map(i,j,0)+mesh->node_friend_element_areas_map(i,j,1)+mesh->node_friend_element_areas_map(i,j,2);

            
            
    //         // if (i<15) 
    //         // {
    //         //     std::cout<<i<<' '<<j<<' '<<ue<<' '<<Ae<<std::endl;
    //         // }

    //         // if (i==5) std::cout<<i<<' '<<mesh->node_friends(i,j)<<' '<<ue<<std::endl;

    //         E += Ae * ue*ue*0.25;
    //         // // E += Ae * (pow(u,2.0) + pow(v,2.0))*0.25;
    //         // // E += ()
    //         // Ai += Ae4;
    //     }

    //     E /= mesh->control_volume_surf_area_map(i);

    // //     // E = 0.5*vel()

    // //     // if (i < 12 || (i > 10267 && i < 10272)) std::cout<<i<<' '<<E<<std::endl;
        

        

    //     Ekin(i) = E;

    // }
};