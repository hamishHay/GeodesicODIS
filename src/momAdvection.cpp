#include "mesh.h"
#include "globals.h"
#include "array1d.h"

#include <Eigen/Dense>

void calculateMomentumAdvection(Globals * globals,
                                Mesh * mesh,
                                Array1D<double> & vel,
                                Array1D<double> & thickness_n,
                                Array1D<double> & advection)
{
    int node_num = globals->node_num;
    int face_num = globals->face_num;
    int vertex_num = globals->vertex_num;

    double rot_rate = globals->angVel.Value();

    // -------------------------------------------------------
    // Step 1: Calculate curl of the velocity field, given at 
    // vertex points
    // -------------------------------------------------------

    Array1D<double> vorticity_v(vertex_num);

    Eigen::Map<Eigen::VectorXd> vel_eig(&vel(0), face_num);
    Eigen::Map<Eigen::VectorXd> vorticity_eig(&vorticity_v(0), mesh->vertex_num);

    // Perform sparse matrix * vector operation
    vorticity_eig = mesh->operatorCurl * vel_eig;


    // -------------------------------------------------------
    // Step 2: Calculate absolute vorticity field, given at
    // vertex points
    // -------------------------------------------------------

    for (int i=0; i<vertex_num; i++) {
        double f = -2 * rot_rate * sin(mesh->vertex_pos_sph(i, 0));
        
        // SWE TEST CASE 2:
        // double lat = mesh->vertex_pos_sph(i, 0);
        // double lon = mesh->vertex_pos_sph(i, 1);
        // double a = 0.3;
        // f = 2 * rot_rate * (-cos(lon)*cos(lat)*sin(a) + sin(lat)*cos(a));

        vorticity_v(i) += f;
    }

    // -------------------------------------------------------
    // Step 3: Calculate potential voritcity field using 
    // thickness field interpolated to vertexes
    // -------------------------------------------------------

    Array1D<double> vertex_thickness(vertex_num);

    for (int i=0; i<vertex_num; i++) {
        double thickness = 0.0;

        for (int j=0; j<3; j++)
        {
            int node_ID = mesh->vertex_nodes(i, j);
            double h = thickness_n(node_ID);

            // Eq. 25, Ringler et al (2010)
            thickness += h * mesh->control_volume_surf_area_map(node_ID) * mesh->vertex_R(i, j);
        }

        thickness /= mesh->vertex_area(i);


        vorticity_v(i) /= thickness;
    }

    // -------------------------------------------------------
    // Step 4: Calculate vorticity from vertii to edges
    // -------------------------------------------------------

    Array1D<double> vorticity_e(face_num); //vorticity at faces

    for (int i=0; i<face_num; i++) {
        int v1 = mesh->face_vertexes(i, 0);
        int v2 = mesh->face_vertexes(i, 1);

        vorticity_e(i) = 0.5 * (vorticity_v(v1) + vorticity_v(v2) );
    }

    // -------------------------------------------------------
    // Step 5: Calculate thickness from node t- edges
    // -------------------------------------------------------

    Array1D<double> thickness_e(face_num); //vorticity at faces

    for (int i = 0; i < face_num; i++)
    {
        int n1 = mesh->face_nodes(i, 0);
        int n2 = mesh->face_nodes(i, 1);

        thickness_e(i) = 0.5 * (thickness_n(n1) + thickness_n(n2));
    }

    // -------------------------------------------------------
    // Step 6: Calculate vorticity flux - note this is very
    // similiar to the coriolis operator
    // -------------------------------------------------------

    for (int i = 0; i < face_num; ++i)
    {
        // This inner loop goes around the faces adjoining the 
        // two control volumes of face i - we must find if either
        // of these control volumes are pentagons here.
        int friend_num = 10;

        int n1, n2;
        n1 = mesh->face_nodes(i, 0);
        n2 = mesh->face_nodes(i, 1);
        
        if (mesh->node_friends(n1, 5) < 0)
            friend_num--;
        if (mesh->node_friends(n2, 5) < 0)
            friend_num--;

        double q_e = vorticity_e(i);

        // vorticity flux at face i
        double F_tang_q = 0.0;
        for (int j = 0; j < friend_num; j++)
        {
            // friend face ID
            int f_ID = mesh->face_interp_friends(i, j);

            // vorticity at edge of friend
            double q_e2 = vorticity_e(f_ID);

            // thickness flux at edge of friend
            double F_e = thickness_e(f_ID) * vel(f_ID);

            // Eq. 49 from Ringler et al (2010)
            //                      w_ee                             l_e                      d_e 
            double coeff = mesh->face_interp_weights(i, j) * mesh->face_len(f_ID) / mesh->face_node_dist(i);
            

            F_tang_q += coeff * F_e * (q_e + q_e2)*0.5;



            // F_tang_q += -coeff * vel(f_ID) * 2 * rot_rate * sin(mesh->face_intercept_pos_sph(i, 0));

            // std::cout << i << ' ' << vel(f_ID) << std::endl;
        }

        advection(i) = F_tang_q;

        

        // advection(i) += coeff * vel(f_ID) 
    }

    // globals->Output->TerminateODIS();

    // -------------------------------------------------------
    // Step 7: Calculate kinetic energy of each control 
    // volume
    // -------------------------------------------------------

    Array1D<double> Ekin(node_num);

    for (int i=0; i<node_num; i++) {
        int friend_num = 6;
        if (mesh->node_friends(i, 5) < 0)
            friend_num--;
        
        double E = 0.0;
        for (int j = 0; j < friend_num; j++)
        {
            int f_ID = mesh->faces(i, j);

            double ue = vel(f_ID);
            double Ae = mesh->face_area(f_ID);

            E += Ae * ue * ue;
        }

        E /= 4.0 * mesh->control_volume_surf_area_map(i);

        Ekin(i) = E;

    }

    // -------------------------------------------------------
    // Step 8: Remove the gradient of the kinetic energy 
    // from the advection term to account for the non-rotational
    // part of momentum advection. (as it appears on the RHS)
    // -------------------------------------------------------

    Eigen::Map<Eigen::VectorXd> advection_eig(&advection(0), face_num);
    Eigen::Map<Eigen::VectorXd> Ek_eig(&Ekin(0), node_num);

    // Perform sparse matrix * vector operation
    advection_eig -= mesh->operatorGradient * Ek_eig;
};