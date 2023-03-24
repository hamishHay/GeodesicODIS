#include "energy.h"
#include "array1d.h"
#include "array2d.h"
#include "globals.h"
#include "mesh.h"
#include "interpolation.h"
#include "gridConstants.h"
#include <math.h>
#include <iostream>
#include <sstream>


void updateEnergy(Globals * globals, double & avg_flux, Array1D<double> & e_flux, Array2D<double> & vel, Array1D<double> & areas)
{
    int node_num, face_num, i;
    double drag_coeff, h, r;
    // double avg_flux;

    drag_coeff = globals->alpha.Value();
    h = globals->h.Value();
    r = globals->radius.Value();

    node_num = globals->node_num;
    face_num = globals->face_num;

    avg_flux = 0.0;

    // UPDATE ENERGY
    switch (globals->fric_type)
    {
      case LINEAR:
          for (i = 0; i<FACE_NUM; i++)
          {
            e_flux(i) = drag_coeff * 1000.0 * h * (vel(i,0)*vel(i,0) + vel(i,1)*vel(i,1));

            avg_flux += e_flux(i)*areas(i);

            // e_flux(i) *= 1./mass(i) * 1000.0 * h; // Convert to flux?



            // flux += 1000.0 * h * drag_coeff * (v(i,0)*v(i,0) + v(i,1)*v(i,1));
        }
        break;

      case QUADRATIC:
          for (i = 0; i<FACE_NUM; i++)
          {
            e_flux(i) = drag_coeff/h * sqrt(vel(i,0)*vel(i,0) + vel(i,1)*vel(i,1))
                                * (vel(i,0)*vel(i,0) + vel(i,1)*vel(i,1));


            avg_flux += e_flux(i)*areas(i);
        }
        break;
    }

    // avg_flux /= face_num;
    // TODO - convert to area weighted average
    avg_flux /= (4* pi * pow(r, 2.0));

};


void calculateKE(Globals * globals, Mesh * mesh, Array1D<double> & ke, Array1D<double> & velocity)
{
  int vertex_num = mesh->vertex_num;
  int face_num = mesh->face_num;
  int node_num = mesh->node_num;
  int i, j;
  double u, v;

  // Array1D<double> KE(node_num);
  Array2D<double> vel_interp(NODE_NUM, 2);

  // interpolateVelocityRBF(globals, mesh, vel_interp, velocity);



  // for (i=0; i<node_num; i++) {
  //   // std::cout<<i<<' '<<pow(vel_interp(i,0),2.0)+pow(vel_interp(i, 1),2.0)<<std::endl;

  //   ke(i) = 0.5*(pow(vel_interp(i,0),2.0)+pow(vel_interp(i, 1),2.0));
  //   // ke(i) = vel_interp(i,1);

  //   // std::cout<<ke(i)<<std::endl;
  //   // ke(i) = vel_interp(i,0);
  // }
  // for (i=0; i<vertex_num; i++){
  //   double total_area=0.0;

  //   double u_vertex = 0.0;
  //   double v_vertex = 0.0;

  //   total_area = mesh->vertex_face_area(i, 0)+mesh->vertex_face_area(i, 1)+mesh->vertex_face_area(i, 2);

  //   // std::cout<<i<<' ';

  //   for (j=0; j<3; j++) {
  //     int f = mesh->vertex_faces(i, j);

  //     double nx, ny;

  //     nx = mesh->face_normal_vec_map(f, 0);
  //     ny = mesh->face_normal_vec_map(f, 1);

  //     u = velocity(f)*mesh->face_normal_vec_map(f, 0);//*mesh->vertex_face_dir(i,(j)%3);
  //     v = velocity(f)*mesh->face_normal_vec_map(f, 1);//*mesh->vertex_face_dir(i,(j)%3);

  //     // u = velocity(f)*mesh->vertex_face_dir(i,j%3);
  //     // v = velocity(f)*mesh->vertex_face_dir(i,j%3);
      
  //     u_vertex += u*u*mesh->vertex_face_area(i, j);
  //     v_vertex += v*v*mesh->vertex_face_area(i, j);

  //     // total_area += mesh->vertex_face_area(i, j);

  //     // std::cout<<velocity(f)<<' ';
  //   }

  //   // std::cout<<std::endl;

  //   u_vertex /= total_area;
  //   v_vertex /= total_area;

  //   // KE_vertex(i) = pow(u_vertex, 2.0) + pow(v_vertex, 2.0);
  //   KE_vertex(i) = u_vertex + v_vertex;

  //   // std::cout<<i<<' '<<u_vertex<<' '<<mesh->vertex_pos_sph(i,0)<<' '<<mesh->vertex_pos_sph(i,1)<<std::endl;
  //   // std::cout<<i<<' '<<pow(u_vertex,2.0) + pow(v_vertex,2.0)<<' '<<mesh->vertex_pos_sph(i,0)<<' '<<mesh->vertex_pos_sph(i,1)<<std::endl;
  // } 

  // for (i=0; i<face_num; i++) {
  //   int v1, v2;

  //   v1 = mesh->face_vertexes(i, 0);
  //   v2 = mesh->face_vertexes(i, 1);

  //   double KE = (KE_vertex(v1) + KE_vertex(v2))*0.5;

  //   // std::cout<<i<<' '<<KE<<' '<<mesh->vertex_pos_sph(i,0)<<' '<<mesh->vertex_pos_sph(i,1)<<std::endl;
  // }

  
};