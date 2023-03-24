#ifndef MESH_H
#define MESH_H

#include "globals.h"
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include "mathRoutines.h"
#include "sphericalHarmonics.h"

//#include <mkl.h>
//#include <mkl_spblas.h>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>

typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMat;
typedef Eigen::Map< SpMat > SpMatMap;
typedef Eigen::Triplet<double> Triplet;
// typedef Eigen::Vector<doub

typedef Eigen::Matrix<double,1,Eigen::Dynamic> DenVector;
typedef Eigen::Map<DenVector> MapVector;

typedef Eigen::Matrix<double, Eigen::Dynamic, 6, Eigen::RowMajor> Vandermonde;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> DenseMat;

class Globals;

//Class stores all coordinate information about the mesh
class Mesh {
private:
  int ReadMeshFile(void);
  int ReadLatLonFile(void);
  int ReadWeightingFile(void);
  int ReadLeastSquaresFile(void);
  int CalcMappingCoords(void);
  int CalcVelocityTransformFactors(void);
  int CalcControlVolumeEdgeLengths(void);
  int CalcControlVolumeEdgeCentres(void);
  int CalcControlVolumeEdgeNormals(void);
  int CalcControlVolumeArea(void);
  int CalcControlVolumeMass(void);
  int CalcControlVolumeVertexR(void);
  int CalcElementAreas(void);
  int CalcCentNodeDists(void);
  int CalcTrigFunctions(void);
  int CalcNodeDists(void);
  int CalcMaxTimeStep(void);
  int CalcControlVolumeInterpMatrix(void);

  int AssignFaces(void);

  int DefineBoundaryCells(void);

  int CalcLegendreFuncs();
  int CalcGradOperatorCoeffs();
  int CalcDivOperatorCoeffs();
  int CalcLaplaceOperatorCoeffs(void);
  int CalcCoriolisOperatorCoeffs(void);
  int CalcLinearDragOperatorCoeffs(void);
  int CalcCurlOperatorCoeffs(void);
  int CalcAdjacencyMatrix(void);
  int CalcRBFInterpMatrix(void);
  int CalcRBFInterpMatrix2(void);
  int CalcFreeSurfaceSolver(void);
  int GeneratePressureSolver(void);
  int GenerateMomentumOperator(void);

public:
//	Mesh(); //constructor
  Mesh(Globals *);
  Mesh(Globals *, int, int, int, int, int);

  std::ostringstream outstring;

  int node_num;
  int face_num;
  int vertex_num;

  // Each row in node_pos_sph represents the lat,lon coords of node. Node ID 5
  // corresponds to node_pos_sph(5,x).
  Array2D<double> node_pos_sph;

  // Same as node_pos_sph but in mapped coords
  Array3D<double> node_pos_map;

  // Each row in node_friends contains the ID (index)
  // to all surrounding friends. E.g, the ID's of the neighbouring nodes to
  // node 4 would be accessed via node_friends[4][0], node_friends[4][1], etc.
  Array2D<int> node_friends;

  Array2D<double> centroid_node_dists_map;

  Array2D<double> node_dists;

  // Each row contains the coordinates of the centroids to the 5 or 6 surrounding
  // cells for a specific node (e.g., node ID 8 corresponds to row 8).
  // Column 0 holds the coordinates for the centroid formed from the parent
  // node, friend 0, and friend 1. Column 1 holds the coordinates for the centroid
  // formed from the parent node, friend 1, and friend 2, etc.
  Array3D<double> centroid_pos_sph;
  Array3D<double> centroid_pos_map;

  // Array to hold mapping factor at at each node and its friends
  Array2D<double> node_m;

  // Array to hold mapping factor at at each centroid
  Array2D<double> centroid_m;

  // Array to hold the cos(alpha) and sin(alpha) velocity conversion coefficients
  // from Lee and Macdonald (2009)
  Array3D<double> node_vel_trans;

  Array2D<double> node_node_RBF;

  // Array to store the edge length of each side to the control volume surrounding
  // the central node
  Array2D<double> control_vol_edge_len;

  // Array to store the midpoint of each side to the control volume surrounding
  // the central node
  Array3D<double> control_vol_edge_centre_pos_map;
  Array3D<double> control_vol_edge_centre_pos_sph;
  // Array3D<double> node_face_normal_vec_map;


  // Array to store the midpoint mapping factor of each side to the control
  // volume surrounding the central node
  Array2D<double> control_vol_edge_centre_m;

  // Array to store the midpoint normal vectors to each side of the control
  // volume surroinding the central node
  Array3D<double> control_vol_edge_normal_map;

  // Array to store the area of the control volume for each node
  Array1D<double> control_volume_surf_area_map;

  Array1D<double> control_volume_mass;

  // Array to store the triangular areas within each subelement.
  Array3D<double> node_friend_element_areas_map;

  Array1D<int> cell_is_boundary;

  // Arrays to store trig functions evaluated at each node
  Array2D<double> trigLat;
  Array2D<double> trigLon;
  Array2D<double> trig2Lat;
  Array2D<double> trig2Lon;
  Array2D<double> trigSqLat;
  Array2D<double> trigSqLon;

  Array2D<int> faces;
  Array2D<int> face_nodes;
  Array2D<int> face_vertexes;
  Array1D<double> face_area;
  // Array2D<int> face_friends;
  Array2D<int> face_dir;
  Array2D<double> face_normal_vec_map;
  Array2D<double> face_normal_vec_xyz;
  Array1D<double> face_node_dist;
  Array1D<double> face_node_dist_r;
  Array2D<double> node_face_arc; // angular distance between node and face center
  Array2D<double> node_face_RBF;
  Array2D<double> face_centre_pos_sph;
  Array2D<double> face_intercept_pos_sph;
  Array2D<double> face_centre_m;
  Array3D<double> face_node_pos_map;
  Array1D<double> face_len;
  Array2D<int> node_face_dir;
  Array3D<double> face_node_vel_trans;
  Array1D<int> face_is_boundary;

  Array2D<int> vertexes;
  Array2D<double> vertex_pos_sph;
  Array2D<double> vertex_R;
  Array2D<int> vertex_nodes;
  Array2D<double> node_R;
  Array2D<int> vertex_faces;
  Array2D<int> vertex_face_dir;
  Array1D<double> vertex_area;
  Array1D<double> vertex_area_r;
  Array1D<double> vertex_sinlat;

  Array2D<int> face_interp_friends;
  Array2D<double> face_interp_weights;

  // ArrayD<int> * faces;

  Array3D<double> trigMLon;
  Array2D<double> sh_matrix;
  Array1D<double> sh_matrix_fort;

  Array1D<double> Pbar_20;
  Array1D<double> Pbar_22;

  int * interpRows;
  int * interpCols;
  double * interpWeights;

  // Handle to the Intel MKL direct sparse matrix solver to solve an elliptical
  // equation in pressure correction
  // _MKL_DSS_HANDLE_t pressureSolverHandle;

  // sparse_matrix_t * interpMatrix;

  SpMat operatorLinearDrag;
  SpMat operatorCoriolis;
  SpMat operatorDivergence;
  SpMat operatorGradient;
  SpMat operatorCurl;
  SpMat node2faceAdj;
  SpMat node2nodeAdj;
  SpMat operatorRBFinterp;
  // SpMat rbfDeriv2;
  SpMat operatorSecondDeriv;
  SpMat operatorDirectionalSecondDeriv;

  SpMat interpMatrix;
  SpMat vandermondeInv;

  // Eigen::ColPivHouseholderQR<Eigen::MatrixXd> * interpSolvers;
  Eigen::FullPivLU<Eigen::MatrixXd> * interpSolvers;
  // Vandermonde * interpVandermonde;
  std::vector<Eigen::MatrixXd> interpVandermonde;
  std::vector<Eigen::MatrixXd> interpVandermondeScalar;
  Eigen::MatrixXd * interpSolversInv;
  // Eigen::ColPivHouseholderQR<Eigen::MatrixXd> * interpSolvers2;
  Eigen::LLT<Eigen::MatrixXd> * interpSolvers2;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> cg;
  // Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> cg;


  // Eigen::Matrix5d * interpSolvers;

  // sparse_matrix_t * operatorMomentum;
  // sparse_matrix_t * operatorLaplacian;
  // sparse_matrix_t * operatorLaplacian2;
  // sparse_matrix_t * operatorGradientX;
  // sparse_matrix_t * operatorGradientY;
  // sparse_matrix_t * operatorGradient;
  // sparse_matrix_t * operatorDivergence;
  // sparse_matrix_t * operatorDivergenceX;
  // sparse_matrix_t * operatorDivergenceY;
  // sparse_matrix_t * operatorCoriolis;
  // sparse_matrix_t * operatorLinearDrag;

  Globals * globals;

};

#endif
