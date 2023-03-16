#pragma once
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include <chrono>
#include <set>
using namespace geometrycentral;
using namespace geometrycentral::surface;

#define USE_DELAUNAY

enum BOUNDARY {CIRCLE, SQUARE};

class SEM {
public:
	// member variables
	ManifoldSurfaceMesh* mesh;
	VertexPositionGeometry* geometry;
//#ifdef USE_DELAUNAY;
//	ManifoldSurfaceMesh mesh_delaunay;
//	VertexPositionGeometry geometry_delaunay;
//#endif
	BOUNDARY boundary_type;
	double param_area;
	double total_area;
	std::vector<int> bdry_vec; // stores the index mapping of fB index to actual vertex index
	std::vector<int> intr_vec; // stores the index mapping of fI index to actual vertex index
	std::vector<double> area_vec; // stores the area of face by face index
	std::vector<double> dual_area_vec; // stores the dual area of vertices
	std::vector<int> sub_vec; // shows the submatrix index in interior or boundary. need to be used with .isBoundary()
	VertexData<Vector3> ORIGINAL; // original geometry data
	VertexData<Vector3> INITIAL;
	DenseMatrix<double> ORIGINAL_MAT;
	DenseMatrix<double> fB_save;
	DenseMatrix<double> fI_save;
	// member functions
	SEM(ManifoldSurfaceMesh* surfaceMesh, VertexPositionGeometry* geo, BOUNDARY boundary_type=CIRCLE); // initialize the parameterization
	DenseMatrix<double> solveInteriorLaplacian(DenseMatrix<double> fB, SparseMatrix<double> LII, SparseMatrix<double> LIB); //use fB and L to get fI
	DenseMatrix<double> solveBoundaryLaplacian(DenseMatrix<double> fI, SparseMatrix<double> LBB, SparseMatrix<double> LBI); //use fI and L to get fB
	std::pair<SparseMatrix<double>, SparseMatrix<double>> subInteriorLaplacian(SparseMatrix<double> L); // returns LII and LIB from a given L
	std::pair<SparseMatrix<double>, SparseMatrix<double>> subBoundaryLaplacian(SparseMatrix<double> L); // returns LBB and LBI from a given L
	SparseMatrix<double> computeLaplacian(); // computes new Laplacian based on parameterization;
	double stepLength(DenseMatrix<double>& grad_v); // compute the current step length
	void computeLaplacianSub(SparseMatrix<double> & LII, SparseMatrix<double> & LBB, SparseMatrix<double> & LIB); // compute sub Laplacian matrix and assemble at the same time
	void updateGeometry(DenseMatrix<double> fB, DenseMatrix<double> fI); // will rewrite vertex position in the geometry member
	void debug(); // debug
	double step();
	double loss(); // using current geometry to compute sigma
	double lossDual();
	void delaunayFlip();

	// for density equalizing gradient flow
	// member 
	SparseMatrix<double> W_fv;
	DenseMatrix<double> rho_v;
	std::shared_ptr<SignpostIntrinsicTriangulation> intTri;
	// function
	DenseMatrix<double> computeRhoV(); // computes RhoV
	std::pair<DenseMatrix<double>, DenseMatrix<double>> gradFlowDebug();
	DenseMatrix<double> gradFlow(double dt); // perform the grad flow of RhoV using previous RhoV
	SparseMatrix<double> computeWfv();
	SparseMatrix<double> computeLaplacianDEM();
	double stepDEMFixBoundary(double dt);
	void computeSEMAll(double eps);
};