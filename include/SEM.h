#pragma once
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include <set>
using namespace geometrycentral;
using namespace geometrycentral::surface;

class SEM {
public:
	// member variables
	ManifoldSurfaceMesh* mesh;
	VertexPositionGeometry* geometry;
	double total_area;
	std::vector<int> bdry_vec; // stores the index mapping of fB index to actial vertex index
	std::vector<int> intr_vec; // stores the index mapping of fI index to actial vertex index
	std::vector<double> area_vec; // stores the area of face by face index
	VertexData<Vector3> ORIGINAL; // original geometry data
	DenseMatrix<double> fB_save;
	DenseMatrix<double> fI_save;
	// member functions
	SEM(ManifoldSurfaceMesh* surfaceMesh, VertexPositionGeometry* geo); // initialize the parameterization
	DenseMatrix<double> solveInteriorLaplacian(DenseMatrix<double> fB, SparseMatrix<double> L); //use fB and L to get fI
	DenseMatrix<double> solveBoundaryLaplacian(DenseMatrix<double> fI, SparseMatrix<double> L); //use fI and L to get fB
	std::pair<SparseMatrix<double>, SparseMatrix<double>> subInteriorLaplacian(SparseMatrix<double> L); // returns LII and LIB from a given L
	std::pair<SparseMatrix<double>, SparseMatrix<double>> subBoundaryLaplacian(SparseMatrix<double> L); // returns LBB and LBI from a given L
	SparseMatrix<double> computeLaplacian(); // computes new Laplacian based on parameterization;
	void updateGeometry(DenseMatrix<double> fB, DenseMatrix<double> fI); // will rewrite vertex position in the geometry member
	void debug(); // debug
	double step();
	double loss(); // using current geometry to compute sigma
};