#include "SEM.h"
#include <algorithm>
#include "geometrycentral/numerical/linear_solvers.h"
#include<Eigen/SparseCholesky>
/*
Implements the Stretch Energy Minimization described in "A Novel Stretch Energy Minimization Algorithm for
Equiareal Parameterizations" https://link.springer.com/article/10.1007/s10915-017-0414-y
*/

/*
handles initialization, 
step 1 to 4 in Algorithm 1
*/
SEM::SEM(ManifoldSurfaceMesh* surfaceMesh, VertexPositionGeometry* geo) {
	this->mesh = surfaceMesh;
	this->geometry = geo;
	this->ORIGINAL = this->geometry->inputVertexPositions;
	//scale the original mesh so it looks proportional to 1 (only for visualization purpose)
	double ratio = 2;
	for (size_t i = 0; i != ORIGINAL.size(); i++)
	{
		ORIGINAL[i] *= ratio;
	}
	// 1. compute total surface area using ORIGINAL (10x in vertex -> 100x in surface area);
	total_area = 0;
	area_vec.resize(mesh->nFaces());
	for (Face f : mesh->faces()) {
		double area = geometry->faceArea(f);
		area_vec[f.getIndex()] = area;
		total_area += area;
	}
	for (auto & i : area_vec) {
		i = i / total_area * PI;
	}
	// 2. initialize boundary mapping
	// 2.1. compute total length
	Halfedge he;
	for (Halfedge hi : mesh->halfedges()) {
		if (!hi.isInterior()) {
			he = hi; // exist at the first boundary halfedge
			break;
		}
	}
	Halfedge hi = he;
	double total_length = 0;
	do {
		bdry_vec.push_back(hi.tipVertex().getIndex());
		total_length += geometry->edgeLength(hi.edge());
		hi = hi.next();
	} while (hi != he);
	// 2.2. map boundary to (cos(theta), sin(theta))
	hi = he;
	DenseMatrix<double> fB(bdry_vec.size(), 2);
	double cur_length = 0;
	int i = 0;
	do {
		cur_length += geometry->edgeLength(hi.edge());
		double theta = cur_length / total_length * 2 * PI;
		fB(i, 0) = cos(theta);
		fB(i, 1) = sin(theta);
		//bdry_vertices(i, 2) = 0;
		hi = hi.next();
		i++;
	} while (hi != he);
	// 3. initialize interior mapping
	// 3.1 construct fI map
	for (int i = 0; i != mesh->nVertices(); i++) {
		if (find(bdry_vec.begin(), bdry_vec.end(), i) == bdry_vec.end()) {
			intr_vec.push_back(i);
		}
	}
	// 3.2 solve Laplacian to get fI
	geometry->requireCotanLaplacian();
	SparseMatrix<double> & L = geometry->cotanLaplacian;
	//std::cout << L << std::endl;
	DenseMatrix<double> fI = solveInteriorLaplacian(fB, L);
	// 4. update geometry
	updateGeometry(fB, fI);
	this->fB_save = fB;
	this->fI_save = fI;
}

/*
solve the equation LII * fI = -LIB * fB given L and fB, 
step 4 and 11
*/
DenseMatrix<double> SEM::solveInteriorLaplacian(DenseMatrix<double> fB, SparseMatrix<double> L)
{
	// 1. construct LII LIB using fB fI indices
	auto subL = subInteriorLaplacian(L);
	SparseMatrix<double> LII = subL.first;
	SparseMatrix<double> LIB = subL.second;
	// 2. solve linear system 
	DenseMatrix<double> rhs = -LIB * fB;
	Eigen::SimplicialLDLT<SparseMatrix<double>> solver(LII);
	DenseMatrix<double> fI = solver.solve(rhs);
	//std::cout << rhs;
	//std::cout << fI;
	return fI;
}
/*
solve the equation LBB * fB = -LBI * fI given L and fI, and centralize and normalize the boundary points
step 8,9,10
*/
DenseMatrix<double> SEM::solveBoundaryLaplacian(DenseMatrix<double> fI, SparseMatrix<double> L) {
	// solve
	auto subL = subBoundaryLaplacian(L);
	SparseMatrix<double> LBB = subL.first;
	SparseMatrix<double> LBI = subL.second;
	DenseMatrix<double> rhs = -LBI * fI;
	Eigen::SimplicialLDLT<SparseMatrix<double>> solver(LBB);
	DenseMatrix<double> fB = solver.solve(rhs);
	// centralize
	auto mean = fB.colwise().mean();
	fB.rowwise() -= mean;
	// uniform
	fB.rowwise().normalize();
	return fB;
};
/*
compute stretch factor and update Laplacian matrix
step 4,5 and a12,13
*/
SparseMatrix<double> SEM::computeLaplacian() {
	// using formula on page 7, equ(6) and following equations
	SparseMatrix<double> L(mesh->nVertices(), mesh->nVertices());
	std::vector<Eigen::Triplet<double>> tripletL;
	int i, j;
	auto w_ijk = [=](Halfedge hk, int i, int j) {
		if (hk.isInterior()) { // if vijk in F(M)
			double area_ijk = area_vec[hk.face().getIndex()];
			int k = hk.next().tipVertex().getIndex();
			Vector3 fi = geometry->vertexPositions[i];
			Vector3 fj = geometry->vertexPositions[j];
			Vector3 fk = geometry->vertexPositions[k];
			double denom = (fi.x - fk.x) * (fj.x - fk.x) + (fi.y - fk.y) * (fj.y - fk.y);
			return denom / area_ijk / 4; // outer half is also divided here
		}
		else return .0;
	};
	for (Edge e : mesh->edges()) {
		i = e.firstVertex().getIndex();
		j = e.secondVertex().getIndex();
		Halfedge hk = e.halfedge();
		double wij = w_ijk(hk, i, j) + w_ijk(hk.twin(), i, j);
		tripletL.push_back(Eigen::Triplet<double>{j, i, -wij});
		tripletL.push_back(Eigen::Triplet<double>{i, j, -wij});
		tripletL.push_back(Eigen::Triplet<double>{j, j, wij});
		tripletL.push_back(Eigen::Triplet<double>{i, i, wij});
	}
	L.setFromTriplets(tripletL.begin(), tripletL.end());
	//std::cout << L;
	return L;
}
/*
returns LII and LBI from a large L
*/
std::pair<SparseMatrix<double>, SparseMatrix<double>> SEM::subInteriorLaplacian(SparseMatrix<double> L) { 
	SparseMatrix<double> LII(intr_vec.size(), intr_vec.size());
	SparseMatrix<double> LIB(intr_vec.size(), bdry_vec.size());
	std::vector<Eigen::Triplet<double>> tripletLII;
	std::vector<Eigen::Triplet<double>> tripletLIB;
	for (size_t i = 0; i != intr_vec.size(); i++)
		for (size_t j = 0; j != bdry_vec.size(); j++)
			if ( L.coeff(intr_vec[i], bdry_vec[j]) != 0) 
				tripletLIB.push_back(Eigen::Triplet<double>{(int)i, (int)j, L.coeff(intr_vec[i], bdry_vec[j])});
	for (size_t i = 0; i != intr_vec.size(); i++)
		for (size_t j = 0; j != intr_vec.size(); j++)
			if ( L.coeff(intr_vec[i], intr_vec[j]) != 0) 
				tripletLII.push_back(Eigen::Triplet<double>{(int)i, (int)j, L.coeff(intr_vec[i], intr_vec[j])});
	LIB.setFromTriplets(tripletLIB.begin(), tripletLIB.end());
	LII.setFromTriplets(tripletLII.begin(), tripletLII.end());
	//std::cout << LII << std::endl << LIB << std::endl;
	return std::make_pair(LII, LIB);
};
/*
returns LBB and LIB from a large L
*/
std::pair<SparseMatrix<double>, SparseMatrix<double>> SEM::subBoundaryLaplacian(SparseMatrix<double> L) {
	SparseMatrix<double> LBB(bdry_vec.size(), bdry_vec.size());
	SparseMatrix<double> LBI(bdry_vec.size(), intr_vec.size());
	std::vector<Eigen::Triplet<double>> tripletLBB;
	std::vector<Eigen::Triplet<double>> tripletLBI;
	for (size_t i = 0; i != bdry_vec.size(); i++)
		for (size_t j = 0; j != intr_vec.size(); j++)
			if ( L.coeff(bdry_vec[i], intr_vec[j]) != 0) 
				tripletLBI.push_back(Eigen::Triplet<double>{(int)i, (int)j, L.coeff(bdry_vec[i], intr_vec[j])});
	for (size_t i = 0; i != bdry_vec.size(); i++)
		for (size_t j = 0; j != bdry_vec.size(); j++)
			if ( L.coeff(bdry_vec[i], bdry_vec[j]) != 0) 
				tripletLBB.push_back(Eigen::Triplet<double>{(int)i, (int)j, L.coeff(bdry_vec[i], bdry_vec[j])});
	LBB.setFromTriplets(tripletLBB.begin(), tripletLBB.end());
	LBI.setFromTriplets(tripletLBI.begin(), tripletLBI.end());
	return std::make_pair(LBB, LBI);
};
// a single step
double SEM::step() {
	auto L = computeLaplacian();
	auto fB = solveBoundaryLaplacian(this->fI_save, L);
	auto fI = solveInteriorLaplacian(fB, L);
	updateGeometry(fB, fI);
	return loss();
};
// compute accumulate difference in mesh using current geometry
double SEM::loss() {
	double total_loss = 0;
	for (Face f : mesh->faces()) {
		total_loss += abs(geometry->faceArea(f) - area_vec[f.getIndex()]);
	}
	return total_loss / mesh->nFaces();
};
/*
update and rewrite the geometry using fI and fB
*/
void SEM::updateGeometry(DenseMatrix<double> fB, DenseMatrix<double> fI) {
	for (size_t i = 0; i != bdry_vec.size(); i++) {
		geometry->vertexPositions[bdry_vec[i]] = Vector3{ fB(i, 0) + 2 ,fB(i, 1) ,0 };
	}
	for (size_t i = 0; i != intr_vec.size(); i++) {
		geometry->vertexPositions[intr_vec[i]] = Vector3{ fI(i, 0) + 2,fI(i, 1) ,0 };
	}
};

void SEM::debug() {
	// 1. return the boudary value for debug and register is as point clouds
	std::cout << "debug";
}