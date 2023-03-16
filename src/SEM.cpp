#include "SEM.h"
#include <algorithm>
#include "geometrycentral/numerical/linear_solvers.h"
#include<Eigen/SparseCholesky>
#include<Eigen/SparseLU>

using namespace std;

/*
Implements the Stretch Energy Minimization described in "A Novel Stretch Energy Minimization Algorithm for
Equiareal Parameterizations" https://link.springer.com/article/10.1007/s10915-017-0414-y
*/

/*
handles initialization, 
step 1 to 4 in Algorithm 1
*/
SEM::SEM(ManifoldSurfaceMesh* surfaceMesh, VertexPositionGeometry* geo, BOUNDARY boundary_type) {
	this->mesh = surfaceMesh;
	this->geometry = geo;
	this->ORIGINAL = this->geometry->inputVertexPositions;
	this->boundary_type = boundary_type;
	rho_v.resize(mesh->nVertices(),1 );
	sub_vec.resize(mesh->nVertices());
	//scale the original mesh so it looks proportional to 1 (only for visualization purpose)
	double ratio = 1;
	for (size_t i = 0; i != ORIGINAL.size(); i++)
	{
		ORIGINAL[i] *= ratio;
	}
	if (boundary_type == CIRCLE)
		param_area = PI;
	else if (boundary_type == SQUARE)
		param_area = 1;
	// 1. compute total surface area using ORIGINAL (10x in vertex -> 100x in surface area);
	total_area = 0;
	area_vec.resize(mesh->nFaces());
	for (Face f : mesh->faces()) {
		double area = geometry->faceArea(f);
		area_vec[f.getIndex()] = area;
		total_area += area;
	}
	for (auto & i : area_vec) {
		i = i / total_area * param_area;
	}
	dual_area_vec.resize(mesh->nVertices());
	for (Vertex v : mesh->vertices()) {
		double area = geometry->vertexDualArea(v);
		dual_area_vec[v.getIndex()] = area;
	}
	for (auto & i : dual_area_vec) {
		i = i / total_area * param_area;
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
	int count_bdry = 0;
	do {
		size_t cur_idx = hi.tipVertex().getIndex();
		sub_vec[cur_idx] = count_bdry++;
		bdry_vec.push_back(cur_idx);
		total_length += geometry->edgeLength(hi.edge());
		hi = hi.next();
	} while (hi != he);
	// 2.2. map boundary to (cos(theta), sin(theta))
	hi = he;
	DenseMatrix<double> fB(bdry_vec.size(), 2);
	double cur_length = 0;
	int i = 0;
	if (boundary_type == CIRCLE) {
		do {  
			cur_length += geometry->edgeLength(hi.edge());
			double theta = cur_length / total_length * 2 * PI;
			fB(i, 0) = cos(theta);
			fB(i, 1) = sin(theta);
			//bdry_vertices(i, 2) = 0;
			hi = hi.next();
			i++;
		} while (hi != he);
	}
	else if (boundary_type == SQUARE) {
		do {
			cur_length += geometry->edgeLength(hi.edge());
			double ratio = cur_length / total_length;
			if (ratio < 0.25) {
				fB(i, 0) = ratio * 4 - 0.5;
				fB(i, 1) = -0.5;
			} else if (ratio < 0.5) {
				fB(i, 0) = 0.5;
				fB(i, 1) = (ratio - 0.25) * 4 - 0.5;
			} else if (ratio < 0.75) {
				fB(i, 0) = 0.5  - (ratio - 0.5) * 4;
				fB(i, 1) = 0.5;
			} else {
				fB(i, 0) = - 0.5 ;
				fB(i, 1) = 0.5 - (ratio - 0.75) * 4;
			}
			hi = hi.next();
			i++;
		} while (hi != he);
	}
	// 3. initialize interior mapping
	// 3.1 construct fI map
	
	int count_intr = 0;
	for (int i = 0; i != mesh->nVertices(); i++) {
		if (find(bdry_vec.begin(), bdry_vec.end(), i) == bdry_vec.end()) {
			intr_vec.push_back(i);
			sub_vec[i] = count_intr++;
		}
	}
	// 3.2 solve Laplacian to get fI
	geometry->requireCotanLaplacian();
	auto subL = subInteriorLaplacian(geometry->cotanLaplacian);
	//SparseMatrix<double> & L = geometry->cotanLaplacian;
	SparseMatrix<double> LII(intr_vec.size(), intr_vec.size());
	SparseMatrix<double> LBB(bdry_vec.size(), bdry_vec.size());
	SparseMatrix<double> LIB(intr_vec.size(), bdry_vec.size());
	computeLaplacianSub(LII, LBB, LIB);
	auto fI = solveInteriorLaplacian(fB, LII, LIB);

	// 4. update geometry
	updateGeometry(fB, fI);
	this->fB_save = fB;
	this->fI_save = fI;
	this->INITIAL = geometry->inputVertexPositions;

	W_fv.resize(mesh->nVertices(), mesh->nFaces());
	std::vector<Eigen::Triplet<double>> tripletW_fv;
	for (Vertex v : mesh->vertices())
	{
		double v_total_area = 0;
		for (Face f : v.adjacentFaces())
			v_total_area += area_vec[f.getIndex()];
		for (Face f : v.adjacentFaces())
			tripletW_fv.push_back({ (int)v.getIndex(), (int)f.getIndex(), area_vec[f.getIndex()] / v_total_area });
	}
	W_fv.setFromTriplets(tripletW_fv.begin(), tripletW_fv.end());
#ifdef USE_DELAUNAY 
	delaunayFlip();
#endif
}

void SEM::delaunayFlip() {
	// overwrite the mesh and geometry ptr to a newly consturcted intrinsic mesh stored inside
	// totally construct a new mesh here because intTri->intrinsicMesh is a unique_ptr
	// this part is very very ugly but now complies with everything
	intTri.reset(new SignpostIntrinsicTriangulation(*mesh, *geometry));
	intTri->flipToDelaunay();
	const std::vector<std::vector<size_t>> intFaces = intTri->intrinsicMesh->getFaceVertexList();
	std::vector<Vector3> intVertices;
	for (int i=0; i!=mesh->nVertices(); i++)
		intVertices.push_back(geometry->vertexPositions[i]);
	SimplePolygonMesh simpleMesh( intFaces, intVertices);
	auto lvals = makeManifoldSurfaceMeshAndGeometry(simpleMesh.polygons, simpleMesh.vertexCoordinates);
	std::unique_ptr<ManifoldSurfaceMesh> mesh_new;
	std::unique_ptr<VertexPositionGeometry> geometry_new;
	std::tie(mesh_new, geometry_new) = std::tuple < std::unique_ptr<ManifoldSurfaceMesh>,
		std::unique_ptr<VertexPositionGeometry>>(std::move(std::get<0>(lvals)),  // mesh
			std::move(std::get<1>(lvals)));
	mesh = mesh_new.release();
	geometry = geometry_new.release();
}

/*
solve the equation LII * fI = -LIB * fB given L and fB, 
step 4 and 11
*/
DenseMatrix<double> SEM::solveInteriorLaplacian(DenseMatrix<double> fB, SparseMatrix<double> LII, SparseMatrix<double> LIB)
{
	// 1. construct LII LIB using fB fI indices
	//auto subL = subInteriorLaplacian(L);
	//SparseMatrix<double> LII = subL.first;
	//SparseMatrix<double> LIB = subL.second;
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
DenseMatrix<double> SEM::solveBoundaryLaplacian(DenseMatrix<double> fI, SparseMatrix<double> LBB, SparseMatrix<double> LBI) {
	// solve
	// auto subL = subBoundaryLaplacian(L);
	//SparseMatrix<double> LBB = subL.first;
	//SparseMatrix<double> LBI = subL.second;
	DenseMatrix<double> rhs = -LBI * fI;
	Eigen::SimplicialLDLT<SparseMatrix<double>> solver(LBB);
	DenseMatrix<double> fB = solver.solve(rhs);
	
	if (boundary_type == CIRCLE) {
		// centralize & normalize
		auto mean = fB.colwise().mean();
		fB.rowwise() -= mean;
		fB.rowwise().normalize();
	}
	else if (boundary_type == SQUARE) {
		for (int i = 0; i != fB.rows(); i++)
		{
			double rad = atan2(fB(i, 1), fB(i, 0));
			if (rad >= - M_PI_4 && rad < M_PI_4) {
				fB(i, 0) = 0.5;
				fB(i, 1) = 0.5 * tan(rad);
			}
			else if (rad >= M_PI_4 && rad < 3 * M_PI_4) {
				fB(i, 0) = 0.5 / tan(rad);
				fB(i, 1) = 0.5;
			}
			else if (rad >= 3 * M_PI_4 || rad < -3 * M_PI_4) {
				fB(i, 0) = -0.5;
				fB(i, 1) = tan(rad) * -0.5;
			}
			else if (rad >= -3 * M_PI_4 && rad < - M_PI_4) {
				fB(i, 0) = -0.5 / tan(rad);
				fB(i, 1) = -0.5;
			}
		}
	}
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
			double nume = (fi.x - fk.x) * (fj.x - fk.x) + (fi.y - fk.y) * (fj.y - fk.y);
			return nume / area_ijk / 4; // outer half is also divided here
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

void SEM::computeLaplacianSub(SparseMatrix<double> & LII, SparseMatrix<double> & LBB, SparseMatrix<double> & LIB) {
	// SparseMatrix<double> LII(intr_vec.size(), intr_vec.size());
	// SparseMatrix<double> LBB(bdry_vec.size(), bdry_vec.size());
	// SparseMatrix<double> LIB(intr_vec.size(), bdry_vec.size());
	std::vector<Eigen::Triplet<double>> tripletLII;
	std::vector<Eigen::Triplet<double>> tripletLBB;
	std::vector<Eigen::Triplet<double>> tripletLIB;

	auto w_ijk = [=](Halfedge hk, int i, int j) {
		if (hk.isInterior()) { // if vijk in F(M)
			double area_ijk = area_vec[hk.face().getIndex()];
			int k = hk.next().tipVertex().getIndex();
			Vector3 fi = geometry->vertexPositions[i];
			Vector3 fj = geometry->vertexPositions[j];
			Vector3 fk = geometry->vertexPositions[k];
			double nume = (fi.x - fk.x) * (fj.x - fk.x) + (fi.y - fk.y) * (fj.y - fk.y);
			return nume / area_ijk / 4; // outer half is also divided here
		}
		else return .0;
	};
	int i, j;
	int sub_i, sub_j;
	for (Edge e : mesh->edges()) {
		i = e.firstVertex().getIndex();
		j = e.secondVertex().getIndex();
		Halfedge hk = e.halfedge();
		double wij = w_ijk(hk, i, j) + w_ijk(hk.twin(), i, j);
		sub_i = sub_vec[i];
		sub_j = sub_vec[j];
		if (mesh->vertex(i).isBoundary() && mesh->vertex(j).isBoundary()) // boundary
		{
			tripletLBB.push_back(Eigen::Triplet<double>{sub_j, sub_i, -wij});
			tripletLBB.push_back(Eigen::Triplet<double>{sub_i, sub_j, -wij});
			tripletLBB.push_back(Eigen::Triplet<double>{sub_j, sub_j, wij});
			tripletLBB.push_back(Eigen::Triplet<double>{sub_i, sub_i, wij});
		}
		else if (!mesh->vertex(i).isBoundary() && !mesh->vertex(j).isBoundary()) // interior
		{
			tripletLII.push_back(Eigen::Triplet<double>{sub_j, sub_i, -wij});
			tripletLII.push_back(Eigen::Triplet<double>{sub_i, sub_j, -wij});
			tripletLII.push_back(Eigen::Triplet<double>{sub_j, sub_j, wij});
			tripletLII.push_back(Eigen::Triplet<double>{sub_i, sub_i, wij});
		}
		else if (!mesh->vertex(i).isBoundary() && mesh->vertex(j).isBoundary()) // one side for LIB
		{
			tripletLIB.push_back(Eigen::Triplet<double>{sub_i, sub_j, -wij});
			tripletLBB.push_back(Eigen::Triplet<double>{sub_j, sub_j, wij}); // count for symmetry
			tripletLII.push_back(Eigen::Triplet<double>{sub_i, sub_i, wij});
		}
		else if (mesh->vertex(i).isBoundary() && !mesh->vertex(j).isBoundary()) // the other side for LIB
		{
			tripletLIB.push_back(Eigen::Triplet<double>{sub_j, sub_i, -wij});
			tripletLBB.push_back(Eigen::Triplet<double>{sub_i, sub_i, wij}); // count for symmetry
			tripletLII.push_back(Eigen::Triplet<double>{sub_j, sub_j, wij});
		}
	}
	LIB.setFromTriplets(tripletLIB.begin(), tripletLIB.end());
	LII.setFromTriplets(tripletLII.begin(), tripletLII.end());
	LBB.setFromTriplets(tripletLBB.begin(), tripletLBB.end());
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
	//auto L = computeLaplacian();
	//auto subL = subBoundaryLaplacian(L);
	//std::cout << "\n" << subL.first << "\n";
	//auto fB = solveBoundaryLaplacian(this->fI_save, L);
	//auto fI = solveInteriorLaplacian(fB, L);
	SparseMatrix<double> LII(intr_vec.size(), intr_vec.size());
	SparseMatrix<double> LBB(bdry_vec.size(), bdry_vec.size());
	SparseMatrix<double> LIB(intr_vec.size(), bdry_vec.size());
	computeLaplacianSub(LII, LBB, LIB);
	auto fB = solveBoundaryLaplacian(fI_save, LBB, LIB.transpose());
	auto fI = solveInteriorLaplacian(fB, LII, LIB);
	//std::cout << "\n" << LBB << "\n";
	updateGeometry(fB, fI);
	return loss();
};

// compute accumulate difference in mesh using current geometry
double SEM::loss() {
	double total_loss = 0;
	for (Face f : mesh->faces()) {
		total_loss += abs(geometry->faceArea(f) - area_vec[f.getIndex()]);
	}
	return total_loss / param_area;
};

// same loss, but using dual area
double SEM::lossDual() {
	double total_loss = 0;
	for (Vertex v : mesh->vertices()) {
		total_loss += abs(geometry->vertexDualArea(v) - dual_area_vec[v.getIndex()]);
	}
	return total_loss / param_area;
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

SparseMatrix<double> SEM::computeWfv() {
	auto triArea = [&](Face f)
	{
		Vector3 v[3];
		int i = 0;
		for (auto vi : f.adjacentVertices())
			v[i++] = geometry->vertexPositions[vi];
		Vector3 v1 = v[1] - v[0];
		Vector3 v2 = v[2] - v[1];
		Vector3 n = cross(v1, v2);
		double sign = n.z < 0 ? 1 : -1;
		//double sign = 1;
		return n.norm()/2 * sign;
	};
	std::vector<Eigen::Triplet<double>> tripletW_fv;
	for (Vertex v : mesh->vertices())
	{
		double v_total_area = 0;
		for (Face f : v.adjacentFaces())
			v_total_area += abs(triArea(f));
		for (Face f : v.adjacentFaces())
			tripletW_fv.push_back({ (int)v.getIndex(), (int)f.getIndex(), triArea(f) / v_total_area });
	}
	W_fv.setFromTriplets(tripletW_fv.begin(), tripletW_fv.end());
	return W_fv;
};

DenseMatrix<double> SEM::computeRhoV() {
	// compute RhoF and left multiply by Wfv
	// DenseMatrix<double> RhoF(mesh->nFaces(), 1);
	auto triArea = [&](Face f)
	{
		Vector3 v[3];
		int i = 0;
		for (auto vi : f.adjacentVertices())
			v[i++] = geometry->vertexPositions[vi];
		Vector3 v1 = v[1] - v[0];
		Vector3 v2 = v[2] - v[0];
		return cross(v1, v2).norm()/2;
	};
	//// might need scaling on small mesh
	//for (Face f: mesh->faces())
	//	RhoF(f.getIndex(), 0) = triArea(f) / area_vec[f.getIndex()];
	geometry->refreshQuantities();
	geometry->requireFaceAreas();
	//RhoF = geometry->faceAreas.toVector();
	//for (int i = 0; i != area_vec.size(); i++)
	//	RhoF(i, 0) = log(RhoF(i, 0) / area_vec[i] ) ; // log scale factor
	for (Vertex v : mesh->vertices()) {
		//double total_area = 0;
		//double total_area_new = 0;
		//for (Face f : v.adjacentFaces())
		//	total_area_new += geometry->faceAreas[f];
		//rho_v(v.getIndex(),0) = log(total_area_new / total_area);
		rho_v(v.getIndex(),0) = log(geometry->vertexDualArea(v) / dual_area_vec[v.getIndex()]);
	}
	//rho_v = W_fv * RhoF;
	// smoothing with implicit step! 
	geometry->requireCotanLaplacian();
	geometry->requireVertexGalerkinMassMatrix();
	auto rho_v_smooth = rho_v;
	SparseMatrix<double> id(mesh->nVertices(), mesh->nVertices());
	id.setIdentity();
	Eigen::SimplicialLDLT<SparseMatrix<double>> solver(geometry->vertexGalerkinMassMatrix + 0.1 * geometry->cotanLaplacian);
	for (int i = 0; i != 1; i++)
		rho_v_smooth = solver.solve(geometry->vertexGalerkinMassMatrix * rho_v);
	rho_v = rho_v_smooth;
	return rho_v_smooth;
};

double SEM::stepLength(DenseMatrix<double>& grad_v) {
	// for a single triangle single grad 
	auto stepLengthPerTri = [&](Face tri){
		// we replace the indexing of subscription using 1->i 2->j 3->k as the example case and rotate 
	    // (actually this step is redundant but I will leave it there since the whole thing is working)
		int i = 0, j = 2, k = 1;
		Vertex v[3];
		int ii = 0;
		for (auto vi : tri.adjacentVertices())
			v[ii++] = vi;
		double xi = geometry->vertexPositions[v[i]].x, yi = geometry->vertexPositions[v[i]].y,
			xj = geometry->vertexPositions[v[j]].x, yj = geometry->vertexPositions[v[j]].y,
			xk = geometry->vertexPositions[v[k]].x, yk = geometry->vertexPositions[v[k]].y;
		double h[3][2]; // subscript comes first, superscript comes second, for superscript [1,2] in paper is [0,1] here
		h[0][0] = grad_v(v[0].getIndex(), 0);
		h[0][1] = grad_v(v[0].getIndex(), 1);
		h[1][0] = grad_v(v[1].getIndex(), 0);
		h[1][1] = grad_v(v[1].getIndex(), 1);
		h[2][0] = grad_v(v[2].getIndex(), 0);
		h[2][1] = grad_v(v[2].getIndex(), 1);

		// solve the quadratic equation when three point becomes collinear
		double A = ((-h[i][1] + h[j][1])*(-h[i][0] + h[k][0]) + (h[i][0] - h[j][0])*(-h[i][1] + h[k][1]));
		double B = (-h[j][1] + h[k][1])*xi + (h[i][1] - h[k][1])*xj + (-h[i][1] + h[j][1])*xk + (h[j][0] - h[k][0])*yi + (-h[i][0] + h[k][0])*yj + (h[i][0] - h[j][0])*yk;
		double C = xj * yi - xk * yi - xi * yj + xk * yj + xi * yk - xj * yk;
		double Delta = B*B-4*A*C;
	
		// return smallest step size
		if (A == 0)
		{
			double temp = -C / B;
			if (temp < 0)
				return 1.0;
			else
				return temp;
		}
		if (Delta < 0) // no real solution:
			return 1.0; // maximum steplength
		else
		{ // take the smallest positive value
		  // geometrically it must be one pos one neg solution for this equation

			double temp1 = (-B + sqrt(Delta)) / 2 / A;
			double temp2 = (-B - sqrt(Delta)) / 2 / A; // is two positive root actually possible??
			if (temp1 < 0)
				temp1 = 1;
			if (temp2 < 0)
				temp2 = 1;
			return min(temp1, temp2); // return the smallest positive root while in 1, if both are negative means we can choose any positive step length
		}
	};

	int face_min = 0;
	double steplength = 1;
	for (auto f : mesh->faces()) {
		if (steplength > stepLengthPerTri(f)) // for debug
			face_min = f.getIndex();
		steplength = min(steplength, stepLengthPerTri(f));
	}
	std::cout << "steplength min at face:" << face_min << "\n";
	return steplength;
};


DenseMatrix<double> SEM::gradFlow(double dt)
{
	//DenseMatrix<double> SEM::gradFlow(double dt) {
	// solve rho_n with rho_n-1
	// setup
	DenseMatrix<double> grad_f(mesh->nFaces(), 3); 
	geometry->requireCotanLaplacian();
	geometry->requireVertexLumpedMassMatrix();
	SparseMatrix<double>& L = geometry->cotanLaplacian;
	SparseMatrix<double>& A = geometry->vertexLumpedMassMatrix;
	//double dt = std::min(rho_v.minCoeff() / rho_v.mean(), rho_v.mean() / rho_v.maxCoeff()) * M_PI;
	// compute rho_vn
	// Eigen::SimplicialLDLT<SparseMatrix<double>> solver(A+dt*L);
	// DenseMatrix<double> rho_vn = solver.solve(A * rho_v);
	// compute grad rho_fn

	for (Face f : mesh->faces())
	{
		Vector3 v[3];
		int v_idx[3];
		int n = 0;
		for (auto vi : f.adjacentVertices())
		{
			v_idx[n] = vi.getIndex();
			v[n++] = geometry->vertexPositions[vi];
		}
		// i: 0; j: 1; k:2
		Vector3 vjk = v[2] - v[1];
		Vector3 vki = v[0] - v[2];
		Vector3 vij = v[1] - v[0];
		Vector3 N = cross(vjk, vki);
		double area = N.norm() / 2;
		Vector3 temp = rho_v(v_idx[0]) * vjk + rho_v(v_idx[1]) * vki + rho_v(v_idx[2]) * vij;
		Vector3 grad = 0.5 / area * cross(N.unit(), rho_v(v_idx[0]) * vjk + rho_v(v_idx[1]) * vki + rho_v(v_idx[2]) * vij);
		grad_f(f.getIndex(), 0) = grad.x;
		grad_f(f.getIndex(), 1) = grad.y;
		grad_f(f.getIndex(), 2) = grad.z;

	}
	computeWfv(); // with direction correction(?), this line recompute W_fv internally
	DenseMatrix<double> grad_v =  W_fv * grad_f;
	// projection on boundary to be tangential
	Vector3 center{ 2, 0, 0 }; //
	for (auto i : bdry_vec)
	{
		Vector3 n = geometry->vertexPositions[i] - center;// normal vector
		Vector3 grad_cur{ grad_v(i, 0), grad_v(i, 1), grad_v(i, 2) };
		grad_cur -= dot(grad_cur, n) * n;
		grad_v(i, 0) = grad_cur.x;
		grad_v(i, 1) = grad_cur.y;
		grad_v(i, 2) = grad_cur.z;
	}
	for (Vertex v : mesh->vertices())
	{
		int i = v.getIndex();
		Vector3 grad_v_cur = Vector3{ grad_v(i,0), grad_v(i, 1), grad_v(i, 2) };
		geometry->vertexPositions[v] = geometry->vertexPositions[v] + dt * grad_v_cur; // the sign depends on 
																					   //geometry->vertexPositions[v] = geometry->vertexPositions[v] + grad_v_cur;
	}
	for (auto i : bdry_vec) // projection after advection
	{
		geometry->vertexPositions[i] = (geometry->vertexPositions[i] - center).unit() + center;
	}
	// rho_v = rho_vn;


	//center /= bdry_vec.size();a
	//for (auto i: bdry_vec)
	//	geometry->vertexPositions[i] = (geometry->vertexPositions[i]-center).unit;
	return grad_v;
};

std::pair<DenseMatrix<double>, DenseMatrix<double>> SEM::gradFlowDebug() {
//DenseMatrix<double> SEM::gradFlow(double dt) {
	// solve rho_n with rho_n-1
	// setup
	DenseMatrix<double> grad_f(mesh->nFaces(), 3); 
	geometry->requireCotanLaplacian();
	geometry->requireVertexLumpedMassMatrix();
	SparseMatrix<double>& L = geometry->cotanLaplacian;
	SparseMatrix<double>& A = geometry->vertexLumpedMassMatrix;
	//double dt = std::min(rho_v.minCoeff() / rho_v.mean(), rho_v.mean() / rho_v.maxCoeff()) * M_PI;
	// compute rho_vn
	// Eigen::SimplicialLDLT<SparseMatrix<double>> solver(A+dt*L);
	// DenseMatrix<double> rho_vn = solver.solve(A * rho_v);
	// compute grad rho_fn
	double dt = 0;
	for (Face f : mesh->faces())
	{
		Vector3 v[3];
		int v_idx[3];
		int n = 0;
		for (auto vi : f.adjacentVertices())
		{
			v_idx[n] = vi.getIndex();
			v[n++] = geometry->vertexPositions[vi];
		}
		// i: 0; j: 1; k:2
		Vector3 vjk = v[2] - v[1];
		Vector3 vki = v[0] - v[2];
		Vector3 vij = v[1] - v[0];
		Vector3 N = cross(vjk, vki);
		double area = N.norm() / 2;
		Vector3 temp = rho_v(v_idx[0]) * vjk + rho_v(v_idx[1]) * vki + rho_v(v_idx[2]) * vij;
		Vector3 grad = 0.5 / area * cross(N.unit(), rho_v(v_idx[0]) * vjk + rho_v(v_idx[1]) * vki + rho_v(v_idx[2]) * vij);
		grad_f(f.getIndex(), 0) = grad.x;
		grad_f(f.getIndex(), 1) = grad.y;
		grad_f(f.getIndex(), 2) = grad.z;
		
	}
	computeWfv(); // with direction correction(?), this line recompute W_fv internally
	DenseMatrix<double> grad_v =  W_fv * grad_f;
	// projection on boundary to be tangential
	Vector3 center{ 2, 0, 0 }; //
	for (auto i : bdry_vec)
	{
		Vector3 n = geometry->vertexPositions[i] - center;// normal vector
		Vector3 grad_cur{ grad_v(i, 0), grad_v(i, 1), grad_v(i, 2) };
		grad_cur -= dot(grad_cur, n) * n;
		grad_v(i, 0) = grad_cur.x;
		grad_v(i, 1) = grad_cur.y;
		grad_v(i, 2) = grad_cur.z;
	}
	for (Vertex v : mesh->vertices())
	{
		int i = v.getIndex();
		Vector3 grad_v_cur = Vector3{ grad_v(i,0), grad_v(i, 1), grad_v(i, 2) };
		geometry->vertexPositions[v] = geometry->vertexPositions[v] + dt * grad_v_cur; // the sign depends on 
		//geometry->vertexPositions[v] = geometry->vertexPositions[v] + grad_v_cur;
	}
	for (auto i : bdry_vec) // projection after advection
	{
		geometry->vertexPositions[i] = (geometry->vertexPositions[i] - center).unit() + center;
	}
	// rho_v = rho_vn;

	
	//center /= bdry_vec.size();a
	//for (auto i: bdry_vec)
	//	geometry->vertexPositions[i] = (geometry->vertexPositions[i]-center).unit;
	return std::make_pair(grad_v, grad_f);
};


double SEM::stepDEMFixBoundary(double dt) {
	double eps = 1e-4;
	//double dt = stepLength(grad_v) * 0.8; // adaptive steplength
	gradFlow(dt);
	delaunayFlip();
	double std_div_mean = sqrt((rho_v.array() - rho_v.mean()).square().sum() / rho_v.size()) / rho_v.mean();
	return std_div_mean;
};

void SEM::debug() {
	// 1. return the boudary value for debug and register is as point clouds
	//std::cout << "debug";
	//std::cout << "v:\n[";
	//for (int i=0; i!=ORIGINAL.size() ; i++)
	//{
	//	std::cout << ORIGINAL[i].x << ", " << ORIGINAL[i].y << ", " << ORIGINAL[i].z << ";\n";
	//}
	//std::cout << "]\n";
	//std::cout << "initial:\n[";
	//for (Vertex v : mesh->vertices())
	//{
	//	std::cout << geometry->vertexPositions[v].x << " " << geometry->vertexPositions[v].y << ";\n";
	//}
	//std::cout << "]\n";
	//std::cout << "f:\n[";
	//for (Face f : mesh->faces())
	//{
	//	for (Vertex v : f.adjacentVertices())
	//	{
	//		std::cout << v.getIndex() + 1 <<" ";
	//	}
	//	std::cout << ";\n";
	//}
	//std::cout << "]\n";
	//std::cout << "population:\n[";
	//for (auto a : area_vec)
	//	std::cout << a << " ";
	//std::cout << "]\n";

	// temporarily changes the vertex value only for debugging 
}


void SEM::computeSEMAll(double eps)
{
	auto begin = std::chrono::steady_clock::now();
	std::cout << "before starting:" << loss() << std::endl;
	double loss = step();
	std::cout << "relative total area diff:" << loss << std::endl;
	int count = 0;
	while (loss > eps && count < 10)
	{
		count++;
		loss = step();
		std::cout << "relative total area diff:" << loss << std::endl;
	}
	if (count == 10)
		std::cout << "reached maximum iteration number!\n";
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-begin;
	std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
}