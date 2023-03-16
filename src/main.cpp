#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/direction_fields.h"
#include "SEM.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "args/args.hxx"
#include "imgui.h"
#include "polyscope/point_cloud.h"
#include "happly.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"


using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh_uptr;
std::unique_ptr<VertexPositionGeometry> geometry_uptr;
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;


// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psMesh;
polyscope::SurfaceMesh *paramMesh;
// Some algorithm parameters
int iter = 10;
double dt = 0.1;
int face_idx = 0;
SEM* sem_ptr;
// Example computation function -- this one computes and registers a scalar
// quantity

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void debug() {
    //sem_ptr->debug();
    //std::vector<std::vector<size_t>> faces = { {0, 1, 2}, {0, 2, 3} };
    //std::vector<Vector3> vertices_pos = { Vector3{0,0,0}, Vector3{2,0,0}, Vector3{1,7,0}, Vector3{-1, 3} };
    //SimplePolygonMesh simpleMesh(faces, vertices_pos);
    //auto lvals = makeManifoldSurfaceMeshAndGeometry(simpleMesh.polygons, simpleMesh.vertexCoordinates);
    //std::unique_ptr<ManifoldSurfaceMesh> mesh_new;
    //std::unique_ptr<VertexPositionGeometry> geometry_new;
    //std::tie(mesh_new, geometry_new) = std::tuple<std::unique_ptr<ManifoldSurfaceMesh>,
    //    std::unique_ptr<VertexPositionGeometry>>(std::move(std::get<0>(lvals)),  // mesh
    //        std::move(std::get<1>(lvals)));
    //delete geometry;
    //delete mesh;
    //mesh = mesh_new.release();
    //geometry = geometry_new.release();
    //geometry->requireCotanLaplacian();
    //std::cout << geometry->cotanLaplacian << std::endl;
    
    // ----- the following part prints the current gradient and scale factor map
    Vector<double> smoothed_rho_v = sem_ptr->computeRhoV();
    //print the scale factor
    Vector<double> rho_v = sem_ptr->rho_v;
    paramMesh->addVertexScalarQuantity("log scale factor", rho_v);
    paramMesh->addVertexScalarQuantity("log scale factor smoothed", smoothed_rho_v);
    auto grad = sem_ptr->gradFlowDebug(); // just for printing
    //std::cout << sem_ptr->W_fv.row(10);
    paramMesh->addVertexVectorQuantity("grad v", grad.first);
    paramMesh->addFaceVectorQuantity("grad f", grad.second);
    
    double m = 0;
    for (int i = 0; i != grad.first.rows(); i++)
    {
        double cur = sqrt(grad.first(i, 0) * grad.first(i, 0) + grad.first(i, 1) * grad.first(i, 1));
        if (cur > m)
            m = cur;
    }
    std::cout << "maximum grad:" << m << "\n";
    // ----- plotting part ends

    // --- this part tests if the step length control works for a simple triangle
    // temporarily changes the vertex value only for debugging 
    /*
    sem_ptr->geometry->vertexPositions[1145] = Vector3{1, 0, 0};
    sem_ptr->geometry->vertexPositions[1143] = Vector3{0, 2, 0};
    sem_ptr->geometry->vertexPositions[1144] = Vector3{-2, 0, 0};
    DenseMatrix<double> grad_temp(3, 2);
    grad_temp << -1, 1, 0, 0.5, -0.4, 0;
    double dt = sem_ptr->stepLength(grad_temp);
    std::cout << dt << "\n";
    std::vector<Vector3> vertices_pos = { Vector3{1,0,0}, Vector3{0,2,0}, Vector3{-2,0,0} };
    std::vector<Vector3> vertices_pos_new = { Vector3{1+dt*grad_temp(0,0),dt*grad_temp(0,1),0}, 
        Vector3{dt*grad_temp(1,0),2+dt*grad_temp(1,1),0}, 
        Vector3{-2+dt*grad_temp(2,0),dt*grad_temp(2,1),0}};
    for (int i = 0; i != 3; i++)
    {
        std::cout << vertices_pos_new[i] << "\n";
    }
    std::vector<std::vector<size_t>> faces = { {0, 1, 2} };
    polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath("old pos"), vertices_pos, faces);
    polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath("new pos"), vertices_pos_new, faces);
    */
}

Vector3 debug1(int face_idx) {
    Face f = mesh->face(face_idx);
    Vector3 fcenter{ 0, 0, 0};
    for (Vertex v : f.adjacentVertices())
        fcenter += sem_ptr->geometry->vertexPositions[v];
    fcenter /= 3;
    //geometry->vertexPositions[10] = fcenter;
    return fcenter;
}

void myCallback() {
    static int countDEM = 0;
    if (ImGui::Button("SEMstep")) {
        static int count = 0;
        sem_ptr->step();
        std::cout << "iteration:" << count++ << "\ttotal error:" << sem_ptr->loss() << "\ttotal error(dual):" << sem_ptr->lossDual()<< std::endl;
        paramMesh->updateVertexPositions(sem_ptr->geometry->vertexPositions);
        polyscope::requestRedraw();
    }
    if (ImGui::Button("DEMstep")) {
        sem_ptr->computeRhoV();
        double error = sem_ptr->stepDEMFixBoundary(dt);
        std::cout << "iteration:" << countDEM++ << "\ttotal error:" << sem_ptr->loss() << "\ttotal error(dual):" << sem_ptr->lossDual()<< std::endl;
        std::vector<Vector3> grad_visual;
        paramMesh->updateVertexPositions(sem_ptr->geometry->vertexPositions);
        polyscope::requestRedraw();
    }
    if (ImGui::Button("Debug")) {
        debug();
        //sem_ptr->debug();
    }
    if (ImGui::Button("Debug1")) {
        polyscope::registerPointCloud("the face", std::vector<Vector3>{debug1(face_idx)});
        //paramMesh->updateVertexPositions(sem_ptr->geometry->vertexPositions);
        polyscope::requestRedraw();
    }
    if (ImGui::Button("AdaptiveDEM")) {
        for (int i = 0; i != iter; i++) {
            sem_ptr->computeRhoV();
            double dt = 0.5 * sem_ptr->stepLength(sem_ptr->gradFlow(0));
            double error = sem_ptr->stepDEMFixBoundary(dt); // TODO: this line is very stupid but for debug I don't care
            std::cout << "step length:" << dt << "\n"; 
            std::cout << "iteration:" << countDEM++ << "\ttotal error:" << sem_ptr->loss() << "\ttotal error(dual):" << sem_ptr->lossDual() << std::endl;
            std::vector<Vector3> grad_visual;
            paramMesh->updateVertexPositions(sem_ptr->geometry->vertexPositions);
            polyscope::requestRedraw();
        }
    }
    if (ImGui::Button("Reset")) {
        sem_ptr->geometry->vertexPositions = sem_ptr->INITIAL;
        paramMesh->updateVertexPositions(sem_ptr->geometry->vertexPositions);
        polyscope::requestRedraw();
    }
    ImGui::SliderInt("iterations", &iter, 1, 1000);
    ImGui::InputDouble("dt", &dt);
    ImGui::InputInt("face to show", &face_idx);
}


void writeOBJMapping(std::string name, SEM* sp) {
    std::ofstream myfile;
    myfile.open (name);
    
    for (Vertex v : mesh->vertices())
        myfile << "v " << sp->ORIGINAL[v.getIndex()].x << " " << sp->ORIGINAL[v.getIndex()].y << " " << sp->ORIGINAL[v.getIndex()].z << "\n";
        //<< sp->geometry->vertexDualArea(v) << " " << sp->geometry->vertexDualArea(v) << " " << sp->geometry->vertexDualArea(v) << "\n";
    for (Vertex v : mesh->vertices())
    {
        myfile << "vn " << sp->geometry->vertexNormals[v].x << " " << sp->geometry->vertexNormals[v].y << " " << sp->geometry->vertexNormals[v].z << "\n";
    }
    double max=0;
    Vector3 center{ 0 };
    for (Vertex v : mesh->vertices())
        center += sp->geometry->vertexPositions[v];
    center /= mesh->nVertices();
    for (Vertex v : mesh->vertices())
        max = std::max(max, (sp->geometry->vertexPositions[v] - center).norm());
    max = max * 2;
    for (Vertex v : mesh->vertices())
        myfile << "vt " << (sp->geometry->vertexPositions[v].x - center.x) / max << " " << (sp->geometry->vertexPositions[v].y - center.y) / max << "\n";
    for (Face f : mesh->faces())
    {
        myfile << "f ";
        for (Vertex v : f.adjacentVertices()) {
            myfile << v.getIndex()+1 << "/"<< v.getIndex()+1 << "/"<< v.getIndex()+1 << " ";
        }
        myfile << "\n";
    }
    myfile.close();
}

void readFiles(std::string filename, std::vector<std::array<double, 3>>& vPos, std::vector<std::array<double, 3>>& pPos, std::vector<std::vector<size_t>> & fInd)
{
    happly::PLYData plyIn(filename);

    vPos = plyIn.getVertexPositions();
    fInd = plyIn.getFaceIndices<size_t>();

    std::vector<double> px = plyIn.getElement("vertex").getProperty<double>("px");
    std::vector<double> py = plyIn.getElement("vertex").getProperty<double>("py");
    std::vector<double> pz = plyIn.getElement("vertex").getProperty<double>("pz");
    pPos.resize(vPos.size());
    for (int i = 0; i != pPos.size(); i++)
    {
        pPos[i][0] = px[i];
        pPos[i][1] = py[i]; // translate to view easier
        pPos[i][2] = pz[i];
    }
}

void writeFiles(std::string filenamein, std::string filenameout, Eigen::MatrixXd & pPos)
{
    happly::PLYData plyIn(filenamein);
    happly::PLYData plyOut;
    plyOut.addFaceIndices(plyIn.getFaceIndices<size_t>());
    plyOut.addVertexPositions(plyIn.getVertexPositions());
    plyOut.getElement("vertex").addProperty<unsigned char>("red", plyIn.getElement("vertex").getProperty<unsigned char>("red"));
    plyOut.getElement("vertex").addProperty<unsigned char>("green", plyIn.getElement("vertex").getProperty<unsigned char>("green"));
    plyOut.getElement("vertex").addProperty<unsigned char>("blue", plyIn.getElement("vertex").getProperty<unsigned char>("blue"));
    std::vector<double> px, py, pz;
    for (int i = 0; i != pPos.rows(); i++) {
        px.push_back(pPos(i, 0));
        py.push_back(pPos(i, 1));
        pz.push_back(pPos(i, 2));
    }
    plyOut.getElement("vertex").addProperty<double>("px", px);
    plyOut.getElement("vertex").addProperty<double>("py", py);
    plyOut.getElement("vertex").addProperty<double>("pz", pz);
    plyOut.write(filenameout);
}

int main(int argc, char **argv) {

    // Configure the argument parser
    //args::ArgumentParser parser("geometry-central & Polyscope example project");
    //args::Positional<std::string> inputFilename(parser, "mesh", "a mesh file.");
    //
    //// Parse args
    //try {
    //  parser.ParseCLI(argc, argv);
    //} catch (args::Help &h) {
    //  std::cout << parser;
    //  return 0;
    //} catch (args::ParseError &e) {
    //  std::cerr << e.what() << std::endl;
    //  std::cerr << parser;
    //  return 1;
    //}

    //// Make sure a mesh name was given
    //
    //std::string name;
    //if (!inputFilename) {
    //  name = "../input/cowhead.obj";
    //  //std::cerr << "Please specify a mesh file as argument" << std::endl;
    //  //return EXIT_FAILURE;
    //}
    //else
    //{
    //    name = args::get(inputFilename);
    //}
  
    std::string filein, fileout;
    double eps = 0.01;
    if (argc > 1)
        filein = argv[1];
    if (argc > 2)
        fileout = argv[2];
    if (argc > 3)
    {
        std::string temp = argv[3];
        eps = stod(temp);
    }



    // --- Load mesh, with our ply format
    std::vector<std::array<double, 3>> vPos, pPos;
    std::vector<std::vector<size_t>> fInd;
    readFiles(filein, vPos, pPos, fInd);
    const std::vector<std::vector<size_t>> polygons = fInd;
    ManifoldSurfaceMesh mesh_ins(polygons);
    VertexPositionGeometry geometry_ins(mesh_ins);
    for (int i = 0; i != vPos.size(); i++) // the most stupid way to read file
    {
        geometry_ins.inputVertexPositions[i] = Vector3{ vPos[i][0], vPos[i][1], vPos[i][2] }; 
    }
    mesh = &mesh_ins;
    geometry = &geometry_ins;
     
    // --- Load mesh, with plain format
    //std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filein);
    //mesh = mesh_uptr.release();
    //geometry = geometry_uptr.release();


    SEM sem(mesh, geometry, CIRCLE);
    sem_ptr = &sem;
    // sem_ptr->computeSEMAll(eps); // this line for generating SEM result

    // Initialize polyscope
    polyscope::init();

    // Set the callback function
    polyscope::state::userCallback = myCallback;
    // Register the mesh with polyscope
    psMesh = polyscope::registerSurfaceMesh("original 3D", sem.ORIGINAL, mesh->getFaceVertexList(), polyscopePermutations(*mesh));    //geometry->inputVertexPositions;
#ifdef USE_DELAUNAY
    paramMesh = polyscope::registerSurfaceMesh("SEM parameterization", geometry->vertexPositions, sem_ptr->intTri->intrinsicMesh->getFaceVertexList());    //geometry->inputVertexPositions;
#else
    paramMesh = polyscope::registerSurfaceMesh("SEM parameterization", geometry->vertexPositions, mesh->getFaceVertexList(), polyscopePermutations(*mesh));    //geometry->inputVertexPositions;
#endif
    polyscope::requestRedraw();

    // --- testing part
    //std::cout << geometry->cotanLaplacian;
    //auto L = sem.computeLaplacian();
    ////std::cout << L << std::endl;
    //auto fB = sem.solveBoundaryLaplacian(sem.fI_save, L);
    //sem.updateGeometry(fB, sem.fI_save);
    //polyscope::registerSurfaceMesh(
    //    polyscope::guessNiceNameFromPath("SEM parameterization"),
    //    geometry->vertexPositions, mesh->getFaceVertexList(),    //geometry->inputVertexPositions
    //    polyscopePermutations(*mesh));

    // --- output
    //Eigen::MatrixXd param(vPos.size(), 3);
    //for (int i = 0; i != vPos.size(); i++)
    //{
    //    param(i, 0) = geometry_ins.inputVertexPositions[i][0];
    //    param(i, 1) = geometry_ins.inputVertexPositions[i][1];
    //    param(i, 2) = geometry_ins.inputVertexPositions[i][2];
    //}
    //writeFiles(filein, fileout, param);
    // ---

    // test 
    //writeOBJMapping("alex_mapped.obj", sem_ptr);
    polyscope::show();
    //delete mesh;
    //delete geometry;
    return EXIT_SUCCESS;
}
