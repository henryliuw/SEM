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

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh_uptr;
std::unique_ptr<VertexPositionGeometry> geometry_uptr;
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psMesh;

// Some algorithm parameters
float param1 = 42.0;
SEM* sem_ptr;
// Example computation function -- this one computes and registers a scalar
// quantity
void GUIstep() {
   static int count = 0;
   std::cout << "iteration:" << count++ << "\taverage error:" << sem_ptr->step() << std::endl;
   polyscope::registerSurfaceMesh(
           polyscope::guessNiceNameFromPath("SEM parameterization"),
           geometry->vertexPositions, mesh->getFaceVertexList(),    //geometry->inputVertexPositions
           polyscopePermutations(*mesh));
   polyscope::requestRedraw();
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {

  if (ImGui::Button("step")) {
    GUIstep();
  }

  ImGui::SliderFloat("param", &param1, 0., 100.);
}

void debug() {
    
}

int main(int argc, char **argv) {

  // Configure the argument parser
  args::ArgumentParser parser("geometry-central & Polyscope example project");
  args::Positional<std::string> inputFilename(parser, "mesh", "a mesh file.");
  
  // Parse args
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help &h) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  // Make sure a mesh name was given
  /*
  if (!inputFilename) {
    std::cerr << "Please specify a mesh file as argument" << std::endl;
    return EXIT_FAILURE;
  }
  */

  // Initialize polyscope
  polyscope::init();
  
  // Set the callback function
  polyscope::state::userCallback = myCallback;

  // Load mesh
  std::string name = "cowhead.obj";
  std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh("../input/" + name);
  mesh = mesh_uptr.release();
  geometry = geometry_uptr.release();

  SEM sem(mesh, geometry);
  sem_ptr = &sem;
  // Register the mesh with polyscope
  psMesh = polyscope::registerSurfaceMesh(
      polyscope::guessNiceNameFromPath("original 3D"),
      sem.ORIGINAL, mesh->getFaceVertexList(),    //geometry->inputVertexPositions
      polyscopePermutations(*mesh));

  // testing part
  //std::cout << geometry->cotanLaplacian;
  //auto L = sem.computeLaplacian();
  ////std::cout << L << std::endl;
  //auto fB = sem.solveBoundaryLaplacian(sem.fI_save, L);
  //sem.updateGeometry(fB, sem.fI_save);
  //polyscope::registerSurfaceMesh(
  //    polyscope::guessNiceNameFromPath("SEM parameterization"),
  //    geometry->vertexPositions, mesh->getFaceVertexList(),    //geometry->inputVertexPositions
  //    polyscopePermutations(*mesh));
  // Give control to the polyscope gui
  polyscope::show();
  delete mesh;
  delete geometry;
  return EXIT_SUCCESS;
}
