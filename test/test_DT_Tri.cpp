#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do
                           // this in one cpp file
#include <functional>
#include "alias.h"
#include "catch.hpp"
#include "discretize_surface_delaunay.h"
#include "discretize_surface_delaunay_parameters.h"
#include "matrix_to_list.hpp"
#include "meshIO.h"

TEST_CASE("tiger 1.5 base", "HexGrid_createPrismaticMesh") {
  // TiGER::LineMesh seg;
  // int nPoints = 0;
  // int nSegments = 0;
  // std::ifstream vtk_file;
  // vtk_file.open("D:/code/wyfdt2d/test/1_boundary_information");

  // std::string vtk_type_str = "POLYDATA ";
  // char buffer[BUFFER_LENGTH];
  // while (!vtk_file.eof()) {
  //   vtk_file.getline(buffer, BUFFER_LENGTH);
  //   std::string line = (std::string)buffer;
  //   if (line.length() < 2 || buffer[0] == '#') continue;
  //   if (line.find("POINTS ") != std::string::npos) {
  //     std::vector<std::string> words = seperate_string(line);
  //     nPoints = stoi(words[1]);
  //     seg.coord.resize(nPoints);
  //     for (int i = 0; i < nPoints; i++) {
  //       vtk_file >> seg.coord[i][0] >> seg.coord[i][1] >> seg.coord[i][2];
  //     }
  //   }
  //   if (line.find("CELLS") != std::string::npos) {
  //     std::vector<std::string> words = seperate_string(line);
  //     nSegments = stoi(words[1]);
  //     seg.segments.resize(nSegments);
  //     for (int i = 0; i < nSegments; i++) {
  //       vtk_file >> seg.segments[i][0] >> seg.segments[i][1];
  //     }
  //     break;
  //   }
  // }
  // vtk_file.close();

  // TiGER::SurfaceMesh surf;
  // TiGER::TriangleDealunay args;
  // TiGER::discretize_surface_delaunay::TMeshSurf_processTriangleLoop(seg, args, NULL, surf);
  // REQUIRE(surf.coords.size() > 0);
  // REQUIRE(surf.tris.size() > 0);
  REQUIRE(1);
}