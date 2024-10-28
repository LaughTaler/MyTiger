#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "matrix_to_list.hpp"
#include "meshIO.h"
#include "alias.h"
#include "discretize_surface_remesh.h"
#include "discretize_surface_remesh_parameters.h"



TEST_CASE("tiger 1.5 base", "test_remesh") {
  TiGER::Mesh m;
  //TiGER::MESHIO::readVTK("/home/ltr/workspace/ti-ger-1.5/build/debugout.vtk", m,
  //                       "surface_id");
  TiGER::MESHIO::readVTK("C:/workspace/tiger1.5/build/ship_in.vtk", m);
//   TiGER::MESHIO::readVTK("/home/ltr/workspace/ti-ger-1.5/file/test_dt.vtk", m,
//                          "surface_id");
  TiGER::SurfaceMesh surf;
  TiGER::matrix_to_list<3>(m.Vertex, surf.coords);
  TiGER::matrix_to_list<3>(m.Topo, surf.tris);
  // TiGER::matrix_to_list<1>(m.Masks, surf.attribute_int);
  for (int i = 0; i < surf.tris.size(); i++) surf.attribute_int.push_back(0);

  TiGER::RemeshParameters args;
  args.setiteration_number(5);


  TiGER::SurfaceMesh surf_out;
  std::vector<std::vector<int>> constrained_edge;
  std::vector<std::vector<int>> conformaing_edge;
  std::function<double(double, double, double)> size_function;
  std::vector<std::vector<int>> constrained_edge_out;

  TiGER::discretize_surface_remesh::TMeshSurf_Opt(
      surf, constrained_edge, conformaing_edge, args, size_function, surf_out,
      constrained_edge_out);
  TiGER::Mesh m_o;
  m_o.Vertex.resize(surf_out.coords.size(), 3);
  m_o.Topo.resize(surf_out.tris.size(), 3);
  for (int i = 0; i < surf_out.coords.size(); ++i) {
    for (int j = 0; j < 3; ++j)
  m_o.Vertex(i, j) = surf_out.coords[i][j];
  }
  for (int i = 0; i < surf_out.tris.size(); ++i) {
    for (int j = 0; j < 3; ++j) m_o.Topo(i, j) = surf_out.tris[i][j];
  }
  TiGER::MESHIO::writeVTK("C:/workspace/tiger1.5/build/ship_out.vtk", m_o);
  REQUIRE(surf_out.coords.size() > 0);
}