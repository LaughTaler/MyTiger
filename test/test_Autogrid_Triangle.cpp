#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "matrix_to_list.hpp"
#include "list_to_matrix.hpp"
#include "meshIO.h"
#include "mesh_data_struct.h"

#include "autogrid.h"
#include "autogrid_parameters.h"
#include "alias.h"



TEST_CASE("tiger 1.5 base", "test_remesh") {
  TiGER::Mesh m, m_out;
  // TiGER::MESHIO::readVTK("/home/ltr/workspace/ti-ger-1.5/build/boxhole.vtk", m,
  //                       "surface_id");
  TiGER::MESHIO::readVTK("C:/workspace/tiger1.5/build/ship_in.vtk", m,
                         "surface_id");

  TiGER::SurfaceMesh surf, surf_out;
  //surf.edges = {{0, 1}};
  TiGER::matrix_to_list<3>(m.Vertex, surf.coords);
  TiGER::matrix_to_list<3>(m.Topo, surf.tris);
  // TiGER::matrix_to_list<1>(m.Masks, surf.attribute_int);

  TiGER::AutogridParameters args;

  args.seteps_rel(1e-3);

  TiGER::discrete_geometry_repair::TMeshSurf_Repair(
      surf,  args, surf_out);

  //std::cout << "constrained edge is :\n";
  //for(int i = 0; i < surf_out.edges.size(); ++i){
  //  std::cout << surf_out.edges[i][0] << " " << surf_out.edges[i][1] << "\n"; 
  //}

  TiGER::list_to_matrix<3>(surf_out.coords, m_out.Vertex);
  TiGER::list_to_matrix<3>(surf_out.tris, m_out.Topo);
  //TiGER::list_to_matrix<3>(surf_out.edges, m_out);

  TiGER::MESHIO::writeVTK("C:/workspace/tiger1.5/build/ship_out.vtk", m_out, "surface_id");

  REQUIRE(surf_out.coords.size() > 0);
}
