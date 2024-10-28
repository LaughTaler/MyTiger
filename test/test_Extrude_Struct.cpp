#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "extrude.h"
#include "extrude_parameters.h"
#include "matrix_to_list.hpp"
#include "../tools/meshIO.h"
#include "alias.h"
#include <igl/readSTL.h>
#include <igl/per_vertex_normals.h>

TEST_CASE("tiger 1.5 base", "test_extrude") {
	TiGER::Mesh m;
	TiGER::MESHIO::readVTK("./file/cubic/cubic.o.vtk", m);
	TiGER::SurfaceMesh surf;
	TiGER::matrix_to_list<3>(m.Vertex, surf.coords);
	TiGER::matrix_to_list<3>(m.Topo, surf.tris);
	TiGER::matrix_to_list<1>(m.Masks, surf.attribute_int);

	TiGER::VolumeMesh vol_mesh;
	Eigen::MatrixXd N;
	TiGER::ExtrudeParameters args;
	TiGER::ExtrudeParameters::GeometricProgressionPara gpp;
	gpp.initial_delta_s = 0.01;
	gpp.growth_rate = 1.1;
	args.setgeo_progre_para(gpp);
	args.setnum_layer(40);
	igl::per_vertex_normals(m.Vertex, m.Topo, N);
	//printf("%d %d %d \n", V.rows(),F.rows(), N.rows());
	//TiGER::extrude::readstltoSurfaceMesh(surface_mesh, V, F, N);
	TiGER::matrix_to_list<3>(N, surf.point_normal);
	TiGER::extrude::HexGrid_Extrude3D(surf, args, vol_mesh);
	REQUIRE(vol_mesh.coords.size()>0);
	REQUIRE(vol_mesh.prisms.size() > 0);
	//MESHIO::writeVTK();
}
