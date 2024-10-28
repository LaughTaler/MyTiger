#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "boundary_layer_mesh.h"
#include "geometry_data_struct.h"

#include "tetrahedral_mesh_parameters.h"
#include "tetrahedral_mesh.h"

#include "matrix_to_list.hpp"
#include "meshIO.h"
#include "alias.h"

TEST_CASE("tiger 1.5 base", "HexGrid_createPrismaticMesh") {

	TiGER::Mesh m;
	TiGER::MESHIO::readVTK("cubic.o.vtk", m,"surface_id");
	TiGER::SurfaceMesh surf;
	TiGER::matrix_to_list<3>(m.Vertex, surf.coords);
	TiGER::matrix_to_list<3>(m.Topo, surf.tris);
	TiGER::matrix_to_list<1>(m.Masks, surf.attribute_int);
	for (auto &i : surf.attribute_int)
		i++;
	TiGER::BoundaryLayerParameters args;
	args.setstep_length(0.01);
	args.setnumber_of_layer(10);
	args.setratio(1.2);
	args.setmultiple_normal(false);
	for (int i = 1; i < 2; i++)
		surf.bc[i] = TiGER::BoundaryCondition::WALL;
	for (int i = 0; i < 1; i++)
		surf.bc[i] = TiGER::BoundaryCondition::WALL; 
	TiGER::VolumeMesh vol;
	TiGER::SurfaceMesh nsurf;
	TiGER::SurfaceMesh nouter;
	std::vector <std::shared_ptr<TiGER::GeometrySurface>> gt;
        std::vector<int> l2g;
	TiGER::boundary_layer::HexGrid_createPrismaticMesh(surf, gt, args, vol, nsurf, nouter,l2g);
	REQUIRE(nsurf.coords.size() > 0);
	REQUIRE(vol.coords.size() > 0);
}


TEST_CASE("tiger 1.5 base multiple", "HexGrid_createPrismaticMesh") {
        TiGER::Mesh m;
        TiGER::MESHIO::readVTK("naca2412.o.vtk", m, "surface_id");
        TiGER::SurfaceMesh surf;
        TiGER::matrix_to_list<3>(m.Vertex, surf.coords);
        TiGER::matrix_to_list<3>(m.Topo, surf.tris);
        TiGER::matrix_to_list<1>(m.Masks, surf.attribute_int);
        for (auto &i : surf.attribute_int) i++;
        TiGER::HybridParameters args;
        args.setstep_length(0.01);
        args.setnumber_of_layer(20);
        args.setratio(1.2);
        args.setmultiple_normal(true);
        args.setaniso_iso_blend(0);
        args.setoutputio(true);
        args.setoptloop(1);
        args.setconstrain(1);
        args.setlayer({2});
        args.setnthread(1);

        for (int i = 0; i < 111; i++) surf.bc[i] = TiGER::BoundaryCondition::WALL;
        TiGER::VolumeMesh vol;
        TiGER::SurfaceMesh nsurf;
        TiGER::SurfaceMesh nouter;
        std::vector<std::shared_ptr<TiGER::GeometrySurface>> gt;
        TiGER::boundary_layer::HexGrid_createHybridMesh(surf, gt, args, vol,
                                                  nsurf);
        TiGER::MESHIO::writeVTK("f6_vol.vtk",vol);
        REQUIRE(nsurf.coords.size() > 0);
        REQUIRE(vol.coords.size() > 0);
}
