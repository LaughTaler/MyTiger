#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "mesh_quality.h"
#include "verdict.h"
#include <iostream>
#include <array>
#include "matrix_to_list.hpp"
#include <igl/readSTL.h>
#include "meshIO.h"
#include "alias.h"

TEST_CASE("tiger 1.5 base", "SurfaceMeshTest") {
    //int pnum = 3;
    //for (int i = 0; i < pnum; i++) {
    //    //double p[3];
    //    //g2d->GetPoint(i, p); 
    //    //surface_mesh.coord.push_back({ p[0], p[1], p[2] });
    //    surface_mesh.coords.push_back({ {0.0, 0.0, 0.0} });
    //    surface_mesh.coords.push_back({ {1.0, 0.0, 0.0} });
    //    surface_mesh.coords.push_back({ {0.0, 1.0, 0.0} });
    //}
    //surface_mesh.tris.push_back({ {0, 1, 2} });
    TiGER::SurfaceMesh surf;
    TiGER::Mesh m;
    TiGER::MESHIO::readVTK("./file/cubic/cubic.o.vtk", m);
    TiGER::matrix_to_list<3>(m.Vertex, surf.coords);
    TiGER::matrix_to_list<3>(m.Topo, surf.tris);
    TiGER::matrix_to_list<1>(m.Masks, surf.attribute_int);
    const TiGER::SurfMeshQualityType type = TiGER::SurfMeshQualityType::TRI_ANGLE;
    std::vector<TiGER::QualityResult> quality;
    TiGER::Tool_SurfQuality(surf, type, quality);

    printf("quality size: %d\n", quality.size());
    for (auto i : quality) {
        std::cout << i.ave_value << std::endl;
        std::cout << i.max_value << std::endl;
        std::cout << i.min_value << std::endl;
    }
    REQUIRE(quality.size()>0);
}

