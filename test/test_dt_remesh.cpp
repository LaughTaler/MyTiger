#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
// #include "discretize_surface.h"
// #include "discretize_surface_parameters.h"
// #include "discretize_curve.h"
// #include "discretize_curve_parameters.h"
// #include "geometry_data_struct.h"
// #include "geometry_interface_occ.h"
// #include "geometry_interface_parameters.h"
// #include "sizing_field_interface.h"
// #include "../tools/meshIO.h"
// #include "../tools/output_info.hpp"
// #include "../tools/test_linemesh.hpp"
// #include "catch.hpp"
// #include "matrix_to_list.hpp"
// #include "list_to_matrix.hpp"
// #include "merge_surface_mesh.hpp"
// #include <chrono>
// #include <thread>
#include "test_api.h"
#include "alias.h"

// TEST_CASE("tiger 1.5 base", "test_remesh") {

// 	TiGER::Geometry Geo;
//   TiGER::readModelParameters args2;
// 	TiGER::Read_Model("../cylinder.stp",Geo,args2);

//         {
//           CurveParametersDimension dimension_args;
//           // dimension_args.setangle(20);
//           dimension_args.setdimension(3);
//           //dimension_args.setdelta_s(0.1);
//           //dimension_args.setmax_dimension(1024);
//           std::vector<LineMesh> segements;
//           TiGER::discretize_curve::Grid1D_Dimension(
//               Geo.curves_, dimension_args, segements);
//           std::vector<std::vector<LineMesh>> sf_boundary;
//           sf_boundary.resize(Geo.surfaces_.size());
//           for (int i = 0; i < Geo.topo_surf_to_curves_.size(); i++) {
//             for (int j = 0; j < Geo.topo_surf_to_curves_[i].size(); j++) {
//               sf_boundary[i].push_back(
//                   segements[Geo.topo_surf_to_curves_[i][j]]);
//             }
//           }

//           IsotropicSurfaceParametersTri surface_args;
//           surface_args.setmUserSpecifiedMaxEdgeLen(1);
//           surface_args.setmUserSpecifiedMinEdgeLen(1);
//           SizingFunction size_field2 = [](const double& x, const double& y,
//                                           const double& z) { return 0.1; };
//           SurfaceMesh sf_mesh2;
//           std::vector<SurfaceMesh> mesh_(Geo.surfaces_.size());
//           //for( int i = 3; i <4; i++ )
//           int i = 3;
//            TiGER::discretize_surface::CADSurf_exportToSTL(Geo.surfaces_[i],
//            sf_boundary[i], surface_args, size_field2,sf_mesh2);
          
//           auto start = std::chrono::high_resolution_clock::now();

//           TiGER::merge_surface_mesh(mesh_, sf_mesh2);

//           auto end = std::chrono::high_resolution_clock::now();
//           std::chrono::duration<double> elapsed = end - start;

//           std::cout << "Elapsed time: " << elapsed.count() << "s\n";
        
//           TiGER::Mesh mesh;
//           TiGER::list_to_matrix<3>(sf_mesh2.coords, mesh.Vertex);
//           TiGER::list_to_matrix<3>(sf_mesh2.tris, mesh.Topo);
//           // mesh.vertex
//           TiGER::MESHIO::writeVTK("../test_dt_remesh_output.vtk", mesh);
//           // REQUIRE(surf_out.coords.size() > 0);
//         }
	
// }

TEST_CASE("surfacemesh test", "[test]")
{
    TiGER::Geometry Geo;
    test_load_geometry("f6.step", Geo);

    REQUIRE(get_geo_size(Geo) != 0);
    // if( get_geo_size(Geo) == 0 )
    // {
    //     return -1;
    // }
    std::cout << "after testing load geo" << std::endl;
    std::cout << "面网格测试\n";
    // 先保存所有
    std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> all_curves_surfaces_(Geo.curves_.size());
    for( int i = 0; i < Geo.topo_surf_to_curves_.size(); i++ )
    {
        for( int j = 0; j < Geo.topo_surf_to_curves_[i].size(); j++ )
        {
            all_curves_surfaces_[Geo.topo_surf_to_curves_[i][j]].push_back(Geo.surfaces_[i]);
        }
    }

    SECTION("surfacemesh test 1"){
        for( int i = 0; i < 53; i++ )
        {
            std::cout << "face i:" << i << std::endl;
            CurveParametersDimension dimension_args;
            // 先生成一个粗的线网格用于获取包围盒对角线长度
            dimension_args.setdimension(10);
            std::vector<LineMesh> segements;
            test_line_mesh_gen(Geo, dimension_args, segements);
            double distance = getDiagonalLength(segements);
            std::cout << "distance" << distance << std::endl;
            for( double mMaxDeviation1 = 0.001; mMaxDeviation1 <= 0.01; mMaxDeviation1 += 0.001 )
            {
                double mUserSpecifiedMaxEdgeLen1 = 0.5 * distance;
                double mUserSpecifiedMinEdgeLen1 = 0.1 * distance;
                std::vector<LineMesh> segements1;
                CurveParametersDimension dimension_args1;
                dimension_args1.setdimension(3);
                dimension_args1.setdelta_s(mUserSpecifiedMaxEdgeLen1);
                dimension_args1.setdeviation(mMaxDeviation1);
                dimension_args1.setuse_database_curvature(true);
                ////////////////////////////////////////////////////////
                std::vector<std::shared_ptr<GeometryCurve>> curves;
                // std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_(
                //     Geo.topo_surf_to_curves_[i].size());
                std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_;
                for( auto curve_id : Geo.topo_surf_to_curves_[i] )
                {
                    curves.push_back(Geo.curves_[curve_id]);
                    surfaces_.push_back(all_curves_surfaces_[curve_id]);
                }
                std::cout << "before line mesh\n";
                test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
                REQUIRE(get_linemesh_size(segements1) != 0);
                // if( get_linemesh_size(segements1) == 0 )
                // {
                //     return -1;
                // }
                std::cout << "after line mesh\n";
                IsotropicSurfaceParametersTri surface_args1;
                surface_args1.setmUserSpecifiedMaxEdgeLen(mUserSpecifiedMaxEdgeLen1);
                surface_args1.setmUserSpecifiedMinEdgeLen(mUserSpecifiedMinEdgeLen1);
                surface_args1.setmMaxDeviation(mMaxDeviation1);
                surface_args1.setmSwapCellWithNoInteriorPoints(false);
                SurfaceMesh surfacemesh;
                test_surf_mesh_gen_dt(Geo.surfaces_[i], surface_args1, segements1, surfacemesh);
                REQUIRE(get_surfacemesh_size(surfacemesh) != 0);
                // if( get_surfacemesh_size(surfacemesh) == 0 )
                // {
                //     // 返回值
                //     return -1;
                // }
            }
        }
    }

    SECTION("surfacemesh test 2")
    {
        for( int i = 0; i <= 53; i++ )
        {
            std::cout << "face i:" << i << std::endl;
            std::vector<LineMesh> segements1;
            CurveParametersDimension dimension_args1;
            // std::cout << dimension_args1.getangle() << std::endl;
            dimension_args1.setdimension(3);
            dimension_args1.setdelta_s(100);
            // std::cout << mMaxDeviation1;
            dimension_args1.setdeviation(0.002);
            // // dimension_args1.setangle(0);
            dimension_args1.setuse_database_curvature(true);
            ////////////////////////////////////////////////////////
            std::vector<std::shared_ptr<GeometryCurve>> curves;
            // std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_(
            //     Geo.topo_surf_to_curves_[i].size());
            std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_;
            for( auto curve_id : Geo.topo_surf_to_curves_[i] )
            {
                curves.push_back(Geo.curves_[curve_id]);
                surfaces_.push_back(all_curves_surfaces_[curve_id]);
            }
            std::cout << "before line mesh\n";
            test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
            getDiagonalLength(segements1);
            double mUserSpecifiedMaxEdgeLen1 = 100;
            double mUserSpecifiedMinEdgeLen1 = 0.2;
            double mMaxDeviation1 = 0.002;
            IsotropicSurfaceParametersTri surface_args1;
            surface_args1.setmUserSpecifiedMaxEdgeLen(mUserSpecifiedMaxEdgeLen1);
            surface_args1.setmUserSpecifiedMinEdgeLen(mUserSpecifiedMinEdgeLen1);
            surface_args1.setmMaxDeviation(mMaxDeviation1);
            surface_args1.setmSwapCellWithNoInteriorPoints(false);
            SurfaceMesh surfacemesh;
            test_surf_mesh_gen_dt(Geo.surfaces_[i], surface_args1, segements1, surfacemesh);
            checkMeshQuality(surfacemesh);
        }
    }
}