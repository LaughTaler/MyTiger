#define CATCH_CONFIG_MAIN // This tells Catch to provide a main() - only do this in one cpp file
// #include "discretize_surface.h"
// #include "discretize_surface_parameters.h"
// #include "discretize_curve.h"
// #include "discretize_curve_parameters.h"
// #include "geometry_data_struct.h"
// #include "geometry_interface_occ.h"
// #include "geometry_interface_parameters.h"
// #include "sizing_field_interface.h"
// #include "geometry_interface_plane.h"
// #include "../tools/meshIO.h"
// #include "../tools/output_info.hpp"
// #include "../tools/test_linemesh.hpp"
// #include "catch.hpp"
// #include "matrix_to_list.hpp"
// #include "list_to_matrix.hpp"
// #include "merge_surface_mesh.hpp"
// #include "MeshMetric.h"
// #include <chrono>
// #include <thread>
#include "test_api.h"
#include "alias.h"
//std::vector<double> size_s;
//int size_cnt = -1;
//double size_f(double x, double y, double z){
//    size_cnt++;
//    return size_s[size_cnt];
//};
TEST_CASE("tiger 1.5 base", "test_remesh")
{
   /* std::vector<std::array<double, 3>> cmp_V;
   /* std::vector<std::array<double, 3>> cmp_V;
    std::vector<std::array<int, 3>> cmp_F;
    std::vector<std::array<int, 4>> cmp_T;
    std::vector<int> cmp_F_geo;
    std::vector<int> cmp_T_geo;
    TiGER::MESHIO::readVTK_WYF_HW("d:/tiger2/test/HW1_CMP.VTK", cmp_V, cmp_F, cmp_T, cmp_F_geo, cmp_T_geo);
    TiGER::VolumeMesh cmp_vol;
    cmp_vol.coords.assign(cmp_V.begin(), cmp_V.end());
    cmp_vol.tetras.assign(cmp_T.begin(), cmp_T.end());*/

    // Reading an isotropic size field from the benchmark mesh
    bool input_file = false;
    if( input_file )
    {
        TiGER::Geometry Geo;
        TiGER::readModelParameters readmodel_args;
        TiGER::Read_Model("d:/tiger2/test/f6_cf/f6.step", Geo, readmodel_args);
        std::string filename = "D:/ti-ger-1.5/ti-ger-1.5/build/0_surface_information";
        std::array<double, 3> sizefiled_info;
        std::vector<std::array<double, 3>> surface_info;
        std::vector<LineMesh> boundary_segments;
        TiGER::IsotropicSurfaceParametersTri iso_args;
        TiGER::read_discrete_surface_info(filename, sizefiled_info, surface_info, boundary_segments, iso_args);
        SizingFunction size_field = [](const double& x, const double& y, const double& z) { return 1200; };
        auto surface = Geo.surfaces_[76];
        // surface->SmoothFace();
        /*iso_args.setmMaxAngle(24);
        iso_args.setmUserSpecifiedMaxEdgeLen(475.356);
        iso_args.setmUserSpecifiedMinEdgeLen(47.5356);*/
        std::array<double, 3> mid_pt, mid_du, mid_dv, mid_dudu, mid_dvdv, mid_dudv, left_mid_pt, right_mid_pt;
        std::array<double, 4> uv_scale = surface->getUVScaleFunction();
        std::array<double, 2> mid_uv = {(uv_scale[1] + uv_scale[0]) / 2, (uv_scale[3] + uv_scale[2]) / 2};
        std::array<double, 2> left_mid_uv = {(mid_uv[0] + uv_scale[0]) / 2, (mid_uv[1] + uv_scale[2]) / 2};
        std::array<double, 2> right_mid_uv = {(mid_uv[0] + uv_scale[1]) / 2, (mid_uv[1] + uv_scale[3]) / 2};
        mid_pt = surface->d0Function(mid_uv);
        left_mid_pt = surface->d0Function(left_mid_uv);
        right_mid_pt = surface->d0Function(right_mid_uv);
        surface->d1Function(mid_uv, mid_du, mid_dv);
        surface->d2Function(mid_uv, mid_dudu, mid_dvdv, mid_dudv);
        TiGER::SurfaceMesh sf_mesh;
        TiGER::discretize_surface::CADSurf_TGrid_AFT(surface, boundary_segments, iso_args, nullptr, sf_mesh);
        TiGER::Mesh mesh;
        TiGER::list_to_matrix<3>(sf_mesh.coords, mesh.Vertex);
        TiGER::list_to_matrix<3>(sf_mesh.tris, mesh.Topo);
        // mesh.vertex
        TiGER::MESHIO::writeVTK("test_AFT_Tri_output.vtk", mesh);
    }
    else
    {
        /*std::ifstream f;
        f.open("D:/weixin/WeChat Files/wxid_qprxgsh4njsk22/FileStorage/File/2024-08/sizes.txt");
        if( !f.is_open() )
        {
            std::cout << "No such file. - "
                      << "filename" << std::endl;
            return;
        }
        string buffer;
        while (getline(f, buffer))
        {
            size_s.push_back(std::stod(buffer));
        }*/
   

        /*std::ifstream f;
        f.open("D:/weixin/WeChat Files/wxid_qprxgsh4njsk22/FileStorage/File/2024-08/sizes.txt");
        if( !f.is_open() )
        {
            std::cout << "No such file. - "
                      << "filename" << std::endl;
            return;
        }
        string buffer;
        while (getline(f, buffer))
        {
            size_s.push_back(std::stod(buffer));
        }*/
   

        TiGER::Geometry Geo;
        TiGER::readModelParameters readmodel_args;
        SizingFunction size_field2 = [](const double& x, const double& y, const double& z) { return 0.1; };

        TiGER::Read_Model("D:/tiger2/test/f6.step", Geo, readmodel_args);
        /*TiGER::Read_Model(
            "D:/CAD_TestSet/D_freesurface/D_freesurface.stp", Geo,
            readmodel_args);*/
        CurveParametersDimension dimension_args;
        CurveParametersGeometric geometric_args;
        
        
        //vector<int> last_linemesh_size(Geo.curves_.size(), 0.0);
        //vector<bool> test_linemesh(Geo.curves_.size(), true);
        //for (double test_angle = 0.0001; test_angle <= 0.001; test_angle *= 2)
        //{
        //    std::cout << test_angle << endl;
        //    //dimension_args.setangle(test_angle);
        //    // TiGER::CurveParameters curve_args;
        //    // dimension_args.setangle(10);
        //   dimension_args.setdeviation(test_angle);
        //    dimension_args.setdimension(2);
        //    dimension_args.setuse_database_curvature(true);
        //    std::vector<LineMesh> segements_222(Geo.curves_.size());
        //    std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_222(Geo.curves_.size());
        //    for( int i = 0; i < Geo.topo_surf_to_curves_.size(); i++ )
        //    {
        //        for( int j = 0; j < Geo.topo_surf_to_curves_[i].size(); j++ )
        //        {
        //            surfaces_222[Geo.topo_surf_to_curves_[i][j]].push_back(Geo.surfaces_[i]);
        //        }
        //    }
        //    TiGER::discretize_curve::discreteCurveDimension(Geo.curves_, surfaces_222, dimension_args, segements_222);
        //    if( test_angle == 0.0001 )
        //    {
        //        for( int i = 0; i < Geo.curves_.size(); i++ )
        //        {
        //            last_linemesh_size[i] = segements_222[i].coord.size();
        //        }
        //    }
        //    else
        //    {
        //        for( int i = 0; i < Geo.curves_.size(); i++ )
        //        {
        //            double coord_diff = last_linemesh_size[i] - segements_222[i].coord.size();
        //            if( coord_diff <= 0.0 )
        //                test_linemesh[i] = false;
        //            std::cout << i<<":"<<coord_diff / (double)(last_linemesh_size[i]-1) << " ";
        //        }
        //        std::cout << endl;
        //        for( int i = 0; i < Geo.curves_.size(); i++ )
        //        {
        //            last_linemesh_size[i] = segements_222[i].coord.size();
        //        }

        //    }

        //}
        //for( int i = 0; i < test_linemesh.size(); i++ )
        //    if( test_linemesh[i] )
        //        std::cout << i << endl;
         //dimension_args.setangle(5);
        TiGER::CurveParameters curve_args;
        //dimension_args.setangle(10);
        //dimension_args.setdeviation(0.01);
        dimension_args.setdimension(10);
        //dimension_args.setdelta_s(0.16);
        //dimension_args.setuse_database_curvature(true);
       // dimension_args.setdelta_s(0.022034142410201196);
        //dimension_args.setmax_dimension(1024);


        std::vector<LineMesh> segements(Geo.curves_.size());
        //for( int i = 0; i < Geo.curves_.size(); i++ )
        //{
        //    // std::cout << i << "th line \n";
        //    TiGER::discretize_curve::Grid1D_Legacy(Geo.curves_[i], size_field2, curve_args, segements[i]);
        //}
        std::vector<LineMesh> segements_geo;
        std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_(Geo.curves_.size());
        for( int i = 0; i < Geo.topo_surf_to_curves_.size(); i++ )
        {
            for( int j = 0; j < Geo.topo_surf_to_curves_[i].size(); j++ )
            {
                surfaces_[Geo.topo_surf_to_curves_[i][j]].push_back(Geo.surfaces_[i]);
            }
        }
        TiGER::discretize_curve::Grid1D_Dimension(Geo.curves_, surfaces_, dimension_args, segements);

        std::vector<std::vector<LineMesh>> sf_boundary;
        sf_boundary.resize(Geo.surfaces_.size());
        for( int i = 0; i < Geo.topo_surf_to_curves_.size(); i++ )
        {
            for( int j = 0; j < Geo.topo_surf_to_curves_[i].size(); j++ )
            {
                sf_boundary[i].push_back(segements[Geo.topo_surf_to_curves_[i][j]]);
            }
        }
  
  
        IsotropicSurfaceParametersTri surface_args;
        //surface_args.setmUserSpecifiedMaxEdgeLen(100);
        //surface_args.setmUserSpecifiedMinEdgeLen(0.0001);
       // surface_args.setmMaxAngle(10);
         //surface_args.setmMaxDeviation(0.002);
        /* surface_args.setOutputInformation(true);
        surface_args.setOutputInformation(true);
        surface_args.setSavepath("");*/
       

        SurfaceMesh sf_mesh2;
        std::vector<SurfaceMesh> mesh_(Geo.surfaces_.size());

        auto start = std::chrono::high_resolution_clock::now();

        /*for( int i = 0; i < 1; i++ )
        {
            std::cout << i << "th face \n";
            TiGER::discretize_surface::CADSurf_TGrid_AFT(Geo.surfaces_[i], sf_boundary[i], surface_args, size_field2,
                                                            mesh_[i]);
            if( !mesh_[i].tris.size() )
                std::cout << i << "th face failed\n";
            else
                std::cout << "element num:" << mesh_[i].tris.size() << std::endl;
        }*/
        // surface_args.setmMaxAngle(20);
        /*surface_args.setmUserSpecifiedMaxEdgeLen(6);
        surface_args.setmUserSpecifiedMinEdgeLen(0.6);*/
        // size_field2 = nullptr;

         //for (int i = 31; i < 32; i++) {
         //  std::cout << i << "th face \n";
         //  TiGER::discretize_surface::discretizeSurfaceTri(
         //      Geo.surfaces_[i], sf_boundary[i], surface_args,
         //      nullptr/*sf.sf[0]*/ ,
         //      mesh_[i]);
         //  if (!mesh_[i].tris.size())
         //      std::cout << i << "th face failed\n";
         //  else
         //    std::cout << "element num:" << mesh_[i].tris.size()<<std::endl;
         //}
         AnisotropicSurfaceParametersTri aniso_args;
         aniso_args.setTopology_priority(false);
         aniso_args.setMaxLayer(100);
          for (int i = 12; i <13; i++) {
            std::cout << i << "th face \n";
            sf_boundary[i][0].boundary_condition.bctype = BoundaryCondition::WALL;
            sf_boundary[i][0].boundary_condition.delta = 0.01;
           //sf_boundary[i][2].boundary_condition.bctype = BoundaryCondition::MATCH_NO_PUSH;
            sf_boundary[i][2].boundary_condition.bctype = BoundaryCondition::WALL;
            sf_boundary[i][2].boundary_condition.delta = 0.01;
            sf_boundary[i][1].boundary_condition.bctype = BoundaryCondition::WALL;
            sf_boundary[i][1].boundary_condition.delta = 0.01;
            sf_boundary[i][3].boundary_condition.bctype = BoundaryCondition::WALL;
            sf_boundary[i][3].boundary_condition.delta = 0.01;
           // sf_boundary[i][3].boundary_condition.bctype = BoundaryCondition::MATCH_NO_PUSH;
            std::vector<std::shared_ptr<TiGER::GeometryCurve>> curves;
            for( int j = 0; j < Geo.topo_surf_to_curves_[i].size(); j++ )
            {
                curves.push_back(Geo.curves_[Geo.topo_surf_to_curves_[i][j]]);
            }
            TiGER::discretize_surface::discretizeSurfaceTriAFLR(Geo.surfaces_[i], curves, aniso_args,
                                                                sf_boundary[i],size_field2 /*sf.sf[0]*/, mesh_[i]);
            if (!mesh_[i].tris.size())
                std::cout << i << "th face failed\n";
            else
              std::cout << "element num:" << mesh_[i].tris.size()<<std::endl;
          }

        TiGER::merge_surface_mesh(mesh_, sf_mesh2);

        // TiGER::discretize_surface::VCADSurf_TGrid_AFT(
        //     Geo.surfaces_, sf_boundary, surface_args, /*sf.sf[0]*/size_field2 ,
        //     sf_mesh2);

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "total element: " << sf_mesh2.tris.size() << endl;
        std::cout << "total element: " << sf_mesh2.tris.size() << endl;
        std::cout << "Elapsed time: " << elapsed.count() << "s\n";

        TiGER::Mesh mesh;
        TiGER::list_to_matrix<3>(sf_mesh2.coords, mesh.Vertex);
        TiGER::list_to_matrix<3>(sf_mesh2.tris, mesh.Topo);
        // mesh.vertex
        TiGER::MESHIO::writeVTK("test_AFT_Tri_output.vtk", mesh);
        // REQUIRE(surf_out.coords.size() > 0);
    }
}

TEST_CASE("linemesh test", "[test1]")
{
    TiGER::Geometry Geo;
    test_load_geometry("f6.step", Geo);
    if( get_geo_size(Geo) == 0 )
    {
        return ;
    }
    REQUIRE(get_geo_size(Geo) != 0);

    std::cout << "after testing load geo" << std::endl;

    SECTION("linemesh test 1")
    {
        for( int i = 1; i <= 10; i++ )
        {
            std::vector<LineMesh> segements;
            double ttl = Geo.curves_[i]->getTotalLenthFunction();
            // SAMPLE1(discreteCurveDimensionParams, delta_s)
            auto samples = sample_discreteCurveDimensionParams_CURVEID_delta_s
                <discreteCurveDimensionParams, int, &discreteCurveDimensionParams::CURVEID, double,
                &discreteCurveDimensionParams::delta_s>( i, i, 1, 0.01 * ttl, 0.1 * ttl, 0.01 * ttl, 19);
            // auto samples = sampleDiscreteCurveDimension(i, i, 1, 3, 3, 10, 0.01 * ttl, 0.1 * ttl, 0.01 * ttl, 10.0,
            //                                             10.0, 10.0, 0.01, 0.01, 0.01, false, 19);
            // for( auto s : samples ){
            //     std::cout << "CURVEID: " << s.CURVEID << ", dimension: " << s.dimension << ", delta_s: " << s.delta_s
            //               << ", angle: " << s.angle << ", deviation: " << s.deviation
            //               << ", use_database_curvature: " << std::boolalpha << s.use_database_curvature << std::endl;
            //     s.angle = 10.0;
            // }
            test_line_mesh_gen2(Geo, samples, segements);
            REQUIRE(get_linemesh_size(segements) != 0);
            // if( get_linemesh_size(segements) == 0 )
            // {
            //     return -1;
            // }
        }
    }

    SECTION("linemesh test 2")
    {
        for( int i = 11; i <= 20; i++ )
        {
            std::vector<LineMesh> segements;
            auto samples = sample_discreteCurveDimensionParams_CURVEID_dimension_angle
                <discreteCurveDimensionParams, int, &discreteCurveDimensionParams::CURVEID, int,
                &discreteCurveDimensionParams::dimension, double, &discreteCurveDimensionParams::angle>(
                    i, i, 1, 3, 10, 1, 5, 40, 5, 21, true
            );
            // auto samples =
            //     sampleDiscreteCurveDimension(i, i, 1, 3, 10, 1, 1, 1, 1, 5.0, 40.0, 5.0, 0.01, 0.01, 0.01, true, 21);
            test_line_mesh_gen2(Geo, samples, segements);
            REQUIRE(get_linemesh_size(segements) != 0);
            // if( get_linemesh_size(segements) == 0 )
            // {
            //     return -1;
            // }
        }
    }

    SECTION("linemesh test 3")
    {
        for( int i = 21; i <= 30; i++ )
        {
            std::vector<LineMesh> segements;
            auto samples = sample_discreteCurveDimensionParams_CURVEID_dimension
                <discreteCurveDimensionParams, int, &discreteCurveDimensionParams::CURVEID, int,
                &discreteCurveDimensionParams::dimension>(
                    i, i, 1, 0, 100, 10, 17
            );
            // auto samples =
            //     sampleDiscreteCurveDimension(i, i, 1, 0, 100, 10, 1, 1, 1, 5.0, 5.0, 1.0, 0.01, 0.01, 0.01, false, 17);
            // for( const auto& s : samples )
            // {
            //     std::cout << samples.size();
            //     std::cout << "CURVEID: " << s.CURVEID << ", dimension: " << s.dimension << ", delta_s: " << s.delta_s
            //               << ", angle: " << s.angle << ", deviation: " << s.deviation
            //               << ", use_database_curvature: " << std::boolalpha << s.use_database_curvature << std::endl;
            // }
            test_line_mesh_gen2(Geo, samples, segements);
            REQUIRE(get_linemesh_size(segements) != 0);
            // if( get_linemesh_size(segements) == 0 )
            // {
            //     return -1;
            // }
        }
    }

    SECTION("linemesh test 4")
    {
        for( int i = 31; i <= 40; i++ )
        {
            std::vector<LineMesh> segements;
            auto samples = sample_discreteCurveDimensionParams_CURVEID_dimension_deviation
                <discreteCurveDimensionParams, int, &discreteCurveDimensionParams::CURVEID, int,
                &discreteCurveDimensionParams::dimension, double, &discreteCurveDimensionParams::deviation>(
                    i, i, 1, 3, 10, 1, 0.001, 0.1, 0.001, 25, true
            );
            // auto samples =
            //     sampleDiscreteCurveDimension(i, i, 1, 3, 10, 1, 1, 1, 1, 5.0, 5.0, 5.0, 0.001, 0.1, 0.001, true, 25);
            test_line_mesh_gen2(Geo, samples, segements);
            REQUIRE(get_linemesh_size(segements) != 0);
            // if( get_linemesh_size(segements) == 0 )
            // {
            //     return -1;
            // }
        }
    }
}

TEST_CASE("surfacemesh test", "[test2]")
{
    TiGER::Geometry Geo;
    test_load_geometry("f6.step", Geo);

    if( get_geo_size(Geo) == 0 )
    {
        return ;
    }
    REQUIRE(get_geo_size(Geo) != 0);
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
                test_surf_mesh_gen(Geo.surfaces_[i], surface_args1, segements1, surfacemesh);
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
            //dimension_args1.setdeviation(0.002);
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
            test_surf_mesh_gen(Geo.surfaces_[i], surface_args1, segements1, surfacemesh);
            checkMeshQuality(surfacemesh);
        }
    }
}

TEST_CASE("tanh_test", "[test3]")
{
    TiGER::Geometry Geo;
    test_load_geometry("f6.step", Geo);
    if( get_geo_size(Geo) == 0 )
    {
        return ;
    }
    REQUIRE(get_geo_size(Geo) != 0);
    std::cout << "after testing load geo" << std::endl;

    SECTION("tanh test 1")
    {
        for( int i = 5; i < 50; i += 5 )
        {
            std::vector<LineMesh> segements;
            auto samples = sample_discreteCurveTanhParams_end_s_dimension
                <discreteCurveTanhParams, double, &discreteCurveTanhParams::end_s, int,
                &discreteCurveTanhParams::dimension>(
                    0.01 * i, 0.5 * i, 0.05 * i, i, i, 5, 7
            );
            // auto samples = sampleDiscreteCurveTanh(0, 0, 1, 0.01 * i, 0.5 * i, 0.05 * i, i, i, 5, 7);
            // for( const auto& s : samples )
            // {
            //     std::cout << samples.size();
            //     std::cout << "begin_s: " << s.begin_s << ", end_s: " << s.end_s << ", dimension: " << s.dimension
            //               << std::endl;
            // }
            test_line_mesh_gen3(Geo, samples, segements);
            REQUIRE(get_linemesh_size(segements) != 0);
            // if( get_linemesh_size(segements) == 0 )
            // {
            //     return -1;
            // }
        }
    }

    SECTION("tanh test 2")
    {
        for( int i = 5; i < 50; i += 5 )
        {
            std::vector<LineMesh> segements;
            auto samples = sample_discreteCurveTanhParams_begin_s_dimension
                <discreteCurveTanhParams, double, &discreteCurveTanhParams::begin_s, int,
                &discreteCurveTanhParams::dimension>(
                    0.01 * i, 0.5 * i, 0.05 * i, i, i, 5, 7
            );
            // auto samples = sampleDiscreteCurveTanh(0.01 * i, 0.5 * i, 0.05 * i, 0, 0, 1, i, i, 5, 7);
            // for( const auto& s : samples )
            // {
            //     std::cout << samples.size();
            //     std::cout << "begin_s: " << s.begin_s << ", end_s: " << s.end_s << ", dimension: " << s.dimension
            //               << std::endl;
            // }
            test_line_mesh_gen3(Geo, samples, segements);
            REQUIRE(get_linemesh_size(segements) != 0);
            // if( get_linemesh_size(segements) == 0 )
            // {
            //     return -1;
            // }
        }
    }

    SECTION("tanh test 3")
    {
        std::vector<LineMesh> segements;
        auto samples = sample_discreteCurveTanhParams_begin_s_end_s
            <discreteCurveTanhParams, double, &discreteCurveTanhParams::begin_s, double,
            &discreteCurveTanhParams::end_s>(
                0.01 * 20, 0.5 * 20, 0.05 * 20, 0.01 * 20, 0.5 * 20, 0.05 * 20, 7
            );
        // auto samples =
        //     sampleDiscreteCurveTanh(0.01 * 20, 0.5 * 20, 0.05 * 20, 0.01 * 20, 0.5 * 20, 0.05 * 20, 20, 20, 1, 7);
        // for( const auto& s : samples )
        // {
        //     std::cout << samples.size();
        //     std::cout << "begin_s: " << s.begin_s << ", end_s: " << s.end_s << ", dimension: " << s.dimension <<
        //     std::endl;
        // }
        test_line_mesh_gen3(Geo, samples, segements);
        REQUIRE(get_linemesh_size(segements) != 0);
        // if( get_linemesh_size(segements) == 0 )
        // {
        //     return -1;
        // }
    }
}

// TEST_CASE("quality test", "[test4]")
// {
//     TiGER::Geometry Geo;
//     test_load_geometry("f6.step", Geo);
//     // if( get_geo_size(Geo) == 0 )
//     // {
//     //     return ;
//     // }
//     REQUIRE(get_geo_size(Geo) != 0);
//     std::cout << "after testing load geo" << std::endl;

//     std::cout << "面网格测试\n";
//     // 先保存所有
//     std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> all_curves_surfaces_(Geo.curves_.size());
//     for( int i = 0; i < Geo.topo_surf_to_curves_.size(); i++ )
//     {
//         for( int j = 0; j < Geo.topo_surf_to_curves_[i].size(); j++ )
//         {
//             all_curves_surfaces_[Geo.topo_surf_to_curves_[i][j]].push_back(Geo.surfaces_[i]);
//         }
//     }

//     SECTION("quality_test mMaxdeviation")
//     {
//         for( int i = 10; i <= 11; i++ )
//         {
//             std::cout << "face i:" << i << std::endl;
//             std::vector<LineMesh> segements1;
//             CurveParametersDimension dimension_args1;
//             // std::cout << dimension_args1.getangle() << std::endl;
//             dimension_args1.setdimension(3);
//             dimension_args1.setdelta_s(100);
//             // std::cout << mMaxDeviation1;
//             dimension_args1.setdeviation(0.002);
//             // // dimension_args1.setangle(0);
//             dimension_args1.setuse_database_curvature(true);
//             ////////////////////////////////////////////////////////
//             std::vector<std::shared_ptr<GeometryCurve>> curves;
//             // std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_(
//             //     Geo.topo_surf_to_curves_[i].size());
//             std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_;
//             for( auto curve_id : Geo.topo_surf_to_curves_[i] )
//             {
//                 curves.push_back(Geo.curves_[curve_id]);
//                 surfaces_.push_back(all_curves_surfaces_[curve_id]);
//             }
//             std::cout << "before line mesh\n";
//             test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
//             double diagonal = getDiagonalLength(segements1);

//             double mUserSpecifiedMaxEdgeLen1 = 100;
//             double mUserSpecifiedMinEdgeLen1 = 0.0;
//             double mMaxDeviation1 = diagonal / 500;
//             std::cout << "before mMaxDeiation:" << mMaxDeviation1 << std::endl;
//             // std::cout << mMaxDeviation1 << std::endl;

//             IsotropicSurfaceParametersTri surface_args1;
//             surface_args1.setmUserSpecifiedMaxEdgeLen(mUserSpecifiedMaxEdgeLen1);
//             surface_args1.setmUserSpecifiedMinEdgeLen(mUserSpecifiedMinEdgeLen1);
//             surface_args1.setmMaxDeviation(mMaxDeviation1);
//             surface_args1.setmSwapCellWithNoInteriorPoints(false);
//             SurfaceMesh surfacemesh;
//             test_surf_mesh_gen(Geo.surfaces_[i], surface_args1, segements1, surfacemesh);
//             std::cout << "first mMaxDeiation:" << mMaxDeviation1 << std::endl;
//             std::cout << "first surface size:" << get_surfacemesh_size(surfacemesh) << std::endl;
//             mMaxDeviation1 /= 2;
//             surface_args1.setmMaxDeviation(mMaxDeviation1);
//             SurfaceMesh surfacemesh2;
//             test_surf_mesh_gen(Geo.surfaces_[i], surface_args1, segements1, surfacemesh2);
//             std::cout << "second mMaxDeiation:" << mMaxDeviation1 << std::endl;
//             std::cout << "second surface size:" << get_surfacemesh_size(surfacemesh2) << std::endl;

//             mMaxDeviation1 /= 10;
//             surface_args1.setmMaxDeviation(mMaxDeviation1);
//             SurfaceMesh surfacemesh3;
//             test_surf_mesh_gen(Geo.surfaces_[i], surface_args1, segements1, surfacemesh3);
//             std::cout << "third mMaxDeiation:" << mMaxDeviation1 << std::endl;
//             std::cout << "third surface size:" << get_surfacemesh_size(surfacemesh3) << std::endl;
//             REQUIRE(get_surfacemesh_size(surfacemesh) <= get_surfacemesh_size(surfacemesh2));
//             REQUIRE(get_surfacemesh_size(surfacemesh2) <= get_surfacemesh_size(surfacemesh3));

//             // if( !(get_surfacemesh_size(surfacemesh) <= get_surfacemesh_size(surfacemesh2)) ||
//             //     !(get_surfacemesh_size(surfacemesh2) <= get_surfacemesh_size(surfacemesh3)) )
//             // {
//             //     return -1;
//             // }
//         }
//     }

//     SECTION("quality_test mMaxAngle")
//     {
//         for( int i = 10; i <= 11; i++ )
//         {
//             std::cout << "face i:" << i << std::endl;
//             std::vector<LineMesh> segements1;
//             CurveParametersDimension dimension_args1;
//             // std::cout << dimension_args1.getangle() << std::endl;
//             dimension_args1.setdimension(3);
//             dimension_args1.setdelta_s(100);
//             // std::cout << mMaxDeviation1;
//             dimension_args1.setdeviation(0.002);
//             // // dimension_args1.setangle(0);
//             dimension_args1.setuse_database_curvature(true);
//             ////////////////////////////////////////////////////////
//             std::vector<std::shared_ptr<GeometryCurve>> curves;
//             // std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_(
//             //     Geo.topo_surf_to_curves_[i].size());
//             std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_;
//             for( auto curve_id : Geo.topo_surf_to_curves_[i] )
//             {
//                 curves.push_back(Geo.curves_[curve_id]);
//                 surfaces_.push_back(all_curves_surfaces_[curve_id]);
//             }
//             std::cout << "before line mesh\n";
//             test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
//             double diagonal = getDiagonalLength(segements1);

//             double mUserSpecifiedMaxEdgeLen1 = 100;
//             double mUserSpecifiedMinEdgeLen1 = 0.2;
//             double mMaxDeviation1 = diagonal / 500;
//             std::cout << "before mMaxDeiation:" << mMaxDeviation1 << std::endl;
//             // std::cout << mMaxDeviation1 << std::endl;

//             IsotropicSurfaceParametersTri surface_args1;
//             surface_args1.setmUserSpecifiedMaxEdgeLen(mUserSpecifiedMaxEdgeLen1);
//             surface_args1.setmUserSpecifiedMinEdgeLen(mUserSpecifiedMinEdgeLen1);
//             surface_args1.setmMaxDeviation(mMaxDeviation1);
//             surface_args1.setmSwapCellWithNoInteriorPoints(false);
//             surface_args1.setmMaxAngle(40);
//             SurfaceMesh surfacemesh;
//             test_surf_mesh_gen(Geo.surfaces_[i], surface_args1, segements1, surfacemesh);
//             std::cout << "first Angle:" << 40 << std::endl;
//             std::cout << "first surface size:" << get_surfacemesh_size(surfacemesh) << std::endl;
//             SurfaceMesh surfacemesh2;
//             test_surf_mesh_gen(Geo.surfaces_[i], surface_args1, segements1, surfacemesh2);
//             std::cout << "second Angle:" << 20 << std::endl;
//             std::cout << "second surface size:" << get_surfacemesh_size(surfacemesh2) << std::endl;

//             SurfaceMesh surfacemesh3;
//             test_surf_mesh_gen(Geo.surfaces_[i], surface_args1, segements1, surfacemesh3);
//             std::cout << "third Angle:" << 10 << std::endl;
//             std::cout << "third surface size:" << get_surfacemesh_size(surfacemesh3) << std::endl;

//             REQUIRE(get_surfacemesh_size(surfacemesh) >= get_surfacemesh_size(surfacemesh2));
//             REQUIRE(get_surfacemesh_size(surfacemesh2) >= get_surfacemesh_size(surfacemesh3));
//             // if( !(get_surfacemesh_size(surfacemesh) <= get_surfacemesh_size(surfacemesh2)) ||
//             //     !(get_surfacemesh_size(surfacemesh2) <= get_surfacemesh_size(surfacemesh3)) )
//             // {
//             //     return -1;
//             // }
//         }
//     }

//     SECTION("quality_test linemesh angle")
//     {
//         for( int i = 0; i <= 53; i++ )
//         {
//             std::cout << "face i:" << i << std::endl;
//             std::vector<LineMesh> segements1;
//             CurveParametersDimension dimension_args1;
//             // std::cout << dimension_args1.getangle() << std::endl;
//             dimension_args1.setdimension(3);
//             dimension_args1.setdelta_s(100);
//             // std::cout << mMaxDeviation1;
//             dimension_args1.setdeviation(0.002);

//             dimension_args1.setuse_database_curvature(true);
//             ////////////////////////////////////////////////////////
//             std::vector<std::shared_ptr<GeometryCurve>> curves;
//             // std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_(
//             //     Geo.topo_surf_to_curves_[i].size());
//             std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_;
//             for( auto curve_id : Geo.topo_surf_to_curves_[i] )
//             {
//                 curves.push_back(Geo.curves_[curve_id]);
//                 surfaces_.push_back(all_curves_surfaces_[curve_id]);
//             }
//             dimension_args1.setangle(40);
//             std::cout << "angle" << 40 << std::endl;
//             test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
//             int lineMeshSize1 = get_linemesh_size(segements1);
//             dimension_args1.setangle(20);
//             std::cout << "angle" << 20 << std::endl;
//             test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
//             int lineMeshSize2 = get_linemesh_size(segements1);
//             dimension_args1.setangle(10);
//             std::cout << "angle" << 10 << std::endl;
//             test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
//             int lineMeshSize3 = get_linemesh_size(segements1);
//             std::cout << lineMeshSize1 << ' ' << lineMeshSize2 << ' ' << lineMeshSize3;
//             REQUIRE(lineMeshSize1 <= lineMeshSize2);
//             REQUIRE(lineMeshSize2 <= lineMeshSize3);
//             // if( !(lineMeshSize1 <= lineMeshSize2 && lineMeshSize2 <= lineMeshSize3) )
//             // {
//             //     return -1;
//             // }
//         }
//     }

//     SECTION("quality_test linemesh deviation")
//     {
//         for( int i = 10; i <= 11; i++ )
//         {
//             std::cout << "face i:" << i << std::endl;
//             std::vector<LineMesh> segements1;
//             CurveParametersDimension dimension_args1;
//             // std::cout << dimension_args1.getangle() << std::endl;
//             dimension_args1.setdimension(3);
//             dimension_args1.setdelta_s(100);
//             // std::cout << mMaxDeviation1;

//             dimension_args1.setangle(40);
//             dimension_args1.setuse_database_curvature(true);
//             ////////////////////////////////////////////////////////
//             std::vector<std::shared_ptr<GeometryCurve>> curves;
//             // std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_(
//             //     Geo.topo_surf_to_curves_[i].size());
//             std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_;
//             for( auto curve_id : Geo.topo_surf_to_curves_[i] )
//             {
//                 // std::cout << "curve_id" << curve_id << std::endl;
//                 curves.push_back(Geo.curves_[curve_id]);
//                 surfaces_.push_back(all_curves_surfaces_[curve_id]);
//                 break;
//             }
//             // 第一轮
//             double deviation1 = curves[0]->getTotalLenthFunction() / 1000;
//             std::cout << "deviation" << deviation1 << std::endl;
//             dimension_args1.setdeviation(deviation1);
//             test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
//             int lineMeshSize1 = get_linemesh_size(segements1);

//             // 第二轮
//             deviation1 /= 2;
//             dimension_args1.setdeviation(deviation1);
//             std::cout << "deviation" << deviation1 << std::endl;
//             test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
//             int lineMeshSize2 = get_linemesh_size(segements1);
//             // 第三轮
//             deviation1 /= 2;
//             dimension_args1.setdeviation(deviation1);
//             std::cout << "deviation" << deviation1 << std::endl;
//             test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
//             int lineMeshSize3 = get_linemesh_size(segements1);
//             std::cout << lineMeshSize1 << ' ' << lineMeshSize2 << ' ' << lineMeshSize3;

//             REQUIRE(lineMeshSize1 <= lineMeshSize2);
//             REQUIRE(lineMeshSize2 <= lineMeshSize3);
//             // if( !(lineMeshSize1 <= lineMeshSize2 && lineMeshSize2 <= lineMeshSize3) )
//             // {
//             //     return -1;
//             // }
//         }
//     }

//     SECTION("quality_test linemesh dimension")
//     {
//         for( int i = 0; i <= 53; i++ )
//         {
//             std::cout << "face i:" << i << std::endl;
//             std::vector<LineMesh> segements1;
//             CurveParametersDimension dimension_args1;
//             dimension_args1.setuse_database_curvature(true);
//             ////////////////////////////////////////////////////////
//             std::vector<std::shared_ptr<GeometryCurve>> curves;
//             // std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_(
//             //     Geo.topo_surf_to_curves_[i].size());
//             std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_;
//             for( auto curve_id : Geo.topo_surf_to_curves_[i] )
//             {
//                 // std::cout << "curve_id" << curve_id << std::endl;
//                 curves.push_back(Geo.curves_[curve_id]);
//                 surfaces_.push_back(all_curves_surfaces_[curve_id]);
//                 // break;
//             }
//             // 第一轮
//             dimension_args1.setdimension(3);
//             std::cout << "dimension" << 3 << std::endl;
//             test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
//             int lineMeshSize1 = get_linemesh_size(segements1);

//             // 第二轮
//             dimension_args1.setdimension(6);
//             std::cout << "dimension" << 6 << std::endl;
//             test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
//             int lineMeshSize2 = get_linemesh_size(segements1);
//             // 第三轮
//             dimension_args1.setdimension(12);
//             std::cout << "dimension" << 12 << std::endl;
//             test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
//             int lineMeshSize3 = get_linemesh_size(segements1);
//             std::cout << lineMeshSize1 << ' ' << lineMeshSize2 << ' ' << lineMeshSize3;

//             REQUIRE(lineMeshSize1 <= lineMeshSize2);
//             REQUIRE(lineMeshSize2 <= lineMeshSize3);
//             // if( !(lineMeshSize1 <= lineMeshSize2 && lineMeshSize2 <= lineMeshSize3) )
//             // {
//             //     return -1;
//             // }
//         }
//     }

//     SECTION("quality_test linemesh delta_s")
//     {
//         for( int i = 0; i <= 53; i++ )
//         {
//             std::cout << "face i:" << i << std::endl;
//             std::vector<LineMesh> segements1;
//             CurveParametersDimension dimension_args1;
//             dimension_args1.setuse_database_curvature(true);
//             ////////////////////////////////////////////////////////
//             std::vector<std::shared_ptr<GeometryCurve>> curves;
//             // std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_(
//             //     Geo.topo_surf_to_curves_[i].size());
//             std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_;
//             for( auto curve_id : Geo.topo_surf_to_curves_[i] )
//             {
//                 // std::cout << "curve_id" << curve_id << std::endl;
//                 curves.push_back(Geo.curves_[curve_id]);
//                 surfaces_.push_back(all_curves_surfaces_[curve_id]);
//                 break;
//             }
//             double linelen = curves[0]->getTotalLenthFunction();
//             double delta_s = linelen;
//             // 第一轮
//             dimension_args1.setdelta_s(delta_s);
//             std::cout << "delta_s" << delta_s << std::endl;
//             test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
//             int lineMeshSize1 = get_linemesh_size(segements1);

//             // 第二轮
//             delta_s /= 2;
//             dimension_args1.setdelta_s(delta_s);
//             std::cout << "delta_s" << delta_s << std::endl;
//             test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
//             int lineMeshSize2 = get_linemesh_size(segements1);
//             // 第三轮
//             delta_s /= 2;
//             dimension_args1.setdelta_s(delta_s);
//             std::cout << "delta_s" << delta_s << std::endl;
//             test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
//             int lineMeshSize3 = get_linemesh_size(segements1);
//             std::cout << lineMeshSize1 << ' ' << lineMeshSize2 << ' ' << lineMeshSize3;

//             REQUIRE(lineMeshSize1 <= lineMeshSize2);
//             REQUIRE(lineMeshSize2 <= lineMeshSize3);
//             // if( !(lineMeshSize1 <= lineMeshSize2 && lineMeshSize2 <= lineMeshSize3) )
//             // {
//             //     return -1;
//             // }
//         }
//     }

//     SECTION("quality_test surfacemesh")
//     {
//         for( int i = 0; i <= 53; i++ )
//         {
//             std::cout << "face i:" << i << std::endl;
//             std::vector<LineMesh> segements1;
//             CurveParametersDimension dimension_args1;
//             // std::cout << dimension_args1.getangle() << std::endl;
//             dimension_args1.setdimension(10);
//             dimension_args1.setdelta_s(100);
//             // std::cout << mMaxDeviation1;
//             dimension_args1.setdeviation(0.002);
//             // // dimension_args1.setangle(0);
//             dimension_args1.setuse_database_curvature(true);
//             ////////////////////////////////////////////////////////
//             std::vector<std::shared_ptr<GeometryCurve>> curves;
//             // std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_(
//             //     Geo.topo_surf_to_curves_[i].size());
//             std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_;
//             for( auto curve_id : Geo.topo_surf_to_curves_[i] )
//             {
//                 curves.push_back(Geo.curves_[curve_id]);
//                 surfaces_.push_back(all_curves_surfaces_[curve_id]);
//             }
//             std::cout << "before line mesh\n";
//             bool judge_line_mesh_gen = test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);

//             REQUIRE(judge_line_mesh_gen);
//             // if( !judge_line_mesh_gen )
//             // {
//             //     return -1;
//             // }
//             double mUserSpecifiedMaxEdgeLen1 = 100;
//             double mUserSpecifiedMinEdgeLen1 = 0.2;
//             double mMaxDeviation1 = 0.002;
//             // std::cout << mMaxDeviation1 << std::endl;

//             IsotropicSurfaceParametersTri surface_args1;
//             surface_args1.setmUserSpecifiedMaxEdgeLen(mUserSpecifiedMaxEdgeLen1);
//             surface_args1.setmUserSpecifiedMinEdgeLen(mUserSpecifiedMinEdgeLen1);
//             surface_args1.setmMaxDeviation(mMaxDeviation1);
//             surface_args1.setmSwapCellWithNoInteriorPoints(false);
//             SurfaceMesh surfacemesh;
//             test_surf_mesh_gen(Geo.surfaces_[i], surface_args1, segements1, surfacemesh);

//             REQUIRE(surfacemesh.tris.size() > 0);
//             // if( surfacemesh.tris.size() == 0 )
//             // {
//             //     return -1;
//             // }

//             vector<double> quality_value;
//             quality_value.push_back(1.00); // skewness
//             quality_value.push_back(2.72e-1);
//             quality_value.push_back(4.07e+5); // aspect
//             quality_value.push_back(8.06e+1);
//             quality_value.push_back(1.63e-4); // angle
//             quality_value.push_back(1.80e+2);
//             quality_value.push_back(3.05e+5); // condition
//             quality_value.push_back(6.38e+1);
//             // 之后改参数
//             // REQUIRE(checkMeshQuality(surfacemesh, quality_value));
//             // if( !checkMeshQuality(surfacemesh, quality_value) )
//             // {
//             //     return -1;
//             // }
//         }
//     }
// }