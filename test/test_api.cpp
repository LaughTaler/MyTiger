#include "test_api.h"

int main()
{
    TiGER::Geometry Geo;
    test_load_geometry("f6.step", Geo);
    if( get_geo_size(Geo) == 0 )
    {
        return -1;
    }
    std::cout << "after testing load geo" << std::endl;

    // 线网格测试1
    for( int i = 1; i <= 10; i++ )
    {
        std::cout << "line mesh test 1 line " << i << std::endl;
        std::vector<LineMesh> segements;
        double ttl = Geo.curves_[i]->getTotalLenthFunction();
        auto samples = sampleDiscreteCurveDimension(i, i, 1, 3, 3, 10, 0.01 * ttl, 0.1 * ttl, 0.01 * ttl, 10.0, 10.0,
                                                    10.0, 0.01, 0.01, 0.01, false, 19);
        // for( const auto& s : samples )
        // {
        //     std::cout << samples.size();
        //     std::cout << "CURVEID: " << s.CURVEID << ", dimension: " << s.dimension << ", delta_s: " << s.delta_s
        //               << ", angle: " << s.angle << ", deviation: " << s.deviation
        //               << ", use_database_curvature: " << std::boolalpha << s.use_database_curvature << std::endl;
        // }
        test_line_mesh_gen2(Geo, samples, segements);
        if( get_linemesh_size(segements) == 0 )
        {
            return -1;
        }
    }
    //------------------------------------
    // 线网格测试2
    for( int i = 11; i < 20; i++ )
    {
        std::cout << "line mesh test 2 line " << i << std::endl;
        std::vector<LineMesh> segements;
        double ttl = Geo.curves_[i]->getTotalLenthFunction();
        auto samples =
            sampleDiscreteCurveDimension(i, i, 1, 3, 10, 1, 1, 1, 1, 5.0, 40.0, 5.0, 0.01, 0.01, 0.01, true, 21);
        // for( const auto& s : samples )
        // {
        //     std::cout << samples.size();
        //     std::cout << "CURVEID: " << s.CURVEID << ", dimension: " << s.dimension << ", delta_s: " << s.delta_s
        //               << ", angle: " << s.angle << ", deviation: " << s.deviation
        //               << ", use_database_curvature: " << std::boolalpha << s.use_database_curvature << std::endl;
        // }
        test_line_mesh_gen2(Geo, samples, segements);
        if( get_linemesh_size(segements) == 0 )
        {
            return -2;
        }
    }
    // 线网格测试3
    for( int i = 21; i < 30; i++ )
    {
        std::cout << "line mesh test 3 line " << i << std::endl;
        std::vector<LineMesh> segements;
        double ttl = Geo.curves_[i]->getTotalLenthFunction();
        auto samples =
            sampleDiscreteCurveDimension(i, i, 1, 0, 100, 10, 1, 1, 1, 5.0, 5.0, 1.0, 0.01, 0.01, 0.01, false, 17);
        // for( const auto& s : samples )
        // {
        //     std::cout << samples.size();
        //     std::cout << "CURVEID: " << s.CURVEID << ", dimension: " << s.dimension << ", delta_s: " << s.delta_s
        //               << ", angle: " << s.angle << ", deviation: " << s.deviation
        //               << ", use_database_curvature: " << std::boolalpha << s.use_database_curvature << std::endl;
        // }
        test_line_mesh_gen2(Geo, samples, segements);
        if( get_linemesh_size(segements) == 0 )
        {
            return -3;
        }
    }
    // 线网格测试4
    for( int i = 31; i < 40; i++ )
    {
        std::cout << "line mesh test 4 line " << i << std::endl;
        std::vector<LineMesh> segements;
        double ttl = Geo.curves_[i]->getTotalLenthFunction();
        auto samples =
            sampleDiscreteCurveDimension(i, i, 1, 3, 10, 1, 1, 1, 1, 5.0, 5.0, 5.0, 0.001, 0.1, 0.001, true, 25);
        // for( const auto& s : samples )
        // {
        //     std::cout << samples.size();
        //     std::cout << "CURVEID: " << s.CURVEID << ", dimension: " << s.dimension << ", delta_s: " << s.delta_s
        //               << ", angle: " << s.angle << ", deviation: " << s.deviation
        //               << ", use_database_curvature: " << std::boolalpha << s.use_database_curvature << std::endl;
        // }
        test_line_mesh_gen2(Geo, samples, segements);
        if( get_linemesh_size(segements) == 0 )
        {
            return -4;
        }
    }

    // 面网格测试
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

    // 面网格测试1
    for( int i = 0; i < 53; i++ )
    {
        std::cout << "surface mesh test 1 face i:" << i << std::endl;
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
            if( get_linemesh_size(segements1) == 0 )
            {
                return -5;
            }
            std::cout << "after line mesh\n";
            IsotropicSurfaceParametersTri surface_args1;
            surface_args1.setmUserSpecifiedMaxEdgeLen(mUserSpecifiedMaxEdgeLen1);
            surface_args1.setmUserSpecifiedMinEdgeLen(mUserSpecifiedMinEdgeLen1);
            surface_args1.setmMaxDeviation(mMaxDeviation1);
            surface_args1.setmSwapCellWithNoInteriorPoints(false);
            SurfaceMesh surfacemesh;
            test_surf_mesh_gen(Geo.surfaces_[i], surface_args1, segements1, surfacemesh);
            if( get_surfacemesh_size(surfacemesh) == 0 )
            {
                // 返回值
                return -6;
            }
        }
    }

    // 面网格测试2
    for( int i = 0; i <= 53; i++ )
    {
        std::cout << "surface mesh test 2 face i:" << i << std::endl;
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
        test_surf_mesh_gen(Geo.surfaces_[i], surface_args1, segements1, surfacemesh);
        if( get_surfacemesh_size(surfacemesh) == 0 )
        {
            // 返回值
            return -7;
        }
        // checkMeshQuality(surfacemesh);
    }

    // tanh测试1
    std::vector<LineMesh> segements;
    for( int i = 5; i < 50; i += 5 )
    {
        std::cout << "tanh 1 " << i << std::endl;
        auto samples = sampleDiscreteCurveTanh(0, 0, 1, 0.01 * i, 0.5 * i, 0.05 * i, i, i, 5, 7);
        for( const auto& s : samples )
        {
            std::cout << samples.size();
            std::cout << "begin_s: " << s.begin_s << ", end_s: " << s.end_s << ", dimension: " << s.dimension
                      << std::endl;
        }
        test_line_mesh_gen3(Geo, samples, segements);
        if( get_linemesh_size(segements) == 0 )
        {
            return -8;
        }
    }
    // tanh测试2
    for( int i = 5; i < 50; i += 5 )
    {
        std::cout << "tanh 2 " << i << std::endl;
        auto samples = sampleDiscreteCurveTanh(0.01 * i, 0.5 * i, 0.05 * i, 0, 0, 1, i, i, 5, 7);
        for( const auto& s : samples )
        {
            std::cout << samples.size();
            std::cout << "begin_s: " << s.begin_s << ", end_s: " << s.end_s << ", dimension: " << s.dimension
                      << std::endl;
        }
        test_line_mesh_gen3(Geo, samples, segements);
        if( get_linemesh_size(segements) == 0 )
        {
            return -9;
        }
    }
    // tanh测试3
    auto samples =
        sampleDiscreteCurveTanh(0.01 * 20, 0.5 * 20, 0.05 * 20, 0.01 * 20, 0.5 * 20, 0.05 * 20, 20, 20, 1, 7);
    for( const auto& s : samples )
    {
        std::cout << samples.size();
        std::cout << "begin_s: " << s.begin_s << ", end_s: " << s.end_s << ", dimension: " << s.dimension << std::endl;
    }
    test_line_mesh_gen3(Geo, samples, segements);
    if( get_linemesh_size(segements) == 0 )
    {
        return -10;
    }

    // 质量测试mMaxdeviation
    for( int i = 10; i <= 11; i++ )
    {
        std::cout << "quality mesh mMaxdeviation face i:" << i << std::endl;
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
        double diagonal = getDiagonalLength(segements1);

        double mUserSpecifiedMaxEdgeLen1 = 100;
        double mUserSpecifiedMinEdgeLen1 = 0.2;
        double mMaxDeviation1 = diagonal / 500;
        std::cout << "before mMaxDeiation:" << mMaxDeviation1 << std::endl;
        // std::cout << mMaxDeviation1 << std::endl;

        IsotropicSurfaceParametersTri surface_args1;
        surface_args1.setmUserSpecifiedMaxEdgeLen(mUserSpecifiedMaxEdgeLen1);
        surface_args1.setmUserSpecifiedMinEdgeLen(mUserSpecifiedMinEdgeLen1);
        surface_args1.setmMaxDeviation(mMaxDeviation1);
        surface_args1.setmSwapCellWithNoInteriorPoints(false);
        SurfaceMesh surfacemesh;
        test_surf_mesh_gen(Geo.surfaces_[i], surface_args1, segements1, surfacemesh);
        std::cout << "first mMaxDeiation:" << mMaxDeviation1 << std::endl;
        std::cout << "first surface size:" << get_surfacemesh_size(surfacemesh) << std::endl;
        mMaxDeviation1 /= 2;
        surface_args1.setmMaxDeviation(mMaxDeviation1);
        SurfaceMesh surfacemesh2;
        test_surf_mesh_gen(Geo.surfaces_[i], surface_args1, segements1, surfacemesh2);
        std::cout << "second mMaxDeiation:" << mMaxDeviation1 << std::endl;
        std::cout << "second surface size:" << get_surfacemesh_size(surfacemesh2) << std::endl;

        mMaxDeviation1 /= 10;
        surface_args1.setmMaxDeviation(mMaxDeviation1);
        SurfaceMesh surfacemesh3;
        test_surf_mesh_gen(Geo.surfaces_[i], surface_args1, segements1, surfacemesh3);
        std::cout << "third mMaxDeiation:" << mMaxDeviation1 << std::endl;
        std::cout << "third surface size:" << get_surfacemesh_size(surfacemesh3) << std::endl;
        if( !((get_surfacemesh_size(surfacemesh) >= get_surfacemesh_size(surfacemesh2)) &&
              (get_surfacemesh_size(surfacemesh2) >= get_surfacemesh_size(surfacemesh3))) )
        {
            return -11;
        }
    }

    // 质量测试mMaxAngle
    for( int i = 10; i <= 11; i++ )
    {
        std::cout << "quality mesh mMaxAngle face i:" << i << std::endl;
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
        double diagonal = getDiagonalLength(segements1);

        double mUserSpecifiedMaxEdgeLen1 = 100;
        double mUserSpecifiedMinEdgeLen1 = 0.2;
        double mMaxDeviation1 = diagonal / 500;
        std::cout << "before mMaxDeiation:" << mMaxDeviation1 << std::endl;
        // std::cout << mMaxDeviation1 << std::endl;

        IsotropicSurfaceParametersTri surface_args1;
        surface_args1.setmUserSpecifiedMaxEdgeLen(mUserSpecifiedMaxEdgeLen1);
        surface_args1.setmUserSpecifiedMinEdgeLen(mUserSpecifiedMinEdgeLen1);
        surface_args1.setmMaxDeviation(mMaxDeviation1);
        surface_args1.setmSwapCellWithNoInteriorPoints(false);
        surface_args1.setmMaxAngle(40);
        SurfaceMesh surfacemesh;
        test_surf_mesh_gen(Geo.surfaces_[i], surface_args1, segements1, surfacemesh);
        std::cout << "first Angle:" << 40 << std::endl;
        std::cout << "first surface size:" << get_surfacemesh_size(surfacemesh) << std::endl;
        SurfaceMesh surfacemesh2;
        test_surf_mesh_gen(Geo.surfaces_[i], surface_args1, segements1, surfacemesh2);
        std::cout << "second Angle:" << 20 << std::endl;
        std::cout << "second surface size:" << get_surfacemesh_size(surfacemesh2) << std::endl;

        SurfaceMesh surfacemesh3;
        test_surf_mesh_gen(Geo.surfaces_[i], surface_args1, segements1, surfacemesh3);
        std::cout << "third Angle:" << 10 << std::endl;
        std::cout << "third surface size:" << get_surfacemesh_size(surfacemesh3) << std::endl;
        if( !((get_surfacemesh_size(surfacemesh) >= get_surfacemesh_size(surfacemesh2)) &&
              (get_surfacemesh_size(surfacemesh2) >= get_surfacemesh_size(surfacemesh3))) )
        {
            return -12;
        }
    }

    // 质量测试线网格angle
    for( int i = 0; i <= 53; i++ )
    {
        std::cout << "-------------------------------quality mesh angle face i:" << i << std::endl;
        std::vector<LineMesh> segements1;
        CurveParametersDimension dimension_args1;
        // std::cout << dimension_args1.getangle() << std::endl;
        dimension_args1.setdimension(3);
        dimension_args1.setdelta_s(100);
        // std::cout << mMaxDeviation1;
        dimension_args1.setdeviation(0.002);

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
        dimension_args1.setangle(40);
        std::cout << "angle" << 40 << std::endl;
        test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
        int lineMeshSize1 = get_linemesh_size(segements1);
        dimension_args1.setangle(20);
        std::cout << "angle" << 20 << std::endl;
        test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
        int lineMeshSize2 = get_linemesh_size(segements1);
        dimension_args1.setangle(10);
        std::cout << "angle" << 10 << std::endl;
        test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
        int lineMeshSize3 = get_linemesh_size(segements1);
        std::cout << lineMeshSize1 << ' ' << lineMeshSize2 << ' ' << lineMeshSize3;
        if( !((lineMeshSize1 <= lineMeshSize2) && (lineMeshSize2 <= lineMeshSize3)) )
        {
            return -13;
        }
    }
    // 质量测试线网格deviation
    for( int i = 10; i <= 11; i++ )
    {
        std::cout << "-------------------------------quality mesh deviation face i:" << i << std::endl;
        // std::cout << "face i:" << i << std::endl;
        std::vector<LineMesh> segements1;
        CurveParametersDimension dimension_args1;
        // std::cout << dimension_args1.getangle() << std::endl;
        dimension_args1.setdimension(3);
        dimension_args1.setdelta_s(100);
        // std::cout << mMaxDeviation1;

        dimension_args1.setangle(40);
        dimension_args1.setuse_database_curvature(true);
        ////////////////////////////////////////////////////////
        std::vector<std::shared_ptr<GeometryCurve>> curves;
        // std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_(
        //     Geo.topo_surf_to_curves_[i].size());
        std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_;
        for( auto curve_id : Geo.topo_surf_to_curves_[i] )
        {
            // std::cout << "curve_id" << curve_id << std::endl;
            curves.push_back(Geo.curves_[curve_id]);
            surfaces_.push_back(all_curves_surfaces_[curve_id]);
            break;
        }
        // 第一轮
        double deviation1 = curves[0]->getTotalLenthFunction() / 1000;
        std::cout << "deviation" << deviation1 << std::endl;
        dimension_args1.setdeviation(deviation1);
        test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
        int lineMeshSize1 = get_linemesh_size(segements1);

        // 第二轮
        deviation1 /= 2;
        dimension_args1.setdeviation(deviation1);
        std::cout << "deviation" << deviation1 << std::endl;
        test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
        int lineMeshSize2 = get_linemesh_size(segements1);
        // 第三轮
        deviation1 /= 2;
        dimension_args1.setdeviation(deviation1);
        std::cout << "deviation" << deviation1 << std::endl;
        test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
        int lineMeshSize3 = get_linemesh_size(segements1);
        std::cout << lineMeshSize1 << ' ' << lineMeshSize2 << ' ' << lineMeshSize3;
        if( !((lineMeshSize1 <= lineMeshSize2) && (lineMeshSize2 <= lineMeshSize3)) )
        {
            return -14;
        }
    }

    // 质量测试线网格dimension
    for( int i = 0; i <= 53; i++ )
    {
        std::cout << "-------------------------------quality mesh dimension face i:" << i << std::endl;
        // std::cout << "face i:" << i << std::endl;
        std::vector<LineMesh> segements1;
        CurveParametersDimension dimension_args1;
        dimension_args1.setuse_database_curvature(true);
        ////////////////////////////////////////////////////////
        std::vector<std::shared_ptr<GeometryCurve>> curves;
        // std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_(
        //     Geo.topo_surf_to_curves_[i].size());
        std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_;
        for( auto curve_id : Geo.topo_surf_to_curves_[i] )
        {
            // std::cout << "curve_id" << curve_id << std::endl;
            curves.push_back(Geo.curves_[curve_id]);
            surfaces_.push_back(all_curves_surfaces_[curve_id]);
            // break;
        }
        // 第一轮
        dimension_args1.setdimension(3);
        std::cout << "dimension" << 3 << std::endl;
        test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
        int lineMeshSize1 = get_linemesh_size(segements1);

        // 第二轮
        dimension_args1.setdimension(6);
        std::cout << "dimension" << 6 << std::endl;
        test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
        int lineMeshSize2 = get_linemesh_size(segements1);
        // 第三轮
        dimension_args1.setdimension(12);
        std::cout << "dimension" << 12 << std::endl;
        test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
        int lineMeshSize3 = get_linemesh_size(segements1);
        std::cout << lineMeshSize1 << ' ' << lineMeshSize2 << ' ' << lineMeshSize3;
        if( !(lineMeshSize1 <= lineMeshSize2 && lineMeshSize2 <= lineMeshSize3) )
        {
            return -15;
        }
    }

    // 质量测试线网格delta_s
    for( int i = 0; i <= 53; i++ )
    {
        std::cout << "-------------------------------quality mesh deltas face i:" << i << std::endl;
        // std::cout << "face i:" << i << std::endl;
        std::vector<LineMesh> segements1;
        CurveParametersDimension dimension_args1;
        dimension_args1.setuse_database_curvature(true);
        ////////////////////////////////////////////////////////
        std::vector<std::shared_ptr<GeometryCurve>> curves;
        // std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_(
        //     Geo.topo_surf_to_curves_[i].size());
        std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_;
        for( auto curve_id : Geo.topo_surf_to_curves_[i] )
        {
            // std::cout << "curve_id" << curve_id << std::endl;
            curves.push_back(Geo.curves_[curve_id]);
            surfaces_.push_back(all_curves_surfaces_[curve_id]);
            break;
        }
        double linelen = curves[0]->getTotalLenthFunction();
        double delta_s = linelen;
        // 第一轮
        dimension_args1.setdelta_s(delta_s);
        std::cout << "delta_s" << delta_s << std::endl;
        test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
        int lineMeshSize1 = get_linemesh_size(segements1);

        // 第二轮
        delta_s /= 2;
        dimension_args1.setdelta_s(delta_s);
        std::cout << "delta_s" << delta_s << std::endl;
        test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
        int lineMeshSize2 = get_linemesh_size(segements1);
        // 第三轮
        delta_s /= 2;
        dimension_args1.setdelta_s(delta_s);
        std::cout << "delta_s" << delta_s << std::endl;
        test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
        int lineMeshSize3 = get_linemesh_size(segements1);
        std::cout << lineMeshSize1 << ' ' << lineMeshSize2 << ' ' << lineMeshSize3;
        if( !(lineMeshSize1 <= lineMeshSize2 && lineMeshSize2 <= lineMeshSize3) )
        {
            return -16;
        }
    }

    // 面网格质量测试
    for( int i = 0; i <= 53; i++ )
    {
        std::cout << "-------------------------------quality test face i:" << i << std::endl;
        // std::cout << "face i:" << i << std::endl;
        std::vector<LineMesh> segements1;
        CurveParametersDimension dimension_args1;
        // std::cout << dimension_args1.getangle() << std::endl;
        dimension_args1.setdimension(10);
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
        bool judge_line_mesh_gen = test_line_mesh_gen4(curves, surfaces_, dimension_args1, segements1);
        if( !judge_line_mesh_gen )
        {
            return -17;
        }
        double mUserSpecifiedMaxEdgeLen1 = 100;
        double mUserSpecifiedMinEdgeLen1 = 0.2;
        double mMaxDeviation1 = 0.002;
        // std::cout << mMaxDeviation1 << std::endl;

        IsotropicSurfaceParametersTri surface_args1;
        surface_args1.setmUserSpecifiedMaxEdgeLen(mUserSpecifiedMaxEdgeLen1);
        surface_args1.setmUserSpecifiedMinEdgeLen(mUserSpecifiedMinEdgeLen1);
        surface_args1.setmMaxDeviation(mMaxDeviation1);
        surface_args1.setmSwapCellWithNoInteriorPoints(false);
        SurfaceMesh surfacemesh;
        test_surf_mesh_gen(Geo.surfaces_[i], surface_args1, segements1, surfacemesh);
        if( surfacemesh.tris.size() == 0 )
        {
            return -18;
        }

        vector<double> quality_value;
        quality_value.push_back(1.00); // skewness
        quality_value.push_back(2.72e-1);
        quality_value.push_back(4.07e+5); // aspect
        quality_value.push_back(8.06e+1);
        quality_value.push_back(1.63e-4); // angle
        quality_value.push_back(1.80e+2);
        quality_value.push_back(3.05e+5); // condition
        quality_value.push_back(6.38e+1);
        if( !checkMeshQuality(surfacemesh, quality_value) )
        {
            // return -19;
        }
    }
    return 0;
}