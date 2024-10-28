// test_load_geometry(filename,geo)
// test_line_mesh_gen(geo,linemesh)
// test_surf_mesh_gen(geo,linemesh,surfacemesh)
// test_mesh_size(surfacemesh,lower_bound,upper_bound)
// test_meshing_time(surfacemesh,lower_bound,upper_bound)
// test_mesh_quality(surfacemesh,lower_bound,upper_bound)
// test_mesh_para(geo,#para,#lower_bound,#upper_bound)

// #pragma once

#include "discretize_surface_STL.h"
#include "discretize_surface.h"
#include "discretize_surface_parameters.h"
#include "discretize_curve.h"
#include "discretize_curve_parameters.h"
#include "geometry_data_struct.h"
#include "geometry_interface_occ.h"
#include "geometry_interface_parameters.h"
#include "sizing_field_interface.h"
#include "geometry_interface_plane.h"
#include "../tools/meshIO.h"
#include "../tools/output_info.hpp"
#include "../tools/test_linemesh.hpp"
#include "catch.hpp"
#include "matrix_to_list.hpp"
#include "list_to_matrix.hpp"
#include "merge_surface_mesh.hpp"
#include "meshIO.h"
#include "tetrahedral_mesh.h"
#include "tetrahedral_mesh_parameters.h"
#include "alias.h"
// #include "MeshMetric.h"
#include <chrono>
#include <thread>
#include <iostream>
#include "mesh_quality.h"
#include <cfloat> // For DBL_MAX and DBL_MIN
#include <functional>
#include <cstdio>
#include <memory>
#include <array>
#include <string>
// #include "spdlog/spdlog.h"
// #include "spdlog/pattern_formatter.h"
bool test_load_geometry(std::string filename, TiGER::Geometry& Geo)
{
#ifdef _MSC_VER
    filename = "Z:/home/sftp1/models/" + filename;
    // filename = "D:/tiger2/test/f6_cf/" + filename;
#else
    filename = "/home/sftp1/models/" + filename;
#endif
    TiGER::readModelParameters readmodel_args;
    std::cout << filename;
    TiGER::Read_Model(filename, Geo, readmodel_args);
    return true; // 假设函数总是返回true，如果有错误处理逻辑请根据实际情况调整
}

bool test_load_vtk(std::string filename, TiGER::Mesh& m)
{
#ifdef _MSC_VER
    filename = "Z:/home/sftp1/models/" + filename;
    // filename = "D:/tiger2/test/f6_cf/" + filename;
#else
    filename = "/home/sftp1/models/" + filename;
#endif
    TiGER::readModelParameters readmodel_args;
    std::cout << filename;
    TiGER::MESHIO::readVTK(filename, m);
    return true; // 假设函数总是返回true，如果有错误处理逻辑请根据实际情况调整
}

bool test_line_mesh_gen(TiGER::Geometry& Geo, CurveParametersDimension& dimension_args,
                        std::vector<LineMesh>& segements)
{

    std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_(Geo.curves_.size());
    for( int i = 0; i < Geo.topo_surf_to_curves_.size(); i++ )
    {
        for( int j = 0; j < Geo.topo_surf_to_curves_[i].size(); j++ )
        {
            surfaces_[Geo.topo_surf_to_curves_[i][j]].push_back(Geo.surfaces_[i]);
        }
    }
    TiGER::discretize_curve::Grid1D_Dimension(Geo.curves_, surfaces_, dimension_args, segements);
    return true;
}
bool test_line_mesh_gen3(TiGER::Geometry& Geo, std::vector<std::shared_ptr<GeometryCurve>>& curves,
                         CurveParametersDimension& dimension_args, std::vector<LineMesh>& segements)
{
    std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_(Geo.curves_.size());

    for( int i = 0; i < curves.size(); i++ )
    {
        surfaces_[i].push_back(Geo.surfaces_[i]);
    }
    TiGER::discretize_curve::Grid1D_Dimension(curves, surfaces_, dimension_args, segements);
    // std::cout << segements.size() << std::endl;
    // std::cout << segements[0].coord.size() << std::endl;
    // std::cout << segements[1].coord.size() << std::endl;
    // std::cout << segements[2].coord.size() << std::endl;
    // std::cout << segements[3].coord.size() << std::endl;
    // std::cout << segements[4].coord.size() << std::endl;
    return true;
}
bool test_line_mesh_gen4(std::vector<std::shared_ptr<GeometryCurve>>& curves,
                         std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>>& surfaces_,
                         CurveParametersDimension& dimension_args, std::vector<LineMesh>& segements)
{
    TiGER::discretize_curve::Grid1D_Dimension(curves, surfaces_, dimension_args, segements);
    // std::cout << segements.size() << std::endl;
    // std::cout << segements[0].coord.size() << std::endl;
    // std::cout << segements[1].coord.size() << std::endl;
    // std::cout << segements[2].coord.size() << std::endl;
    // std::cout << segements[3].coord.size() << std::endl;
    // std::cout << segements[4].coord.size() << std::endl;
    return true;
}
bool test_surf_mesh_gen2(TiGER::Geometry& Geo, IsotropicSurfaceParametersTri& surface_args,
                         std::vector<LineMesh>& linemesh, SurfaceMesh& surfacemesh)
{
    // args.setangle(dimension_angle);
    // args.setdelta_s(dimension_delta_s);
    // args.setdimension(dimension_dimension);
    // args.setmax_dimension(dimension_maxdimension);

    SizingFunction size_field = [](const double& x, const double& y, const double& z) { return 1; };
    std::vector<std::vector<LineMesh>> sf_boundary;
    sf_boundary.resize(Geo.surfaces_.size());
    std::cout << Geo.topo_surf_to_curves_.size() << std::endl;
    std::cout << "before for" << std::endl;
    for( int i = 0; i < Geo.topo_surf_to_curves_.size(); i++ )
    {
        // std::cout << Geo.topo_surf_to_curves_[i].size() << std::endl;
        for( int j = 0; j < Geo.topo_surf_to_curves_[i].size(); j++ )
        {
            sf_boundary[i].push_back(linemesh[Geo.topo_surf_to_curves_[i][j]]);
        }
    }
    std::vector<SurfaceMesh> sf_mesh(sf_boundary.size());
    std::cout << "after sf_mesh" << std::endl;
    for( int k = 0; k < sf_boundary.size(); k++ )
    {

        std::cout << "jisuan k" << k << std::endl;
        TiGER::discretize_surface::CADSurf_TGrid_AFT(Geo.surfaces_[k], sf_boundary[k], surface_args, size_field,
                                                        sf_mesh[k]);
    }
    SurfaceMesh sf_output;
    TiGER::merge_surface_mesh(sf_mesh, sf_output);
    surfacemesh = sf_output;
    return true;
}
bool test_surf_mesh_gen(std::shared_ptr<TiGER::GeometrySurface>& GeoSurface,
                        IsotropicSurfaceParametersTri& surface_args, std::vector<LineMesh>& linemesh,
                        SurfaceMesh& surfacemesh)
{
    SizingFunction size_field = [](const double& x, const double& y, const double& z) { return 1; };
    std::cout << "before for" << std::endl;
    TiGER::discretize_surface::CADSurf_TGrid_AFT(GeoSurface, linemesh, surface_args, size_field, surfacemesh);
    std::cout << "after sf_mesh" << std::endl;
    return true;
}

bool test_surf_mesh_gen_dt(std::shared_ptr<TiGER::GeometrySurface>& GeoSurface,
                           IsotropicSurfaceParametersTri& surface_args, std::vector<LineMesh>& linemesh,
                           SurfaceMesh& surfacemesh)
{
    SizingFunction size_field = [](const double& x, const double& y, const double& z) { return 1; };
    std::cout << "before for" << std::endl;
    TiGER::discretize_surface::CADSurf_exportToSTL(GeoSurface, linemesh, surface_args, size_field, surfacemesh);
    std::cout << "after sf_mesh" << std::endl;
    return true;
}

bool test_surfacemesh_size(SurfaceMesh& surfacemesh, int lower_bound, int upper_bound)
{
    int sumSize = surfacemesh.tris.size() + surfacemesh.quads.size();

    if( sumSize >= lower_bound && sumSize <= upper_bound )
    {
        return true;
    }
    return false;
}
// 返回几何中弧线数量
int get_geo_size(TiGER::Geometry& Geo)
{

    return Geo.curves_.size();
}

int get_linemesh_size(std::vector<LineMesh>& linemesh)
{
    int summ = 0;
    for( int i = 0; i < linemesh.size(); i++ )
    {
        summ += linemesh[i].coord.size();
    }
    return summ;
}
int get_surfacemesh_size(SurfaceMesh& surfacemesh)
{
    return surfacemesh.tris.size() + surfacemesh.quads.size();
}
bool isLargerMesh(SurfaceMesh& surfacemesh1, SurfaceMesh& surfacemesh2)
{
    int sumSize1 = surfacemesh1.tris.size() + surfacemesh1.quads.size();
    int sumSize2 = surfacemesh2.tris.size() + surfacemesh2.quads.size();
    if( sumSize1 > sumSize2 )
    {
        return true;
    }
    return false;
}

struct discreteCurveDimensionParams
{
    int CURVEID;
    int dimension;
    double delta_s;
    double angle;
    double deviation;
    bool use_database_curvature;
    int mask;

    discreteCurveDimensionParams(int CURVEID_val = 0, int dimension_val = 3, double delta_s_val = 1,
                                 double angle_val = 5.0, double deviation_val = 0.01,
                                 bool use_database_curvature_val = false)
        : CURVEID(CURVEID_val)
        , dimension(dimension_val)
        , delta_s(delta_s_val)
        , angle(angle_val)
        , deviation(deviation_val)
        , use_database_curvature(use_database_curvature_val)
    {
    }

    discreteCurveDimensionParams(bool judge, int CURVEID_val = 0, int dimension_val = 3, double delta_s_val = 1,
                                 double angle_val = 5.0, double deviation_val = 0.01,
                                 bool use_database_curvature_val = true)
        : CURVEID(CURVEID_val)
        , dimension(dimension_val)
        , delta_s(delta_s_val)
        , angle(angle_val)
        , deviation(deviation_val)
        , use_database_curvature(use_database_curvature_val)
    {
    }
};
#define SAMPLE_ASTRUCT(value1, value2, value3, value4, value5, value6)                                                 \
    template <typename T1, typename T2, typename T3, typename T4, typename T5>                                         \
    std::vector<discreteCurveDimensionParams> sample_##value1##_##value2##_##value3##_##value4##_##value5##_##value6(  \
        T1 min1, T1 max1, T1 distance1, T2 min2, T2 max2, T2 distance2, T3 min3, T3 max3, T3 distance3, T4 min4,       \
        T4 max4, T4 distance4, T5 min5, T5 max5, T5 distance5, bool includeTrueInBoolParam, int mask)                  \
    {                                                                                                                  \
        std::vector<discreteCurveDimensionParams> ans;                                                                 \
        for( T1 i = min1; i <= max1; i += distance1 )                                                                  \
        {                                                                                                              \
            for( T2 j = min2; j <= max2; j += distance2 )                                                              \
            {                                                                                                          \
                for( T3 k = min3; k <= max3; k += distance3 )                                                          \
                {                                                                                                      \
                    for( T4 l = min4; l <= max4; l += distance4 )                                                      \
                    {                                                                                                  \
                        for( T5 m = min5; m <= max5; m += distance5 )                                                  \
                        {                                                                                              \
                            if( includeTrueInBoolParam )                                                               \
                            {                                                                                          \
                                for( bool n : {false, true} )                                                          \
                                {                                                                                      \
                                    discreteCurveDimensionParams a;                                                    \
                                    a.value1 = i;                                                                      \
                                    a.value2 = j;                                                                      \
                                    a.value3 = k;                                                                      \
                                    a.value4 = l;                                                                      \
                                    a.value5 = m;                                                                      \
                                    a.value6 = n;                                                                      \
                                    a.mask = mask;                                                                     \
                                    ans.push_back(a);                                                                  \
                                }                                                                                      \
                            }                                                                                          \
                            else                                                                                       \
                            {                                                                                          \
                                discreteCurveDimensionParams a;                                                        \
                                a.value1 = i;                                                                          \
                                a.value2 = j;                                                                          \
                                a.value3 = k;                                                                          \
                                a.value4 = l;                                                                          \
                                a.value5 = m;                                                                          \
                                a.value6 = false;                                                                      \
                                a.mask = mask;                                                                         \
                                ans.push_back(a);                                                                      \
                            }                                                                                          \
                        }                                                                                              \
                    }                                                                                                  \
                }                                                                                                      \
            }                                                                                                          \
        }                                                                                                              \
        return ans;                                                                                                    \
    }
SAMPLE_ASTRUCT(CURVEID, dimension, delta_s, angle, deviation, use_database_curvature)
std::vector<discreteCurveDimensionParams> sampleDiscreteCurveDimension(
    int curveidMin, int curveidMax, int curveidDistance, int dimensionMin, int dimensionMax, int dimensionDistance,
    double delta_sMin, double delta_sMax, double delta_sDistance, double angleMin, double angleMax,
    double angleDistance, double deviationMin, double deviationMax, double deviationDistance,
    bool includeTrueInuse_database_curvature, int mask)
{
    auto samples =
        sample_CURVEID_dimension_delta_s_angle_deviation_use_database_curvature<int, int, double, double, double>(
            curveidMin, curveidMax, curveidDistance,       // CURVEID: min, max, distance
            dimensionMin, dimensionMax, dimensionDistance, // dimension: min, max, distance 1
            delta_sMin, delta_sMax, delta_sDistance,       // delta_s: min, max, distance 2
            angleMin, angleMax, angleDistance,             // angle: min, max, distance 4
            deviationMin, deviationMax, deviationDistance, // deviation: min, max, distance 8
            includeTrueInuse_database_curvature, // includeTrueInBoolParam: true to include both false and true, false
                                                 // to include only false 16
            mask                                 // 用于选择参数的掩码
        );

    return samples;
}

bool test_line_mesh_gen2(TiGER::Geometry& Geo,
                         std::vector<discreteCurveDimensionParams>& discreteCurveDimensionParamsVector,
                         std::vector<LineMesh>& segements)
{

    std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> surfaces_(Geo.curves_.size());
    for( int i = 0; i < Geo.topo_surf_to_curves_.size(); i++ )
    {
        for( int j = 0; j < Geo.topo_surf_to_curves_[i].size(); j++ )
        {
            surfaces_[Geo.topo_surf_to_curves_[i][j]].push_back(Geo.surfaces_[i]);
        }
    }
    std::vector<std::shared_ptr<GeometryCurve>> singleCurves_(1);
    // singleCurves_[0].
    std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>> singleSurface(1);
    CurveParametersDimension dimension_args;
    for( int i = 0; i < discreteCurveDimensionParamsVector.size(); i++ )
    {
        singleCurves_[0] = Geo.curves_[discreteCurveDimensionParamsVector[i].CURVEID];
        singleSurface[0] = surfaces_[discreteCurveDimensionParamsVector[i].CURVEID];
        int mask = discreteCurveDimensionParamsVector[i].mask;
        if( mask & 1 )
        {
            dimension_args.setdimension(discreteCurveDimensionParamsVector[i].dimension);
        }
        mask >>= 1;
        if( mask & 1 )
        {
            dimension_args.setdelta_s(discreteCurveDimensionParamsVector[i].delta_s);
        }
        mask >>= 1;
        if( mask & 1 )
        {
            dimension_args.setangle(discreteCurveDimensionParamsVector[i].angle);
        }
        mask >>= 1;
        if( mask & 1 )
        {
            dimension_args.setdeviation(discreteCurveDimensionParamsVector[i].deviation);
        }
        mask >>= 1;
        if( mask & 1 )
        {
            dimension_args.setuse_database_curvature(discreteCurveDimensionParamsVector[i].use_database_curvature);
        }
        mask >>= 1;
        // std::cout << "this round" << " dimension: " << dimension_args.getdimension()
        //           << ", delta_s: " << dimension_args.getdelta_s() << ", angle: " << dimension_args.getangle()
        //           << ", deviation: " << dimension_args.getdeviation() << ", use_database_curvature: " <<
        //           std::boolalpha
        //           << dimension_args.getuse_database_curvature() << std::endl;
        TiGER::discretize_curve::Grid1D_Dimension(singleCurves_, singleSurface, dimension_args, segements);
    }

    return true;
}

struct discretizeSurfaceTriParams
{
    int SurfaceID;
    double mUserSpecifiedMaxEdgeLen;
    double mUserSpecifiedMinEdgeLen;
    double mMaxAngle;
    double mMaxDeviation;
    bool mSwapCellWithNoInteriorPoints;
    int mask;
};
#define SAMPLE_DISCRETIZESURFACETRI(value1, value2, value3, value4, value5, value6)                                    \
    template <typename T1, typename T2, typename T3, typename T4, typename T5>                                         \
    std::vector<discretizeSurfaceTriParams> sample_##value1##_##value2##_##value3##_##value4##_##value5##_##value6(    \
        T1 min1, T1 max1, T1 distance1, T2 min2, T2 max2, T2 distance2, T3 min3, T3 max3, T3 distance3, T4 min4,       \
        T4 max4, T4 distance4, T5 min5, T5 max5, T5 distance5, bool includeTrueInBoolParam, int mask)                  \
    {                                                                                                                  \
        std::vector<discretizeSurfaceTriParams> ans;                                                                   \
        for( T1 i = min1; i <= max1; i += distance1 )                                                                  \
        {                                                                                                              \
            for( T2 j = min2; j <= max2; j += distance2 )                                                              \
            {                                                                                                          \
                for( T3 k = min3; k <= max3; k += distance3 )                                                          \
                {                                                                                                      \
                    for( T4 l = min4; l <= max4; l += distance4 )                                                      \
                    {                                                                                                  \
                        for( T5 m = min5; m <= max5; m += distance5 )                                                  \
                        {                                                                                              \
                            if( includeTrueInBoolParam )                                                               \
                            {                                                                                          \
                                for( bool n : {false, true} )                                                          \
                                {                                                                                      \
                                    discretizeSurfaceTriParams a;                                                      \
                                    a.value1 = i;                                                                      \
                                    a.value2 = j;                                                                      \
                                    a.value3 = k;                                                                      \
                                    a.value4 = l;                                                                      \
                                    a.value5 = m;                                                                      \
                                    a.value6 = n;                                                                      \
                                    a.mask = mask;                                                                     \
                                    ans.push_back(a);                                                                  \
                                }                                                                                      \
                            }                                                                                          \
                            else                                                                                       \
                            {                                                                                          \
                                discretizeSurfaceTriParams a;                                                          \
                                a.value1 = i;                                                                          \
                                a.value2 = j;                                                                          \
                                a.value3 = k;                                                                          \
                                a.value4 = l;                                                                          \
                                a.value5 = m;                                                                          \
                                a.value6 = false;                                                                      \
                                a.mask = mask;                                                                         \
                                ans.push_back(a);                                                                      \
                            }                                                                                          \
                        }                                                                                              \
                    }                                                                                                  \
                }                                                                                                      \
            }                                                                                                          \
        }                                                                                                              \
        return ans;                                                                                                    \
    }
SAMPLE_DISCRETIZESURFACETRI(SurfaceID, mUserSpecifiedMaxEdgeLen, mUserSpecifiedMinEdgeLen, mMaxAngle, mMaxDeviation,
                            mSwapCellWithNoInteriorPoints)
std::vector<discretizeSurfaceTriParams> sampleDiscretizeSurfaceTri(
    int surfaceidMin, int surfaceidMax, int surfaceidDistance, int mUserSpecifiedMaxEdgeLenMin,
    int mUserSpecifiedMaxEdgeLenMax, int mUserSpecifiedMaxEdgeLenDistance, double mUserSpecifiedMinEdgeLenMin,
    double mUserSpecifiedMinEdgeLenMax, double mUserSpecifiedMinEdgeLenDistance, double mMaxAngleMin,
    double mMaxAngleMax, double mMaxAngleDistance, double mMaxDeviationMin, double mMaxDeviationMax,
    double mMaxDeviationDistance, bool includeTrueInmSwapCellWithNoInteriorPoints, int mask)
{
    auto samples =
        sample_SurfaceID_mUserSpecifiedMaxEdgeLen_mUserSpecifiedMinEdgeLen_mMaxAngle_mMaxDeviation_mSwapCellWithNoInteriorPoints<
            int, int, double, double, double>(
            surfaceidMin, surfaceidMax, surfaceidDistance, // surfaceid: min, max, distance
            mUserSpecifiedMaxEdgeLenMin, mUserSpecifiedMaxEdgeLenMax,
            mUserSpecifiedMaxEdgeLenDistance, // dimension: min, max, distance 1
            mUserSpecifiedMinEdgeLenMin, mUserSpecifiedMinEdgeLenMax,
            mUserSpecifiedMinEdgeLenDistance,                          // delta_s: min, max, distance 2
            mMaxAngleMin, mMaxAngleMax, mMaxAngleDistance,             // angle: min, max, distance 4
            mMaxDeviationMin, mMaxDeviationMax, mMaxDeviationDistance, // deviation: min, max, distance 8
            includeTrueInmSwapCellWithNoInteriorPoints, // includeTrueInBoolParam: true to include both false and true,
                                                        // false to include only false 16
            mask                                        // 用于选择参数的掩码
        );

    return samples;
}
bool test_surf_mesh_gen2(TiGER::Geometry& Geo,
                         std::vector<discretizeSurfaceTriParams>& discretizeSurfaceTriParamsVector,
                         std::vector<LineMesh>& linemesh, SurfaceMesh& surfacemesh)
{
    SizingFunction size_field = [](const double& x, const double& y, const double& z) { return 1; };
    std::vector<std::vector<LineMesh>> sf_boundary;
    sf_boundary.resize(Geo.surfaces_.size());
    std::cout << Geo.topo_surf_to_curves_.size() << std::endl;
    for( int i = 0; i < Geo.topo_surf_to_curves_.size(); i++ )
    {
        // std::cout << Geo.topo_surf_to_curves_[i].size() << std::endl;
        for( int j = 0; j < Geo.topo_surf_to_curves_[i].size(); j++ )
        {
            sf_boundary[i].push_back(linemesh[Geo.topo_surf_to_curves_[i][j]]);
        }
    }
    std::vector<SurfaceMesh> sf_mesh(sf_boundary.size());

    //////////////////////////////////////////////////////////////////////////////////////

    for( int i = 0; i < discretizeSurfaceTriParamsVector.size(); i++ )
    {
        IsotropicSurfaceParametersTri surface_args;
        int mask = discretizeSurfaceTriParamsVector[i].mask;
        std::cout << mask << std::endl;
        if( mask & 1 )
        {
            surface_args.setmUserSpecifiedMaxEdgeLen(discretizeSurfaceTriParamsVector[i].mUserSpecifiedMaxEdgeLen);
        }
        mask >>= 1;
        std::cout << mask << std::endl;
        if( mask & 1 )
        {
            surface_args.setmUserSpecifiedMinEdgeLen(discretizeSurfaceTriParamsVector[i].mUserSpecifiedMinEdgeLen);
        }
        mask >>= 1;
        std::cout << mask << std::endl;
        if( mask & 1 )
        {
            surface_args.setmMaxAngle(discretizeSurfaceTriParamsVector[i].mMaxAngle);
        }
        mask >>= 1;
        std::cout << mask << std::endl;
        if( mask & 1 )
        {
            surface_args.setmMaxDeviation(discretizeSurfaceTriParamsVector[i].mMaxDeviation);
        }
        std::cout << mask << std::endl;
        mask >>= 1;
        if( mask & 1 )
        {
            surface_args.setmSwapCellWithNoInteriorPoints(
                discretizeSurfaceTriParamsVector[i].mSwapCellWithNoInteriorPoints);
        }
        mask >>= 1;
        std::cout << mask << std::endl;

        TiGER::discretize_surface::CADSurf_TGrid_AFT(Geo.surfaces_[discretizeSurfaceTriParamsVector[i].SurfaceID],
                                                        sf_boundary[discretizeSurfaceTriParamsVector[i].SurfaceID],
                                                        surface_args, size_field,
                                                        sf_mesh[discretizeSurfaceTriParamsVector[i].SurfaceID]);
    }
    // SurfaceMesh sf_output;
    // TiGER::merge_surface_mesh(sf_mesh, sf_output);
    // surfacemesh = sf_output;
    return true;
}
struct discreteCurveTanhParams
{
    double begin_s;
    double end_s;
    int dimension;
    int mask;

    discreteCurveTanhParams(double begin_s_val = 0, double end_s_val = 0, int dimension_val = 20)
        : begin_s(begin_s_val)
        , end_s(end_s_val)
        , dimension(dimension_val)
    {
    }
};
#define SAMPLE_TANH(value1, value2, value3)                                                                            \
    template <typename T1, typename T2, typename T3>                                                                   \
    std::vector<discreteCurveTanhParams> sample_##value1##_##value2##_##value3(                                        \
        T1 min1, T1 max1, T1 distance1, T2 min2, T2 max2, T2 distance2, T3 min3, T3 max3, T3 distance3, int mask)      \
    {                                                                                                                  \
        std::vector<discreteCurveTanhParams> ans;                                                                      \
        for( T1 i = min1; i <= max1; i += distance1 )                                                                  \
        {                                                                                                              \
            for( T2 j = min2; j <= max2; j += distance2 )                                                              \
            {                                                                                                          \
                for( T3 k = min3; k <= max3; k += distance3 )                                                          \
                {                                                                                                      \
                    discreteCurveTanhParams a;                                                                         \
                    a.value1 = i;                                                                                      \
                    a.value2 = j;                                                                                      \
                    a.value3 = k;                                                                                      \
                    a.mask = mask;                                                                                     \
                    ans.push_back(a);                                                                                  \
                }                                                                                                      \
            }                                                                                                          \
        }                                                                                                              \
        return ans;                                                                                                    \
    }
SAMPLE_TANH(begin_s, end_s, dimension)
std::vector<discreteCurveTanhParams> sampleDiscreteCurveTanh(double begin_sMin, double begin_sMax,
                                                             double begin_sDistance, double end_sMin, double end_sMax,
                                                             double end_sDistance, int dimensionMin, int dimensionMax,
                                                             int dimensionDistance, int mask)
{
    auto samples = sample_begin_s_end_s_dimension<double, double, int>(
        begin_sMin, begin_sMax, begin_sDistance,       // begin_s: min, max, distance 1
        end_sMin, end_sMax, end_sDistance,             // end_s: min, max, distance 2
        dimensionMin, dimensionMax, dimensionDistance, // dimension: min, max, distance 4
        mask                                           // 用于选择参数的掩码
    );

    return samples;
}

bool test_line_mesh_gen3(TiGER::Geometry& Geo, std::vector<discreteCurveTanhParams>& discreteCurveTanhParamsVector,
                         std::vector<LineMesh>& segements)
{
    CurveParametersTanh Tanh_args;
    for( int i = 0; i < discreteCurveTanhParamsVector.size(); i++ )
    {
        int mask = discreteCurveTanhParamsVector[i].mask;
        std::cout << mask << std::endl;
        if( mask & 1 )
        {
            Tanh_args.setbegin_s(discreteCurveTanhParamsVector[i].begin_s);
        }
        mask >>= 1;
        std::cout << mask << std::endl;
        if( mask & 1 )
        {
            Tanh_args.setend_s(discreteCurveTanhParamsVector[i].end_s);
        }
        mask >>= 1;
        std::cout << mask << std::endl;
        if( mask & 1 )
        {
            Tanh_args.setdimension(discreteCurveTanhParamsVector[i].dimension);
        }
        mask >>= 1;
        // std::cout << mask << std::endl;
        // std::cout << "this round" << " begin_s: " << Tanh_args.getbegin_s() << ", end_s: " << Tanh_args.getend_s()
        //           << ", dimension: " << Tanh_args.getdimension() << std::endl;
        TiGER::discretize_curve::Grid1D_Tanh(Geo.curves_, Tanh_args, segements);
    }

    return true;
}
// std::vector<double> sampleParameters(double min, double max, double sampleDistance)
// {
//     std::vector<double> samples;

//     // 确保最大值大于最小值，并且步长为正数
//     if( max > min && sampleDistance > 0 )
//     {
//         for( double value = min; value <= max; value += sampleDistance )
//         {
//             samples.push_back(value);
//         }
//     }

//     return samples;
// }
// bool checkMeshQuality(SurfaceMesh& surfacemesh, vector<double> value)
bool checkMeshQuality(SurfaceMesh& surfacemesh)
{
    const TiGER::SurfMeshQualityType type = TiGER::SurfMeshQualityType::TRI_ANGLE;
    std::vector<TiGER::QualityResult> quality;
    TiGER::Tool_SurfQuality(surfacemesh, type, quality);

    printf("quality size: %ld\n", quality.size());
    std::cout << quality[0].ave_value << std::endl;
    std::cout << quality[0].min_value << std::endl;
    std::cout << quality[0].max_value << std::endl;
    std::cout << quality[0].min_index << std::endl;
    std::cout << quality[0].max_index << std::endl;
    std::cout << quality[0].min_index_2 << std::endl;
    std::cout << quality[0].max_index_2 << std::endl;
    // for( auto i : quality )
    // {
    //     if( i.ave_value < value[0] )
    //     {
    //         return false;
    //     }
    //     if( i.max_value < value[1] )
    //     {
    //         return false;
    //     }
    //     if( i.min_value < value[2] )
    //     {
    //         return false;
    //     }
    // }
    return true;
}
// 包围盒对角线长度
double getDiagonalLength(std::vector<LineMesh>& segements)
{
    double minx = DBL_MAX, miny = DBL_MAX, minz = DBL_MAX;
    double maxx = -DBL_MAX, maxy = -DBL_MAX, maxz = -DBL_MAX;
    std::cout << "nihao3";
    for( int j = 0; j < segements.size(); j++ )
    {
        for( int k = 0; k < segements[j].coord.size(); k++ )
        {
            if( segements[j].coord[k][0] > maxx )
            {
                maxx = segements[j].coord[k][0];
            }
            if( segements[j].coord[k][1] > maxy )
            {
                maxy = segements[j].coord[k][1];
            }
            if( segements[j].coord[k][2] > maxz )
            {
                maxz = segements[j].coord[k][2];
            }

            if( segements[j].coord[k][0] < minx )
            {
                minx = segements[j].coord[k][0];
            }
            if( segements[j].coord[k][1] < miny )
            {
                miny = segements[j].coord[k][1];
            }
            if( segements[j].coord[k][2] < minz )
            {
                minz = segements[j].coord[k][2];
            }
        }
    }
    std::cout << minx << ' ' << miny << ' ' << minz << std::endl
              << maxx << ' ' << maxy << ' ' << maxz << std::endl
              << std::endl;
    double distance =
        std::sqrt((maxx - minx) * (maxx - minx) + (maxy - miny) * (maxy - miny) + (maxz - minz) * (maxz - minz));
    std::cout << "distance" << distance << std::endl;
    return distance;
}

bool checkMeshQuality(SurfaceMesh& surfacemesh, vector<double> value)
{
    bool quality_ret = true;

    const TiGER::SurfMeshQualityType skewness_type = TiGER::SurfMeshQualityType::TRI_EQUIANGLE_SKEWNESS;
    std::vector<TiGER::QualityResult> skewness_quality;
    TiGER::Tool_SurfQuality(surfacemesh, skewness_type, skewness_quality);

    printf("skewness quality size: %d\n", skewness_quality.size());

    for( auto i : skewness_quality )
    {
        std::cout << "skewness max: " << i.max_value;
        if( i.max_value < value[0] )
        {
            std::cout << " < ";
        }
        else
        {
            std::cout << " > ";
            quality_ret = false;
        }
        std::cout << value[0] << " ave_value :" << i.ave_value;
        if( i.ave_value < value[1] )
        {
            std::cout << " < ";
        }
        else
        {
            std::cout << " > ";
            quality_ret = false;
        }
        std::cout << value[1] << std::endl;
    }

    const TiGER::SurfMeshQualityType aspect_type = TiGER::SurfMeshQualityType::TRI_ASPECT;
    std::vector<TiGER::QualityResult> aspect_quality;
    TiGER::Tool_SurfQuality(surfacemesh, aspect_type, aspect_quality);

    printf("aspect quality size: %d\n", aspect_quality.size());

    for( auto i : aspect_quality )
    {
        std::cout << "aspect max: " << i.max_value;
        if( i.max_value < value[2] )
        {
            std::cout << " < ";
        }
        else
        {
            std::cout << " > ";
            quality_ret = false;
        }
        std::cout << value[2] << " ave_value :" << i.ave_value;
        if( i.ave_value < value[3] )
        {
            std::cout << " < ";
        }
        else
        {
            std::cout << " > ";
            quality_ret = false;
        }
        std::cout << value[3] << std::endl;
    }

    const TiGER::SurfMeshQualityType angle_type = TiGER::SurfMeshQualityType::TRI_ANGLE;
    std::vector<TiGER::QualityResult> angle_quality;
    TiGER::Tool_SurfQuality(surfacemesh, angle_type, angle_quality);

    printf("angle quality size: %d\n", angle_quality.size());

    for( auto i : angle_quality )
    {
        std::cout << "angle min: " << i.min_value;
        if( i.min_value > value[4] )
        {
            std::cout << " > ";
        }
        else
        {
            std::cout << " < ";
            quality_ret = false;
        }
        std::cout << value[4] << " max :" << i.max_value;
        if( i.max_value < value[5] )
        {
            std::cout << " < ";
        }
        else
        {
            std::cout << " > ";
            quality_ret = false;
        }
        std::cout << value[5] << std::endl;
    }

    const TiGER::SurfMeshQualityType condition_type = TiGER::SurfMeshQualityType::TRI_CONDITION;
    std::vector<TiGER::QualityResult> condition_quality;
    TiGER::Tool_SurfQuality(surfacemesh, condition_type, condition_quality);

    printf("condition quality size: %d\n", condition_quality.size());

    for( auto i : condition_quality )
    {
        std::cout << "condition max: " << i.max_value;
        if( i.max_value < value[6] )
        {
            std::cout << " < ";
        }
        else
        {
            std::cout << " > ";
            quality_ret = false;
        }
        std::cout << value[6] << " ave_value :" << i.ave_value;
        if( i.ave_value < value[7] )
        {
            std::cout << " < ";
        }
        else
        {
            std::cout << " > ";
            quality_ret = false;
        }
        std::cout << value[7] << std::endl;
    }

    return quality_ret;
    // return true;
}

// bool isLineMeshEmpty(LineMesh& linemesh)
// {
//     if( linemesh.coord.size() == 0 )
//     {
//         return true;
//     }
//     return false;
// }
// bool isSurfaceMeshEmpty(SurfaceMesh& surfacemesh)
// {
//     if( surfacemesh.coords.size() == 0 )
//     {
//         return true;
//     }
//     return false;
// }

#define SAMPLE1(value1, value2)                                                                                        \
    template <typename T, typename MemberType, MemberType T::*MemberPtr>                                               \
    std::vector<T> sample_##value1##_##value2(MemberType min, MemberType max, MemberType distance, int mask,           \
                                              bool judge = false)                                                      \
    {                                                                                                                  \
        std::vector<T> ans;                                                                                            \
        for( MemberType i = min; i <= max; i += distance )                                                             \
        {                                                                                                              \
            T a = T();                                                                                                 \
            a.*MemberPtr = i;                                                                                          \
            a.mask = mask;                                                                                             \
            ans.push_back(a);                                                                                          \
            if( judge )                                                                                                \
            {                                                                                                          \
                T b = T(true);                                                                                         \
                b.*MemberPtr = i;                                                                                      \
                b.mask = mask;                                                                                         \
                ans.push_back(b);                                                                                      \
            }                                                                                                          \
        }                                                                                                              \
        return ans;                                                                                                    \
    }

#define SAMPLE2(value1, value2, value3)                                                                                \
    template <typename T, typename MemberType1, MemberType1 T::*MemberPtr1, typename MemberType2,                      \
              MemberType2 T::*MemberPtr2>                                                                              \
    std::vector<T> sample_##value1##_##value2##_##value3(MemberType1 min1, MemberType1 max1, MemberType1 distance1,    \
                                                         MemberType2 min2, MemberType2 max2, MemberType2 distance2,    \
                                                         int mask, bool judge = false)                                 \
    {                                                                                                                  \
        std::vector<T> ans;                                                                                            \
        for( MemberType1 i = min1; i <= max1; i += distance1 )                                                         \
        {                                                                                                              \
            for( MemberType2 j = min2; j <= max2; j += distance2 )                                                     \
            {                                                                                                          \
                T a = T();                                                                                             \
                a.*MemberPtr1 = i;                                                                                     \
                a.*MemberPtr2 = j;                                                                                     \
                a.mask = mask;                                                                                         \
                ans.push_back(a);                                                                                      \
                if( judge )                                                                                            \
                {                                                                                                      \
                    T b = T(true);                                                                                     \
                    b.*MemberPtr1 = i;                                                                                 \
                    b.*MemberPtr2 = j;                                                                                 \
                    b.mask = mask;                                                                                     \
                    ans.push_back(b);                                                                                  \
                }                                                                                                      \
            }                                                                                                          \
        }                                                                                                              \
        return ans;                                                                                                    \
    }

#define SAMPLE3(value1, value2, value3, value4)                                                                        \
    template <typename T, typename MemberType1, MemberType1 T::*MemberPtr1, typename MemberType2,                      \
              MemberType2 T::*MemberPtr2, typename MemberType3, MemberType3 T::*MemberPtr3>                            \
    std::vector<T> sample_##value1##_##value2##_##value3##_##value4(                                                   \
        MemberType1 min1, MemberType1 max1, MemberType1 distance1, MemberType2 min2, MemberType2 max2,                 \
        MemberType2 distance2, MemberType3 min3, MemberType3 max3, MemberType3 distance3, int mask,                    \
        bool judge = false)                                                                                            \
    {                                                                                                                  \
        std::vector<T> ans;                                                                                            \
        for( MemberType1 i = min1; i <= max1; i += distance1 )                                                         \
        {                                                                                                              \
            for( MemberType2 j = min2; j <= max2; j += distance2 )                                                     \
            {                                                                                                          \
                for( MemberType3 k = min3; k <= max3; k += distance3 )                                                 \
                {                                                                                                      \
                    T a = T();                                                                                         \
                    a.*MemberPtr1 = i;                                                                                 \
                    a.*MemberPtr2 = j;                                                                                 \
                    a.*MemberPtr3 = k;                                                                                 \
                    a.mask = mask;                                                                                     \
                    ans.push_back(a);                                                                                  \
                    if( judge )                                                                                        \
                    {                                                                                                  \
                        T b = T(true);                                                                                 \
                        b.*MemberPtr1 = i;                                                                             \
                        b.*MemberPtr2 = j;                                                                             \
                        b.*MemberPtr3 = k;                                                                             \
                        b.mask = mask;                                                                                 \
                        ans.push_back(b);                                                                              \
                    }                                                                                                  \
                }                                                                                                      \
            }                                                                                                          \
        }                                                                                                              \
        return ans;                                                                                                    \
    }

SAMPLE2(discreteCurveDimensionParams, CURVEID, delta_s)
SAMPLE3(discreteCurveDimensionParams, CURVEID, dimension, angle)
SAMPLE2(discreteCurveDimensionParams, CURVEID, dimension)
SAMPLE3(discreteCurveDimensionParams, CURVEID, dimension, deviation)
SAMPLE2(discreteCurveTanhParams, end_s, dimension)
SAMPLE2(discreteCurveTanhParams, begin_s, dimension)
SAMPLE2(discreteCurveTanhParams, begin_s, end_s)