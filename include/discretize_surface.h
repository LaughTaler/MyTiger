#include "geometry_data_struct.h"
#include "./sizing_field_interface.h"
#include "discretize_surface_parameters.h"
#include "./mesh_data_struct.h"
#include <memory>
#include <vector>
#include "alias.h"
/**
 * @brief 曲面离散函数集 
 * 
 */
namespace TiGER{
namespace discretize_surface{


    /**
     * @brief 离散单张曲面生成三角形网格
     * 
     * @param[in] surface 几何曲面
     * @param[in] boundary_segmenets 边界曲线
     * @param[in] sizefield 尺寸场
     * @param[out] output_mesh 输出的网格, 
     * 拓扑存储在output_mesh.tris，
     * 点存储在output_mesh.coords,
     * 以及output_mesh.uv，
     */

    Tiger_API void CADSurf_TGrid_AFT(

        std::shared_ptr<GeometrySurface>& surface,
        const std::vector<LineMesh>& boundary_segments,
        const IsotropicSurfaceParametersTri &iso_args,
        const SizingFunction& sizefield,
        SurfaceMesh& output_mesh
    );


        /**
     * @brief 离散单张曲面生成三角形网格
     * 
     * @param[in] surface 几何曲面
     * @param[in] boundary_segmenets 边界曲线
     * @param[in] sizefield 尺寸场
     * @param[in] bcs 边界条件
     * @param[out] output_mesh 输出的网格, 
     * 拓扑存储在output_mesh.tris，
     * 点存储在output_mesh.coords,
     * 以及output_mesh.uv，
     */
    Tiger_API void CADSurf_TGrid_AFLR(
        std::shared_ptr<GeometrySurface>& surface,
        const std::vector<std::shared_ptr<TiGER::GeometryCurve>>& curves,
        const AnisotropicSurfaceParametersTri& aniso_args,
        std::vector<LineMesh>& boundary_segments,
        const SizingFunction& sizefield,
        SurfaceMesh& output_mesh
    );

        /**
     * @brief 离散单张曲面生成四边形网格
     * 
     * @param[in] surface 几何曲面
     * @param[in] boundary_segmenets 边界曲线
     * @param[in] sizefield 尺寸场
     * @param[in] bcs 边界条件
     * @param[out] output_mesh 输出的网格, 
     * 拓扑存储在output_mesh.tris，
     * 点存储在output_mesh.coords,
     * 以及output_mesh.uv，
     */
    Tiger_API void CADSurf_QGrid_AFLR(
        const GeometrySurface& surface,
        const std::vector<LineMesh>& boundary_segments,
        const std::vector<SurfaceBoundaryCondition> bcs,
        const SizingFunction& sizefield,
        SurfaceMesh& output_mesh
    );
    /**
     * @brief 离散所有曲面生成三角形网格
     * 
     * @param[in] surface 几何曲面
     * @param[in] boundary_segments 裁剪曲线
     * @param[in] sizefield 尺寸场
     * @param[out] output_mesh 输出的网格, 
     * 拓扑存储在output_mesh.tris，
     * 点存储在output_mesh.coords,
     * 以及output_mesh.uv，
     * uv值面id存储在output_mesh.uv_region 
     */
    Tiger_API void VCADSurf_TGrid_AFT(
        const std::vector<std::shared_ptr<GeometrySurface>>& surfaces,
        const std::vector< std::vector< LineMesh> >& boundary_segments,
        const IsotropicSurfaceParametersTri& iso_args,
        const SizingFunction& sizefield,
        SurfaceMesh& output_mesh
    );
    
    /**
     * @brief 离散单张虚曲面生成三角形网格
     * 
     * @param[in] surfaces 一组几何曲面
     * @param[in] curves 一组几何曲线
     * @param[in] topo_surface_to_curves 曲面是由哪些曲线构成
     * @param[in] virtual_trimmed_curve 虚曲面的裁剪曲线
     * @param[in] sizefield 尺寸场
     * @param[out] output_mesh 输出的网格, 
     * 拓扑存储在output_mesh.tris，
     * 点存储在output_mesh.coords,
     * 以及output_mesh.uv，
     * uv值面id存储在output_mesh.uv_region 
     */
    Tiger_API void VSurf_TMesh_meshTriangulately(
        const std::vector<GeometrySurface>& surfaces,
        const std::vector<GeometryCurve>& curves,
        const std::vector<std::vector<int>>& topo_surface_to_curves,
        const std::vector<LineMesh> virtual_trimmed_curve,
        const SizingFunction& sizefield,
        SurfaceMesh& output_mesh
    );

    
    /**
     * @brief 离散单张曲面生成四边形为主网格
     * 
     * @param[in] surface 几何曲面
     * @param[in] trimmed_curve 裁剪曲线
     * @param[in] sizefield 尺寸场
     * @param[out] output_mesh 输出的网格, 拓扑存储在output_mesh.tris 和output_mesh.quads，点存储在output_mesh.coords 
     */
    Tiger_API void CADSurf_QGrid_Paving(
        const GeometrySurface& surface,
        const std::vector<LineMesh>& trimmed_curve,
        const SizingFunction& sizefield,
        SurfaceMesh& output_mesh
    );

    /**
     * @brief 离散所有曲面生成四边形为主网格
     * 
     * @param[in] surface 几何曲面
     * @param[in] trimmed_curve 裁剪曲线
     * @param[in] sizefield 尺寸场
     * @param[out] output_mesh 输出的网格, 拓扑存储在output_mesh.tris 和output_mesh.quads，点存储在output_mesh.coords 
     */
    Tiger_API void discretizeSurfacesQuad(
        const std::vector<GeometrySurface>& surfaces,
        const std::vector<LineMesh>& trimmed_curve,
        const SizingFunction& sizefield,
        SurfaceMesh& output_mesh
    );
    
    /**
     * @brief 离散单张虚曲面生成四面形为主网格
     * 
     * @param[in] surfaces 一组几何曲面
     * @param[in] curves 一组几何曲线
     * @param[in] topo_surface_to_curves 曲面与曲线的对应关系
     * @param[in] virtual_trimmed_curve 虚曲面的裁剪曲线
     * @param[in] sizefield 尺寸场
     * @param[out] output_mesh 输出的网格, 拓扑存储在output_mesh.tris 和output_mesh.quads，点存储在output_mesh.coords 
     */
    Tiger_API void VSurf_QMesh_Paving(
        const std::vector<GeometrySurface>& surfaces,
        const std::vector<GeometryCurve>& curves,
        const std::vector<std::vector<int>>& topo_surface_to_curves,
        const std::vector<LineMesh> virtual_trimmed_curve,
        const SizingFunction& sizefield,
        SurfaceMesh& output_mesh
    );


    /**
     * @brief 单张曲面生成结构四边形网格 
     *  
     *    __ __ __ __ __ __
     *   |     a           |
     * d |                 | b
     *   |                 |
     *   |__ __ __ __ __ __|
     *        c
     * @param[in] surface 支撑曲面
     * @param[in] discrete_boundary 离散边界曲线，共4条，如上图所示，非相邻边的拓扑需一致
     * 
     * @param[out] output_mesh 输出的网格, 拓扑存储output_mesh.quads，点存储在output_mesh.coords 
     */
    Tiger_API void CADSurf_QGrid_Submap(
        const GeometrySurface& surface,
        const std::array<LineMesh, 4>& discrete_boundary,
        SurfaceMesh& output_mesh);


    /**
     * @brief 离散虚曲面生成结构四边形网格 
     *  
     *    __ __ __ __ __ __
     *   |     a           |
     * d |                 | b
     *   |                 |
     *   |__ __ __ __ __ __|
     *        c
     * @param[in] surfaces 支撑曲面
     * @param[in] curves 一组几何曲线
     * @param[in] topo_surface_to_curves 曲面与曲线的对应关系
     * @param[in] discrete_boundary 离散边界曲线，共4条，如上图所示，非相邻边的拓扑需一致
     * 
     * @param[out] output_mesh 输出的网格, 拓扑存储output_mesh.quads，点存储在output_mesh.coords 
     */
   Tiger_API void VSurf_QGrid_Submap(
        const std::vector<GeometrySurface>& surfaces,
        const std::vector<GeometryCurve>& curves,
        const std::vector<std::vector<int>>& topo_surface_to_curves,
        const std::array<LineMesh,4>& discrete_boundary,
        SurfaceMesh& output_mesh
    );
} // discretize_curve
} // TiGER