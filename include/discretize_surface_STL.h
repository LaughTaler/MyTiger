#include "geometry_data_struct.h"
#include "./sizing_field_interface.h"
#include "discretize_surface_parameters.h"
#include "./mesh_data_struct.h"
#include "./definitions.h"
#include "alias.h"
#include <memory>
#include <vector>

/**
 * @brief 曲面离散函数集 
 * 
 */
namespace TiGER{
namespace discretize_surface{


  /**
     * @brief 离散单张曲面生成三角形网格，利用四叉树加密得到初始采样点，再使用DT2D，得到面网格，最后用remesh精炼
     * 
     * @param[in] surface 几何曲面
     * @param[in] boundary_segmenets 边界曲线
     * @param[in] sizefield 尺寸场
     * @param[out] output_mesh 输出的网格, 
     * 拓扑存储在output_mesh.tris，
     * 点存储在output_mesh.coords,
     * 以及output_mesh.uv，
     */
Tiger_API void CADSurf_exportToSTL(
    std::shared_ptr<GeometrySurface>& surface, const std::vector<LineMesh>& boundary_segments,
    const IsotropicSurfaceParametersTri& iso_args, const SizingFunction& sizefield, SurfaceMesh& output_mesh);

} // discretize_curve
} // TiGER