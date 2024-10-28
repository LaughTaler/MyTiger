#include <array>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "./definitions.h"
#include "./mesh_data_struct.h"
#include "./sizing_field_interface.h"
#include "./tetrahedral_mesh_parameters.h"
namespace TiGER
{
enum AttributeType
{
    point_id = 0,
    sizing_value,
    surface_id,
    color,
    point_normal,
    facet_normal
};

/**
 * @brief 体网格生成模块命名空间
 */
namespace tetrahedral_mesh
{
/**
 * @brief 约束体网格边界恢复
 *
 * @param surface_mesh 输入的表面网格网格拓扑
 * @param args 输入的控制参数
 * @param volume_mesh 输出的体网格拓扑
 */
Tiger_API void constrainDelaunay(const SurfaceMesh& surface_mesh, const TetrahedraParameters& args,
                                 VolumeMesh& volume_mesh, const SizingFunction& func);

/**
 * @brief 约束体网格边界恢复
 *
 * @param surface_mesh 输入的表面网格网格拓扑
 * @param args 输入的控制参数
 * @param volume_mesh 输出的体网格拓扑
 */
Tiger_API void TetGrid_Delaunay_QuadInputs(const SurfaceMesh& surface_mesh, const TetrahedraParameters& args,
                                           VolumeMesh& volume_mesh, const SizingFunction& func);

/**
 * @brief 各向异性的体网格优化
 *
 * @param volume_mesh_in 输入的体网格拓扑
 * @param args 输入的控制参数
 * @param volume_mesh_out 输出的体网格拓扑
 * @param func 输入的张量尺寸场
 */
Tiger_API void TetGrid_Opt(const VolumeMesh& volume_mesh_in, const TetrahedraParameters& args,
                           VolumeMesh& volume_mesh_out,
                           std::function<std::array<double, 9>(std::array<double, 3>)>& func);
} // namespace tetrahedral_mesh
} // namespace TiGER
