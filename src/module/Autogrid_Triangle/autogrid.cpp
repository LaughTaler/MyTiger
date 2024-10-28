#include <array>
#include <vector>

#include "mesh_data_struct.h"
#include "autogrid.h"
#include "alias.h"
#include "Autogrid_base.h"
#include "list_to_matrix.hpp"
#include "matrix_to_list.hpp"



/**
 * @brief 修复脏离散三角形网格，得到水密的无自相交的面网格
 *
 * @param[in] surf_in 输入三角形表面网格
 * @param[in] local_eps 输入局部容差
 * @param[in] args 输入字符串控制参数
 * @param[out] surf_out 输出三角形表面网格
 */
void TiGER::discrete_geometry_repair::TMeshSurf_Repair(const TiGER::SurfaceMesh& surf_in,
                         const TiGER::AutogridParameters& args,
                         TiGER::SurfaceMesh& surf_out, ConstrainDelaunayFunc constrain_delaunay) {
  autogrid::MeshT mesh, mesh_out;

  TiGER::list_to_matrix<3>(surf_in.coords, mesh.V);
  TiGER::list_to_matrix<3>(surf_in.tris, mesh.T);
  if (surf_in.edges.size() > 0) TiGER::list_to_matrix<2>(surf_in.edges, mesh.E);
  if(surf_in.point_attribute_double.size() == surf_in.coords.size()){
    TiGER::list_to_matrix<1>(surf_in.point_attribute_double, mesh.local_eps);
  }
  if (surf_in.attribute_int.size() > 0)
    TiGER::list_to_matrix<1>(surf_in.attribute_int, mesh.M);
  std::string str_args = args.toCommandLineString();
  autogrid::autogrid(mesh, str_args, mesh_out, constrain_delaunay);
  TiGER::matrix_to_list<3>(mesh_out.V, surf_out.coords);
  TiGER::matrix_to_list<3>(mesh_out.T, surf_out.tris);
  TiGER::matrix_to_list<2>(mesh_out.E, surf_out.edges);
  if (mesh_out.M.rows() > 0) {
    for (int i = 0; i < mesh_out.M.rows(); ++i) {
      surf_out.attribute_int.push_back(mesh_out.M(i, 0));
    }
  }
}
