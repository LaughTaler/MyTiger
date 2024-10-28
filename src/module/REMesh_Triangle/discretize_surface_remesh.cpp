#include "discretize_surface_remesh.h"
#include "geometry_data_struct.h"
#include "api_remesh.h"
#include "alias.h"
#include "discretize_curve_parameters.h"
#include <iostream>

namespace TiGER {
namespace discretize_surface_remesh {

    void TMeshSurf_Opt(
    const SurfaceMesh &surfmesh_in,
    const std::vector<std::vector<int>> &constrained_edge,
    const std::vector<std::vector<int>> &conformaing_edge,
    const RemeshParameters &args,
    std::function<double(const double&, const double&, const double&)> &size_function,
    SurfaceMesh &surfmesh_out,
    std::vector<std::vector<int>> &constrained_edge_out) {
  std::string str_args = args.toCommandLineString();
  std::cout << str_args << "\n";
  std::vector<std::array<double, 3>> vertexs_in = surfmesh_in.coords;
  std::vector<std::array<int, 3>> facets_in = surfmesh_in.tris;
  std::vector<int> facets_mark = surfmesh_in.attribute_int;

  std::vector<std::array<double, 3>> vertexs_out;
  std::vector<std::array<int, 3>> facets_out;
  std::vector<int> marks_out;

  std::function<double(double, double, double)> size_function2 =
      [&](double x, double y, double z) { return size_function(x, y, z); };
  TriangleRemesh(vertexs_in, facets_in, facets_mark, constrained_edge,
                 conformaing_edge, str_args, size_function2, vertexs_out,
                 facets_out, marks_out, constrained_edge_out);
  surfmesh_out.coords = vertexs_out;
  surfmesh_out.tris = facets_out;
  surfmesh_out.attribute_int = marks_out;
}


void TMeshSurf_Opt(
    const SurfaceMesh &surfmesh_in,
    const std::vector<std::vector<int>> &constrained_edge,
    const std::vector<std::vector<int>> &conformaing_edge,
    const RemeshParameters &args,
    std::function<double(double, double, double)> &size_function,
    SurfaceMesh &surfmesh_out,
    std::vector<std::vector<int>> &constrained_edge_out) {
  std::string str_args = args.toCommandLineString();
  std::cout << str_args << "\n";
  std::vector<std::array<double, 3>> vertexs_in = surfmesh_in.coords;
  std::vector<std::array<int, 3>> facets_in = surfmesh_in.tris;
  std::vector<int> facets_mark = surfmesh_in.attribute_int;

  std::vector<std::array<double, 3>> vertexs_out;
  std::vector<std::array<int, 3>> facets_out;
  std::vector<int> marks_out;
  TriangleRemesh(vertexs_in, facets_in, facets_mark, constrained_edge,
                 conformaing_edge, str_args, size_function, vertexs_out,
                 facets_out, marks_out, constrained_edge_out);
  surfmesh_out.coords = vertexs_out;
  surfmesh_out.tris = facets_out;
  surfmesh_out.attribute_int = marks_out;
}

void TMeshSurf_Opt(
    const SurfaceMesh &surfmesh_in,
    const std::vector<GeometrySurface> &surfaces_in,
    const std::vector<std::vector<int>> &constrained_edge,
    const std::vector<std::vector<int>> &conforming_edge,
    const RemeshParameters &args,
    std::function<double(double, double, double)> &size_function,
    SurfaceMesh &surfmesh_out,
    std::vector<std::vector<int>> &constrained_edge_out) {
  std::string str_args = args.toCommandLineString();
  std::vector<std::array<double, 3>> vertexs_in = surfmesh_in.coords;
  std::vector<std::array<int, 3>> facets_in = surfmesh_in.tris;
  std::vector<int> facets_mark = surfmesh_in.attribute_int;

  std::vector<std::array<double, 3>> vertexs_out;
  std::vector<std::array<int, 3>> facets_out;
  std::vector<int> marks_out;
  TriangleRemesh(vertexs_in, facets_in, facets_mark, constrained_edge,
                 conforming_edge, str_args, size_function, vertexs_out,
                 facets_out, marks_out, constrained_edge_out);
  std::vector<bool> vis(vertexs_out.size(), false);
  for(int i = 0; i < facets_out.size(); ++i){
    for(int j = 0; j < 3; ++j){
      int v = facets_out[i][j];
      if(vis[v]) continue;
      vis[v] = 1;
      std::array<double, 2> res_uv;
      std::array<double, 4> uv_scale = surfaces_in[marks_out[i]].getUVScaleFunction();
      std::array<double, 2> hint_uv = {(uv_scale[0] + uv_scale[1]) * 0.5, (uv_scale[2] + uv_scale[3]) * 0.5};
      double eps=1e-3;
      surfaces_in[marks_out[i]].projectFunction(vertexs_out[v], eps, hint_uv, res_uv);
      vertexs_out[v] = surfaces_in[marks_out[i]].d0Function(res_uv);
    }
  }
  surfmesh_out.coords = vertexs_out;
  surfmesh_out.tris = facets_out;
  surfmesh_out.attribute_int = marks_out;
}

}  // namespace discretize_surface_remesh
}  // namespace TiGER