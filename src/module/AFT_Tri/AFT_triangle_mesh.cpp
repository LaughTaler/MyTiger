#include "../tools/merge_surface_mesh.hpp"
#include "../tools/output_info.hpp"
#include "./meshIO.h"
#include "Args.h"
#include "discretize_curve.h"
#include "discretize_surface.h"
#include "discretize_surface_remesh.h"
#include "discretize_surface_remesh_parameters.h"
#include "igl/remove_duplicate_vertices.h"
#include "n_api.h"
#include "../include/alias.h"

namespace TiGER {
namespace discretize_curve {

void Grid1D_Legacy(const std::shared_ptr<TiGER::GeometryCurve> curve,
                   const SizingFunction& sizefield,
                   const TiGER::CurveParameters& args,
                   TiGER::LineMesh& segments) {
  SmeshGen::LineData line;
  SmeshGen::DiscreteSingleCurve(curve, sizefield, line);
  for (int j = 0; j < line.coord.size(); j++) {
    segments.coord.push_back(line.coord[j]);
  }
 /* for (int j = 0; j < line.point_normal.size(); j++) {
    segments.point_normal.push_back(line.point_normal[j]);
  }
  for (int j = 0; j < line.regions.size(); j++) {
    segments.regions.push_back(line.regions[j]);
  }*/
  return;
}

void discretePartCurve(const std::shared_ptr<TiGER::GeometryCurve> curve, const std::array<double, 3> StartPt,
                       const std::array<double, 3> EndPt, const SizingFunction& sizefield,
                       const TiGER::CurveParameters& args, TiGER::LineMesh& segments)
{
    SmeshGen::LineData line;
    segments.coord.clear();
    segments.point_normal.clear();
    segments.regions.clear();

    SmeshGen::DiscreteSinglePartCurve(curve, StartPt, EndPt, sizefield, line);
    
    for( int j = 0; j < line.coord.size(); j++ )
    {
        segments.coord.push_back(line.coord[j]);
    }
    for( int j = 0; j < line.point_normal.size(); j++ )
    {
        segments.point_normal.push_back(line.point_normal[j]);
    }
    for( int j = 0; j < line.regions.size(); j++ )
    {
        segments.regions.push_back(line.regions[j]);
    }
    return;
}




void discreteCurveDimension(
    const std::vector<std::shared_ptr<GeometryCurve>>& curves,
    const std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>>& surfaces,
    const CurveParametersDimension& args, std::vector<LineMesh>& segements) {
  segements.clear();
  segements.resize(curves.size());
  int dimension = args.getdimension();
  double delta_s = args.getdelta_s();
  double angle = args.getangle();
  double deviation = args.getdeviation();
  int max_dimension = args.getmax_dimension();
  bool database = args.getuse_database_curvature();
  for (int i = 0; i < curves.size(); i++) {
    SmeshGen::LineData line;
      if( !database )
      {
          SmeshGen::discreteCurve_dimension_test(curves[i], dimension, delta_s, max_dimension, angle, deviation, line);
      }
    else
    {
      if (i>=surfaces.size()||surfaces[i].empty()) {
        std::cout << "please input topo surfaces to "<< i <<" th curve if you want to use database curvature\n";
        SmeshGen::discreteCurve_dimension_test(curves[i], dimension, delta_s, max_dimension, angle, deviation, line);
      }
      else
      {
          SmeshGen::discreteCurve_dimension_test(curves[i], surfaces[i], dimension, delta_s, max_dimension, angle,
                                                 deviation, line);
      }
    }
    for (int j = 0; j < line.coord.size(); j++) {
      segements[i].coord.push_back(line.coord[j]);
    }
    for (int j = 0; j < line.point_normal.size(); j++) {
      segements[i].point_normal.push_back(line.point_normal[j]);
    }
    for (int j = 0; j < line.regions.size(); j++) {
      segements[i].regions.push_back(line.regions[j]);
    }
  }
  return;
};

void Grid1D_Tanh(
    const std::vector<std::shared_ptr<GeometryCurve>>& curves,
    const TiGER::CurveParametersTanh& args,
    std::vector<TiGER::LineMesh>& segements) {
  segements.resize(curves.size());
  int dimension = args.getdimension();
  double begin_s = args.getbegin_s();
  double end_s = args.getend_s();
  for (int i = 0; i < curves.size(); i++) {
    SmeshGen::LineData line;
    SmeshGen::discreteCurve_tanh_test(curves[i], dimension, begin_s, end_s,
                                      line);
    for (int j = 0; j < line.coord.size(); j++) {
      segements[i].coord.push_back(line.coord[j]);
        segements[i].segments.push_back(std::array<int, 2>{j, (j+1) >= line.coord.size()?0:(j+1)});
    }
    for (int j = 0; j < line.point_normal.size(); j++) {
      segements[i].point_normal.push_back(line.point_normal[j]);
    }
    for (int j = 0; j < line.regions.size(); j++) {
      segements[i].regions.push_back(line.regions[j]);
    }
  }
}

void Grid1D_Geometric(
    const std::vector<std::shared_ptr<GeometryCurve>>& curves,
    const TiGER::CurveParametersGeometric& args,
    std::vector<TiGER::LineMesh>& segements) {
  segements.resize(curves.size());
  int dimension = args.getdimension();
  double begin_s = args.getbegin_s();
  double end_s = args.getend_s();
  for (int i = 0; i < curves.size(); i++) {
    SmeshGen::LineData line;
    SmeshGen::discreteCurve_geometric_test(curves[i], dimension, begin_s, end_s,
                                           line);
    for (int j = 0; j < line.coord.size(); j++) {
      segements[i].coord.push_back(line.coord[j]);
    }
    for (int j = 0; j < line.point_normal.size(); j++) {
      segements[i].point_normal.push_back(line.point_normal[j]);
    }
    for (int j = 0; j < line.regions.size(); j++) {
      segements[i].regions.push_back(line.regions[j]);
    }
  }
}
}  // namespace discretize_curve
namespace discretize_surface {
void CADSurf_TGrid_AFT(std::shared_ptr<GeometrySurface>& surface,
                          const std::vector<LineMesh>& boundary_segments,
                          const IsotropicSurfaceParametersTri& iso_args,
                          const SizingFunction& sizefield,
                          SurfaceMesh& output_mesh) {
   std::cout << "2024/8/29\n";
  // static int sf_id;
  bool debug = iso_args.getOutputInformation();
  if (debug)
    output_discrete_surface_info(iso_args.getSavepath(), surface,
                                 boundary_segments, iso_args, sizefield);
  SmeshGen::SurfaceData sf;
  SmeshGen::Args configure;
  configure.Transition_ratio = iso_args.getExcessRate();
  std::vector<SmeshGen::LineData> boundary;
  boundary.resize(boundary_segments.size());
  for (int i = 0; i < boundary_segments.size(); i++) {
    for (int j = 0; j < boundary_segments[i].coord.size(); j++)
      boundary[i].coord.push_back(boundary_segments[i].coord[j]);
  }
  double min_len = iso_args.getmUserSpecifiedMinEdgeLen();
  double max_len = iso_args.getmUserSpecifiedMaxEdgeLen();
  double angle = iso_args.getmMaxAngle();
  double deviation = iso_args.getmMaxDeviation();
  bool swapcell = iso_args.getmSwapCellWithNoInteriorPoints();
  // configure.remesh_time = 0;
  //configure.savepath = "D:/meshgen_win/test/smesh_wrong_information/";
  SmeshGen::CADSurf_TGrid_AFT(surface, boundary, min_len, max_len, angle,
                                 deviation, swapcell, sizefield, sf, configure);
  if (!sf.coord.size() || !sf.tris.size()) return;
  SurfaceMesh surfmesh_in;
  for (auto& pt : sf.coord) {
    surfmesh_in.coords.push_back(pt);
  }
  for (auto& tri : sf.tris) {
    surfmesh_in.tris.push_back(tri);
  }
  
  std::vector<std::vector<int>> constrained_edge=sf.boundary_segment;
  std::vector<std::vector<int>> conforming_edge;
  RemeshParameters args;
  args.setsplit_num(0);
  args.setcollapse_num(0);
  args.setb_sizefunction(true);
  args.setiteration_number(5);
  args.setb_no_inter_check(true);
  args.setb_nofeature(true);
  std::function<double(double, double, double)> size_function = sizefield;
  SurfaceMesh surfmesh_out;
  std::vector<std::vector<int>> constrained_edge_out;
  //std::vector<GeometrySurface> surfaces;
  //surfaces.push_back(*surface);
  TiGER::discretize_surface_remesh::TMeshSurf_Opt(surfmesh_in, constrained_edge, conforming_edge, args,
                                                   size_function,surfmesh_out, constrained_edge_out);
  output_mesh = surfmesh_out;

  return;
}

void VCADSurf_TGrid_AFT(
    const std::vector<std::shared_ptr<GeometrySurface>>& surfaces,
    const std::vector<std::vector<LineMesh>>& boundary_segments,
    const IsotropicSurfaceParametersTri& iso_args,
    const SizingFunction& sizefield, SurfaceMesh& output_mesh) {
  std::vector<SmeshGen::SurfaceData> sf(surfaces.size());

  SmeshGen::Args configure;
  std::vector<std::vector<SmeshGen::LineData>> boundary;
  boundary.resize(boundary_segments.size());
  for (int l = 0; l < boundary_segments.size(); l++) {
    boundary[l].resize(boundary_segments[l].size());
    for (int i = 0; i < boundary_segments[l].size(); i++) {
      for (int j = 0; j < boundary_segments[l][i].coord.size(); j++)
        boundary[l][i].coord.push_back(boundary_segments[l][i].coord[j]);
    }
  }
  double min_len = iso_args.getmUserSpecifiedMinEdgeLen();
  double max_len = iso_args.getmUserSpecifiedMaxEdgeLen();
  double angle = iso_args.getmMaxAngle();
  double deviation = iso_args.getmMaxDeviation();
  bool swapcell = iso_args.getmSwapCellWithNoInteriorPoints();
  auto surfaces_internal = (surfaces);

  // void VCADSurf_TGrid_AFT(
  // 		std::vector<std::shared_ptr<TiGER::GeometrySurface>>& surfaces_,
  // 		std::vector<std::vector<LineData>>& boundary_segments,
  // 		std::function<double(double, double, double)>sizefunction_,
  // 		std::vector<SurfaceData>& output_meshs,
  // 		Args& configure
  // 	)

  SmeshGen::VCADSurf_TGrid_AFT(surfaces_internal, boundary, min_len, max_len,
                                  angle, deviation, swapcell, sizefield, sf,
                                  configure);

  // TODO： 去重
  // igl::remove_duplicate_vertices();

  std::vector<TiGER::SurfaceMesh> all_surface(sf.size());
  for (int i = 0; i < sf.size(); i++) {
    for (auto& pt : sf[i].coord) {
      all_surface[i].coords.push_back(pt);
    }
    for (auto& tri : sf[i].tris) {
      all_surface[i].tris.push_back(tri);
    }
  }
  TiGER::merge_surface_mesh(all_surface, output_mesh);

  return;
}

void CADSurf_TGrid_AFLR(
    std::shared_ptr<GeometrySurface>& surface,
    const std::vector<std::shared_ptr<TiGER::GeometryCurve>>& curves,
    const AnisotropicSurfaceParametersTri& aniso_args,
    std::vector<LineMesh>& boundary_segments, const SizingFunction& sizefield,
    SurfaceMesh& output_mesh)
{
    if( aniso_args.getMaxLayer() <= 0 )
    {
        double max_edge_len = 0;

        auto get_len = [](const std::array<double, 3>& a, const std::array<double, 3>& b) {
            double len = 0.0;
            for( int i = 0; i < 3; ++i )
            {
                len += (a[i] - b[i]) * (a[i] - b[i]);
            }
            return std::sqrt(len);
        };

        for( int i = 0; i < boundary_segments.size(); ++i )
        {
            for( int j = 0; j < boundary_segments[i].coord.size(); ++j )
            {
                max_edge_len =
                    std::max(max_edge_len, get_len(boundary_segments[i].coord[j], boundary_segments[i].coord[0]));
            }
        }

        IsotropicSurfaceParametersTri surface_args;
        surface_args.setmUserSpecifiedMaxEdgeLen(max_edge_len * 0.1);
        return CADSurf_TGrid_AFT(surface, boundary_segments, surface_args, nullptr, output_mesh);
    }
    std::vector<int> layer_num;
    std::vector<int> match;
    std::vector<double> first_step;
  for( int i = 0; i < boundary_segments.size(); i++ )
        if( boundary_segments[i].boundary_condition.bctype == BoundaryCondition::WALL )
        {
            layer_num.push_back(i);
            first_step.push_back(boundary_segments[i].boundary_condition.delta);
        }
        else if( boundary_segments[i].boundary_condition.bctype ==
                 TiGER::BoundaryCondition::MATCH_NO_PUSH )
          match.push_back(i);
  SmeshGen::SurfaceData sf;
  SmeshGen::Args configure;
  configure.savepath = "D:/meshgen_win/test/aflr_debug/";
  configure.Topology_priority = aniso_args.getTopology_priority();
  configure.angle_threshold = aniso_args.getAngle_threshold();
  std::vector<SmeshGen::LineData> boundary;
  boundary.resize(boundary_segments.size());
  for (int i = 0; i < boundary_segments.size(); i++) {
    for (int j = 0; j < boundary_segments[i].coord.size(); j++)
      boundary[i].coord.push_back(boundary_segments[i].coord[j]);
  }
  SmeshGen::discretizeSurfaceTri_AFLR(surface, curves, layer_num, match, first_step,
                                      boundary, aniso_args.getMaxLayer(),
                                      aniso_args.getgrowth_ratio(),sizefield,sf, configure);
  if (!sf.coord.size() || !sf.tris.size()) return;
  for (auto& id : match)
  {
      boundary_segments[id].coord.clear();
      for( auto &pt : boundary[id].coord )
          boundary_segments[id].coord.push_back(pt);
  }
  for (auto& pt : sf.coord) {
    output_mesh.coords.push_back(pt);
  }
  for (auto& tri : sf.tris) {
    output_mesh.tris.push_back(tri);
  }
  return;
}
}  // namespace discretize_surface
}  // namespace TiGER
