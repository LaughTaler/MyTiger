// #define CATCH_CONFIG_MAIN // This tells Catch to provide a main() - only do
// this in one cpp file
#pragma once

#include <igl/readSTL.h>

#include <cfloat> // For DBL_MAX and DBL_MIN
#include <chrono>
#include <thread>
#include <iostream>
#include <array>
#include <iostream>
#include <memory>
#include "alias.h"
#include "spdlog/spdlog.h"
#include "list_to_matrix.hpp"
#include "tiger.h"
#include "tiger_parameters.h"
#include "../tools/meshIO.h"
#include "../tools/output_info.hpp"
#include "../tools/test_linemesh.hpp"
#include "../tools/predicates.h"
#include "geom_func.h"
#include "catch.hpp"
#include "CLI11.hpp"
// #include "test_api.h"
using namespace TiGER;
bool test_load_geometry(std::string filename, TiGER::Geometry& Geo)
{
#ifdef _MSC_VER
    filename = "Z:/home/sftp1/models/" + filename;
#else
    filename = "/home/sftp1/models/" + filename;
    // filename = "/home/zxy/test/" + filename;
#endif
    TiGER::readModelParameters readmodel_args;
    std::cout << filename<<endl;
    std::cout << Geo.curves_.size()<<endl;
    TiGER::Read_Model(filename, Geo, readmodel_args);
    std::cout << Geo.curves_.size()<<endl;
    return true; // 假设函数总是返回true，如果有错误处理逻辑请根据实际情况调整
}

void outputSurfaceMesh(const TiGER::SurfaceMesh &surf, const std::string &filename)
{
    TiGER::Mesh m;
    TiGER::list_to_matrix<3>(surf.coords, m.Vertex);
    TiGER::list_to_matrix<3>(surf.tris, m.Topo);
    TiGER::list_to_matrix<1>(surf.regions, m.Masks);
    TiGER::MESHIO::writeVTK(filename, m);
}


int main(int argc, char* argv[]) {
  
  TiGER::Parameters args;
  TiGER::Geometry Geo;
  bool if_load = test_load_geometry("cylinder.stp", Geo);
  spdlog::info("load success...\n");
  
  spdlog::info("size field...\n");
  TiGER::SizingManager sf;
  double size_scale_factor = 1;

  TiGER::CylinderSourceParameters args1;
  std::array<std::array<double, 3>, 2> xyz1 = { 
        std::array<double, 3>{0.0, 0.0, -1.0}, 
        std::array<double, 3>{0.0, 0.0, 2.0} 
  };
  std::array<double,2> r={0.5,0.5};
  TiGER::size_field::addCylinderSource(xyz1, r, args1, sf);
  // std::cout<< "formal size:"<< static_cast<TiGER::cylinderSource*>(sf.sf[0].get())->getsize(0.0, 0.0, 1) << std::endl;
  // std::cout<< "formal size:"<< static_cast<TiGER::cylinderSource*>(sf.sf[0].get())->getsize(0.0, 0.0, 3) << std::endl;
  static_cast<TiGER::cylinderSource*>(sf.sf[0].get())->args.setsourcetype(0);
  // static_cast<TiGER::cylinderSource*>(sf.sf[0].get())->args.setspacingtype(0);
  std::cout<< "constant size:"<< static_cast<TiGER::cylinderSource*>(sf.sf[0].get())->getsize(0.0, 0.0, 1) << std::endl;
  std::cout<< "constant size:"<< static_cast<TiGER::cylinderSource*>(sf.sf[0].get())->getsize(0.0, 0.5, 1) << std::endl;
  std::cout << "constant size:" << static_cast<TiGER::cylinderSource*>(sf.sf[0].get())->getsize(0.0, 4.0, 1)<< std::endl;
  std::cout << "constant size:" << static_cast<TiGER::cylinderSource*>(sf.sf[0].get())->getsize(0.0, 4.0, 3) << std::endl;
  // static_cast<TiGER::cylinderSource*>(sf.sf[0].get())->args.setspacingtype(1);
  // std::cout<< "parametric size:"<< static_cast<TiGER::cylinderSource*>(sf.sf[0].get())->getsize(0.0, 0.0, 1) << std::endl;
  // std::cout<< "parametric size:"<< static_cast<TiGER::cylinderSource*>(sf.sf[0].get())->getsize(0.0, 0.0, 3) << std::endl;
  // static_cast<TiGER::cylinderSource*>(sf.sf[0].get())->args.setspacingtype(2);
  // std::cout<< "axis size:"<< static_cast<TiGER::cylinderSource*>(sf.sf[0].get())->getsize(0.0, 0.0, 1) << std::endl;
  // std::cout<< "axis size:"<< static_cast<TiGER::cylinderSource*>(sf.sf[0].get())->getsize(0.0, 0.0, 3) << std::endl;
  // static_cast<TiGER::cylinderSource*>(sf.sf[0].get())->args.setspacingtype(3);
  // std::cout<< "center size:"<< static_cast<TiGER::cylinderSource*>(sf.sf[0].get())->getsize(0.0, 0.0, 1) << std::endl;
  // std::cout<< "center size:"<< static_cast<TiGER::cylinderSource*>(sf.sf[0].get())->getsize(0.0, 0.0, 3) << std::endl;

  TiGER::CubicSourceParameters args2;
  std::array<std::array<double, 3>, 8> xyz2 = { 
    std::array<double, 3>{0.0, 0.0, 0.0}, 
    std::array<double, 3>{0.0, 2.0, 0.0},
    std::array<double, 3>{2.0, 2.0, 0.0},
    std::array<double, 3>{2.0, 0.0, 0.0},
    std::array<double, 3>{0.0, 0.0, 2.0}, 
    std::array<double, 3>{0.0, 2.0, 2.0},
    std::array<double, 3>{2.0, 2.0, 2.0},
    std::array<double, 3>{2.0, 0.0, 2.0}
  };
  // std::array<std::array<double, 3>, 8> xyz2 = { 
  //   std::array<double, 3>{-1.0, -1.0, -1.0}, 
  //   std::array<double, 3>{-1.0, 1.0, -1.0},
  //   std::array<double, 3>{1.0, 1.0, -1.0},
  //   std::array<double, 3>{1.0, -1.0, -1.0},
  //   std::array<double, 3>{-1.0, -1.0, 4.0}, 
  //   std::array<double, 3>{-1.0, 1.0, 4.0},
  //   std::array<double, 3>{1.0, 1.0, 4.0},   
  //   std::array<double, 3>{1.0, -1.0, 4.0},
  // };
  TiGER::size_field::addCubicSource(xyz2, args2, sf);
  static_cast<TiGER::cubicSource*>(sf.sf[1].get())->args.setsourcetype(1);
  // std::cout<< "formal size:"<< static_cast<TiGER::cubicSource*>(sf.sf[0].get())->getsize(0.0, 0.0, 1) << std::endl;
  // std::cout<< "formal size:"<< static_cast<TiGER::cubicSource*>(sf.sf[0].get())->getsize(2.0, 0.0, 7) << std::endl;
  // static_cast<TiGER::cubicSource*>(sf.sf[0].get())->args.setsourcetype(1);
  // static_cast<TiGER::cubicSource*>(sf.sf[0].get())->args.setspacingtype(0);
  // std::cout << "constant size:" << static_cast<TiGER::cubicSource*>(sf.sf[0].get())->getsize(0.0, 0.0, 1) << std::endl;
  // std::cout << "constant size:" << static_cast<TiGER::cubicSource*>(sf.sf[0].get())->getsize(2.0, 0.0, 1) << std::endl;
  // std::cout << "constant size:" << static_cast<TiGER::cubicSource*>(sf.sf[0].get())->getsize(2.0, 0.0, 7)
  //           << std::endl;
  // std::cout << "constant size:" << static_cast<TiGER::cubicSource*>(sf.sf[0].get())->getsize(2.0, 2.0, 7)
  //           << std::endl;

  TiGER::PointSourceParameters args3;
  std::array<double, 3> xyz3 = {0.28, 1.97, 1.5};
  TiGER::size_field::addPointSource(xyz3, args3, sf);
  static_cast<TiGER::pointSource*>(sf.sf[2].get())->args.setsourcetype(1);

  TiGER::LineSourceParameters args4;
  std::array<std::array<double, 3>, 2> xyz4 = {std::array<double, 3>{1.38, -1.44, 1.0},
                                              std::array<double, 3>{1.38, -1.44, 8.0}};
  TiGER::size_field::addLineSource(xyz4, args4, sf);
  static_cast<TiGER::lineSource*>(sf.sf[3].get())->args.setsourcetype(1);

  TiGER::SphereSourceParameters args5;
  std::array<double, 3> xyz5 = {1.20, -1.62, 8};
  std::array<double, 1> r1 = {1.2};
  TiGER::size_field::addSphereSource(xyz5, r1, args5, sf);
  static_cast<TiGER::sphereSource*>(sf.sf[4].get())->args.setsourcetype(1);

  std::cout<< TiGER::size_field::SizingFunction_getSizeAtPoint(0, 0, 0, sf) << std::endl;

  spdlog::info("Add source succeed");
  spdlog::info(sf.sf.size());

  spdlog::info("Start discretize...\n");
  std::vector<LineMesh> segements;
  // SizingFunction size_field =[](const double& x, const double& y,const double& z){return 0.2;};
  SizingFunction size_field = [&](const double& x, const double& y,
                                 const double& z) { return 0.2*static_cast<TiGER::cylinderSource*>(sf.sf[0].get())->getsize(x, y, z); };
  TiGER::CurveParametersTanh args_curve;
  args_curve.setdimension(10);
  spdlog::info("Start discretize curves...\n");
  TiGER::discretize_curve::Grid1D_Tanh(Geo.curves_, args_curve, segements);
  spdlog::info("End discretize curves...\n");
  
  spdlog::info(segements[0].coord.size() );
  
  std::vector<std::vector<LineMesh>> sf_boundary;
  sf_boundary.resize(Geo.surfaces_.size());
  for (int i = 0; i < Geo.topo_surf_to_curves_.size(); i++) {
    for (int j = 0; j < Geo.topo_surf_to_curves_[i].size(); j++) {
      sf_boundary[i].push_back(segements[Geo.topo_surf_to_curves_[i][j]]);
    }
  }
  std::vector<SurfaceMesh> sf_mesh(sf_boundary.size());

  for (int k = 0; k < sf_boundary.size(); k++) {
    spdlog::info("discretize surfaces...\n");
    std::cout << k << std::endl;
    spdlog::info("discretize surfaces...\n");
    std::cout << k << std::endl;
    TiGER::discretize_surface::CADSurf_TGrid_AFT(
        Geo.surfaces_[k], sf_boundary[k], args, size_field, sf_mesh[k]);
  }
  TiGER::SurfaceMesh sf_output;
  TiGER::merge_surface_mesh(sf_mesh, sf_output);

  spdlog::info("discretize success...\n");

  outputSurfaceMesh(sf_output, "surface_mesh_source.vtk");


  // 尺寸场
  // spdlog::info("size field...\n");
  // TiGER::SizingManager sf;
  // spdlog::info(sf.sf.size());
  // spdlog::info(sf_output.coords.size());
  // spdlog::info(sf_output.tris.size());
  // spdlog::info(sf_output.regions.size());
  // spdlog::info("before setSizeField\n");
  // TiGER::size_field::SizingFunction_setUniformSize(1);
  // spdlog::info("after setSizeField\n");
  // spdlog::info(sf.sf.size());

  

  

  // TiGER::BackGroundParameters bkg_args;
  // double size_hmin = -1;
  // double size_hmax = -1;
  // double size_theta = 0.1;
  // double size_beta = 1.2;
  // double size_proximity = 2;
  // bkg_args.sethmin(size_hmin);
  // bkg_args.sethmax(size_hmax);
  // bkg_args.settheta(size_theta);
  // bkg_args.setbeta(size_beta);
  // bkg_args.setproximity(size_proximity);
  // TiGER::size_field::setSizeField(sf_output, bkg_args, sf);


  

  
  // spdlog::info("Using size field to remesh...\n");
  // TiGER::RemeshParameters args_remesh;
  // double remesh_epsilon = 1e-0;
  // int remesh_iterater = 5;
  // args_remesh.setiteration_number(remesh_iterater);
  // args_remesh.setb_sizefunction(false);
  // args_remesh.seteps_rel(remesh_epsilon);  
  // std::vector<std::vector<int>> constrained_edge;
  // std::vector<std::vector<int>> conformaing_edge;
  // std::vector<std::vector<int>> constrained_edge_out;
  // std::function<double(double, double, double)> scaled_size_function = [&](double x, double y, double z)
  // {
  //     // return size_scale_factor * TiGER::size_field::SizingFunction_getSizeAtPoint(x, y, z, sf);
  //     // return size_scale_factor * static_cast<TiGER::sphereSource*>(sf.sf[0].get())->getsize(x, y, z); 
  //     return size_scale_factor;
  //     // return size_scale_factor;
  // };
  // SurfaceMesh remesh_surf_out;
  // TiGER::discretize_surface_remesh::TMeshSurf_Opt(
  //     sf_output, constrained_edge, conformaing_edge, args_remesh,
  //     scaled_size_function, remesh_surf_out, constrained_edge_out);
  


  // outputSurfaceMesh(remesh_surf_out, "surface_mesh2.vtk");

  // 生成四面体
  TiGER::VolumeMesh vol;
  TiGER::tetrahedral_mesh::TetGrid_Constrained(sf_output, args, vol, NULL);
  // args.setnthread(1);
  // args.setnthread(1);
  // 检查质量
  const TiGER::VolMeshQualityType type = TiGER::VolMeshQualityType::TET_VOLUME;
  std::vector<TiGER::QualityResult> quality;
  TiGER::Tool_VolQuality(vol, type, quality);
  for (auto& i : quality) {
    spdlog::info(i.ave_value);
    spdlog::info(i.max_value);
    spdlog::info(i.min_value);
  }

  // 转换成mesh
  TiGER::Mesh m1;
  TiGER::list_to_matrix<3>(vol.coords, m1.Vertex);
  TiGER::list_to_matrix<4>(vol.tetras, m1.Topo);
  TiGER::MESHIO::writeVTK("vol_mesh.vtk", m1);
}