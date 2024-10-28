#include "boundary_layer_mesh.h"

#include <exception>
#include <iostream>
#include <limits>
#include <vector>
#include <spdlog/spdlog.h>

#include "./geometry_data_struct.h"
#include "tetrahedral_mesh.h"
#include "vmesh.h"
#include "alias.h"

namespace TiGER {
namespace boundary_layer {

void convertToIntPointer(const SurfaceMesh& mesh, int** trisPtr,
                         int** regionsPtr);
void convertToDoublePointer(const SurfaceMesh& mesh, double** coordPtr);

void HexGrid_createHybridMesh(
    const SurfaceMesh& boundary_mesh,
    const std::vector<std::shared_ptr<TiGER::GeometrySurface>>& surfaces,
    const HybridParameters& args, VolumeMesh& vol_mesh,
    SurfaceMesh& out_boundary_mesh) {
  SurfaceMesh tris_boundary_mesh;
  VolumeMesh bdvol_mesh;
  std::vector<int> l2g;
  HexGrid_createPrismaticMesh(boundary_mesh, surfaces, args, bdvol_mesh,
                        out_boundary_mesh, tris_boundary_mesh, l2g);

  VolumeMesh tetra_mesh;
  SizingFunction func = [](const double& x, const double& y, const double& z) {
    return std::numeric_limits<double>::max();
  };
  tetrahedral_mesh::TetGrid_Constrained(tris_boundary_mesh, args, tetra_mesh,
                                      func);
  if (std::vector<std::array<double, 3>>(
          tetra_mesh.coords.begin(),
          tetra_mesh.coords.begin() + tris_boundary_mesh.coords.size()) !=
      tris_boundary_mesh.coords) {
    throw std::runtime_error("tetra failed, non-constrain mesh");
  } else {
    spdlog::info("coordinate matched.");
  }
  // for (int i = 0; i < tetra_mesh.tetras.size(); i++) {
  //   vol_mesh.tetras.push_back(std::array<int, 4>{});
  // }
  vol_mesh = bdvol_mesh;
  for (int i = tris_boundary_mesh.coords.size(); i < tetra_mesh.coords.size();
       i++) {
    vol_mesh.coords.push_back(tetra_mesh.coords[i]);
  }
  for (int i = 0; i < tetra_mesh.tetras.size(); i++) {
    std::array<int, 4> tetra;
    bool b = false;
    for (int k = 0; k < 4; k++) {
      
      int v = -1;
      if (tetra_mesh.tetras[i][k] < tris_boundary_mesh.coords.size()) {
        v = l2g[tetra_mesh.tetras[i][k]];
        b = true;
      } else {
        v = tetra_mesh.tetras[i][k] + bdvol_mesh.coords.size() -
            tris_boundary_mesh.coords.size();
      }
      tetra[k] = (v);
    }
    vol_mesh.tetras.push_back(tetra);
  }
}

void HexGrid_createPrismaticMesh(
    const SurfaceMesh& boundary_mesh,
    const std::vector<std::shared_ptr<TiGER::GeometrySurface>>& surfaces,
    const BoundaryLayerParameters& args, VolumeMesh& vol_mesh,
    SurfaceMesh& out_boundary_mesh, SurfaceMesh& tris_boundary_mesh,
    std::vector<int>& tris_boundary_point_l2g) {

  double* pdSNC; /* 曲面网格节点坐标，coord[3*i], coord[3*i+1],
                    coord[3*i+2]为第i个节点x,y,z坐标 **/
  convertToDoublePointer(boundary_mesh, &pdSNC);

  // 实现 HexGrid_createPrismaticMesh 函数的代码
  /* ------------------------- 输入参数 --------------------------**/
  int nSN = boundary_mesh.coords.size(); /* 曲面网格边界节点数目 **/
  int* pnSFFm;                           /* 曲面网格单元节点编号 **/
  int* pnSFTp = nullptr; /* 曲面网格单元类型。当前仅支持三角形单元 **/
  int* pnSFPt; /* 曲面网格单元所在几何面编号 ,从1开始 **/
  convertToIntPointer(boundary_mesh, &pnSFFm, &pnSFPt);
  int nSF = boundary_mesh.tris.size(); /* 曲面网格单元数目 **/
  std::vector<int> pnFT(boundary_mesh.bc.rbegin()->first+1 ,
                        0); /*  几何面类型： 0为远场； 1 为物面； 2为对称面
                               ,3为不长边界层的面,4为周期性面 **/
  for (auto k : boundary_mesh.bc) {

      if (k.first < 0) {
          throw std::runtime_error("negative bc face id");
      }
    switch (k.second) {
      case TiGER::BoundaryCondition::WALL:
        pnFT[k.first] = 1;
        break;
      case TiGER::BoundaryCondition::SYMMETRY:
        pnFT[k.first] = 2;
        break;
      case TiGER::BoundaryCondition::FARFIELD:
        pnFT[k.first] = 0;
        break;
      case TiGER::BoundaryCondition::ADJACENT_GRID:
        pnFT[k.first] = 5;
        break;
      case TiGER::BoundaryCondition::MATCH_NO_PUSH:
        pnFT[k.first] = 3;
        break;
      default:
        break;
    }
    std::cout << k.first << " " << k.second << " " << pnFT[k.first]
              << std::endl;
    std::cout.flush();
  }
  int nLN = args.getnumber_of_layer(); /* 边界层层数 **/
  double dLen = args.getstep_length(); /* 边界层第一层厚度 **/
  double dRto = args.getratio();       /* 边界层厚度增长因子 **/
  double bisostop = args.getaniso_iso_blend(); /* 各向同性停止**/
  /* ------------------------- 输出参数 -------------------------**/
  double* ppdMNC; /* 体网格节点坐标 **/
  int pnMN;       /* 体网格节点数目 **/
  int* ppnMEFm;   /* 体网格单元节点编号 **/
  int* ppnMETp; /* 体网格单元类型。当前支持单元类型：三棱柱，金字塔 **/
  int pnME; /* 体网格单元数目 **/
  /* ------------------------- 边界层网格顶面层面网格参数
   * -------------------------**/
  double* ppdSNC0; /* 顶面网格节点坐标，coord[3*i], coord[3*i+1],
                      coord[3*i+2]为第i个节点x,y,z坐标 **/
  int pnSN0;       /* 顶面网格边界节点数目 **/
  int pnSEO;       /* 顶面网格单元数目 **/
  int* ppnSFTpO; /* 顶面网格单元类型。当前仅支持三角形单元 **/
  int* ppnSFFmO; /* 体网格单元数目 **/
  int* l2g;      /* local点id到global点id的映射 */

  double* point_sizing;  /*  顶面网格点目标尺寸 */

  /* ------------------------- 其他参数 -------------------------**/
  bool b_have_pyramid = true; /*是否有金字塔**/
  bool b_use_multiple_normals =
      args.getmultiple_normal(); /*是否启用多法向 缺省为false**/
  bool b_output_io_file =
      args.getoutputio(); /*是否将api的输入和输出都输出到文件系统中（仅用于DEBUG）缺省为false**/
  std::string filename;   /*几何文件名，缺省为virtualmesh**/

  std::array<double, 12> per_matrix = {
      0}; /*周期性面控制矩阵,前9位为旋转矩阵 m00，m01，m02 ....
             ，后三位为位移向量xyz**/
  int num_boundary_face;  /////* 边界面网格数量 ***/
  int*
      boundary_t_mesh;  /////*
                        /// 边界面网格，注意每四个为一组，而且注意如果为三角形，最后一项为-1，id从0开始
                        ///***/
  int* boundary_face;  /////* 边界面网格对应的面id，长度为num_boundary_face ***/

  TiGER::API_Gen_Boundary_Mesh(
      /* ------------------------- 输入参数 --------------------------**/
      pdSNC,  /* 曲面网格节点坐标，coord[3*i], coord[3*i+1],
                 coord[3*i+2]为第i个节点x,y,z坐标 **/
      nSN,    /* 曲面网格边界节点数目 **/
      pnSFFm, /* 曲面网格单元节点编号 **/
      pnSFTp, /* 曲面网格单元类型。当前仅支持三角形单元 **/
      pnSFPt, /* 曲面网格单元所在几何面编号 ,从1开始 **/
      nSF,    /* 曲面网格单元数目 **/
      pnFT.data(), /*  几何面类型： 0为远场； 1 为物面； 2为对称面
                      ,3为不长边界层的面,4为周期性面.5为临近单元 **/
      nLN,         /* 边界层层数 **/
      dLen,        /* 边界层第一层厚度 **/
      dRto,        /* 边界层厚度增长因子 **/
      bisostop, /* 各向同性停止**/
      /* ------------------------- 输出参数 -------------------------**/
      &ppdMNC,  /* 体网格节点坐标 **/
      &pnMN,    /* 体网格节点数目 **/
      &ppnMEFm, /* 体网格单元节点编号 **/
      &ppnMETp, /* 体网格单元类型。当前支持单元类型：三棱柱，金字塔 **/
      &pnME, /* 体网格单元数目 **/
      /* ------------------------- 边界层网格顶面层面网格参数
         -------------------------**/
      &ppdSNC0,  /* 顶面网格节点坐标，coord[3*i], coord[3*i+1],
                    coord[3*i+2]为第i个节点x,y,z坐标 **/
      &pnSN0,    /* 顶面网格边界节点数目 **/
      &pnSEO,    /* 顶面网格单元数目 **/
      &ppnSFTpO, /* 顶面网格单元类型。当前仅支持三角形单元 **/
      &ppnSFFmO, /* 顶面网格单元节点编号 **/
      &l2g,      /* local点id到global点id的映射 */
      &point_sizing, /* 点尺寸 */
      /////* ------------------------- 边界信息 ------------------------- ***/
      &num_boundary_face,  /////* 边界面网格数量 ***/
      &boundary_t_mesh,    /////*
                         /// 边界面网格，注意每四个为一组，而且注意如果为三角形，最后一项为-1，id从0开始
                         ///***/
      &boundary_face,  /////* 边界面网格对应的面id，长度为num_boundary_face ***/
      /* ------------------------- 其他参数 -------------------------**/
      b_have_pyramid,         /*是否有金字塔**/
      b_use_multiple_normals, /*是否启用多法向 缺省为false**/

      b_output_io_file, /*是否将api的输入和输出都输出到文件系统中（仅用于DEBUG）缺省为false**/

      filename, /*几何文件名，缺省为virtualmesh**/

      per_matrix /*周期性面控制矩阵,前9位为旋转矩阵 m00，m01，m02 ....
                    ，后三位为位移向量xyz**/
  );
  tris_boundary_point_l2g.resize(pnSN0);
  for (int i = 0; i < pnSN0; i++) {
    tris_boundary_point_l2g[i] = l2g[i];
  }

  /* 体网格赋值模块 **/
  int count = 0;
  for (int i = 0; i < pnSN0; i++) {
    tris_boundary_mesh.coords.push_back(
        {ppdSNC0[count], ppdSNC0[count + 1], ppdSNC0[count + 2]});
    count += 3;
  }
  count = 0;
  for (int i = 0; i < pnSEO; i++) {
    tris_boundary_mesh.tris.push_back(
        {ppnSFFmO[count], ppnSFFmO[count + 1], ppnSFFmO[count + 2]});
    count += 3;
  }
  tris_boundary_mesh.point_attribute_double.resize(tris_boundary_mesh.coords.size());
  for (int i = 0; i < tris_boundary_mesh.point_attribute_double.size(); i++) {
    tris_boundary_mesh.point_attribute_double[i] = point_sizing[i];
  }
  /* 最外层网格赋值模块 **/
  count = 0;
  for (int i = 0; i < pnMN; i++) {
    vol_mesh.coords.push_back(
        {ppdMNC[count], ppdMNC[count + 1], ppdMNC[count + 2]});
    count += 3;
  }
  count = 0;
  for (int i = 0; i < pnME; i++) {
    if (ppnMETp[i] == TiGER::EntityTopology::TETRAHEDRON) {
      vol_mesh.tetras.push_back({ppnMEFm[count], ppnMEFm[count + 1],
                                 ppnMEFm[count + 2], ppnMEFm[count + 3]});
      count += 4;
    }
    if (ppnMETp[i] == TiGER::EntityTopology::PRISM) {
      vol_mesh.prisms.push_back({ppnMEFm[count], ppnMEFm[count + 1],
                                 ppnMEFm[count + 2], ppnMEFm[count + 3],
                                 ppnMEFm[count + 4], ppnMEFm[count + 5]});
      count += 6;
    }
    if (ppnMETp[i] == TiGER::EntityTopology::PYRAMID) {
      vol_mesh.pyramids.push_back({ppnMEFm[count], ppnMEFm[count + 1],
                                   ppnMEFm[count + 2], ppnMEFm[count + 3],
                                   ppnMEFm[count + 4]});
      count += 5;
    }
  }

  out_boundary_mesh.coords = vol_mesh.coords;
  for (int i = 0; i < num_boundary_face; i++) {
    if (boundary_t_mesh[4 * i + 3] == -1) {
      out_boundary_mesh.tris.push_back({boundary_t_mesh[4 * i + 0],
                                        boundary_t_mesh[4 * i + 1],
                                        boundary_t_mesh[4 * i + 2]});
      out_boundary_mesh.attribute_int.push_back(boundary_face[i]);
    }
  }
  for (int i = 0; i < num_boundary_face; i++) {
    if (boundary_t_mesh[4 * i + 3] != -1) {
      out_boundary_mesh.quads.push_back(
          {boundary_t_mesh[4 * i + 0], boundary_t_mesh[4 * i + 1],
           boundary_t_mesh[4 * i + 2], boundary_t_mesh[4 * i + 3]});
      out_boundary_mesh.attribute_int.push_back(boundary_face[i]);
    }
  }
}

void HexGrid_Legacy (const SurfaceMesh& boundary_mesh,
                     const std::vector<BoundaryCondition>& boundary_conditions,
                     const std::vector<GeometrySurface>& surfaces,
                     const BoundaryLayerParameters& args, VolumeMesh& vol_mesh,
                     SurfaceMesh& out_boundary_mesh) {
  // 实现 HexGrid_Legacy  函数的代码
}

void convertToDoublePointer(const SurfaceMesh& mesh, double** coordPtr) {
  // Convert coord to double*
  *coordPtr = new double[mesh.coords.size() * 3];
  int index = 0;
  for (const auto& point : mesh.coords) {
    (*coordPtr)[index++] = point[0];
    (*coordPtr)[index++] = point[1];
    (*coordPtr)[index++] = point[2];
  }
}

void convertToIntPointer(const SurfaceMesh& mesh, int** trisPtr,
                         int** regionsPtr) {
  // Convert tris to int*
  *trisPtr = new int[mesh.tris.size() * 3];
  int index = 0;
  for (const auto& tri : mesh.tris) {
    (*trisPtr)[index++] = tri[0] + 1;
    (*trisPtr)[index++] = tri[1] + 1;
    (*trisPtr)[index++] = tri[2] + 1;
  }

  // Convert regions to int*
  if (mesh.attribute_int.size() == mesh.tris.size()) {
    *regionsPtr = new int[mesh.attribute_int.size()];
    for (size_t i = 0; i < mesh.attribute_int.size(); ++i) {
      (*regionsPtr)[i] = mesh.attribute_int[i];
    }
  }
  if (mesh.regions.size() == mesh.tris.size()) {
    *regionsPtr = new int[mesh.regions.size()];
    for (size_t i = 0; i < mesh.regions.size(); ++i) {
      (*regionsPtr)[i] = mesh.regions[i];
    }
  }
}

}  // namespace boundary_layer
}  // namespace TiGER
