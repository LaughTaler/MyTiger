/**
 * @author :yhf
*/
#include "../include/mesh_data_struct.h"
#include "./list_to_matrix.hpp"
#include "mesh_repair.h"
#include "binary_tree.hpp"
#include "igl/boundary_loop.h"
#include "igl/remove_duplicate_vertices.h"
#include <unordered_map>
#include <iostream>
#include"mesh_repair.h"
#include "alias.h"
namespace TiGER {

bool arePointsEqual(const std::array<double, 3>& p1,
                    const std::array<double, 3>& p2, double tol = 1e-6) {
  return std::abs(p1[0] - p2[0]) < tol && std::abs(p1[1] - p2[1]) < tol &&
         std::abs(p1[2] - p2[2]) < tol;
}

int findset(std::vector<int>& belong, const int x) {
  if (belong[x] != x) {
    belong[x] = findset(belong, belong[x]);  // ·��ѹ��
  }
  return x;
}
void unionset(std::vector<int>& belong, const int a, const int b) {
  int root_a = findset(belong, a);
  int root_b = findset(belong, b);
  if (root_a != root_b) {
    belong[root_b] = root_a;
  }
}

int detectduplicatepoints(const std::vector<std::array<double, 3>>& points_in,
                          const double& eps, std::vector<int>& bcj_pre) {
  int count = 0;
  bcj_pre.resize(points_in.size());
  std::vector<int> belong(points_in.size());
  for (int i = 0; i < points_in.size(); i++) {
    belong[i] = i;  // ��ʼ�� ÿ���ڵ㶼����Ϊһ��
  }

  TiGER::common::BinaryTree<Eigen::Matrix3d> envelop_tree_;
  Eigen::MatrixXd V(points_in.size(), 3);
  for (int i = 0; i < points_in.size(); i++) {
    V.row(i) =
        Eigen::RowVector3d(points_in[i][0], points_in[i][1], points_in[i][2]);
  }
  envelop_tree_.init(V);
  for (int i = 0; i < points_in.size(); i++) {
    Eigen::Vector3d v1(points_in[i][0], points_in[i][1], points_in[i][2]);
    std::vector<size_t> result;
    envelop_tree_.query(V.row(i), eps, result);
    for (int j = 0; j < result.size();
         j++) {  // result �����ĵĵ�box�루�ݲΧbox���ཻ
      // �ݲΧ�ڵĵ�
      Eigen::Vector3d v2(points_in[result[j]][0], points_in[result[j]][1],
                         points_in[result[j]][2]);
      double dist = (v1 - v2).norm();
      if (dist <= eps) {
        if (belong[i] == i && result[j] > i &&
            belong[result[j]] == result[j]) {  // ��ֹ�ظ�����
          // std::cout<<i<<" "<<result[j]<<" "<<dist<<endl;
          unionset(belong, i, result[j]);  // �ϲ���ǰһ���ڵ�ļ���
          count++;
        }
      }
    }
  }
  for (int i = 0; i < points_in.size(); i++) {
    if (belong[i] == i) {
      bcj_pre[i] = -1;
    } else {
      bcj_pre[i] = belong[i];
    }
  }

  return count;
}


void merge_surface_mesh(const std::vector<SurfaceMesh>& input,
                        SurfaceMesh& output) {
    std::vector<std::array<double, 3>> pt_global;
    for (int i = 0; i < input.size(); i++)
        for (auto pt : input[i].coords) {
            pt_global.push_back(pt);
            // output.coords.push_back(pt);
        }
    int cnt = 0;
    std::vector<std::array<int, 3>> tri_global;
    std::vector<int> surface_id;
    for (int i = 0; i < input.size(); i++) {
        for (auto tri : input[i].tris) {
            std::array<int, 3> new_tri = { tri[0] + cnt, tri[1] + cnt, tri[2] + cnt };
            tri_global.push_back(new_tri);
            surface_id.push_back(i);         // output.tris.push_back(new_tri);
        }
        cnt += input[i].coords.size();
    }
    // return;
    double eps = 1e-4;
    std::vector<int> bcj_pre(pt_global.size(), -1);
    TiGER::repair::Tool_detectPointDuplicates(pt_global, eps, bcj_pre);
    std::unordered_map<int, int> pt_index;
    int index = 0;

    for (int i = 0; i < bcj_pre.size(); i++)
    {
        if (bcj_pre[i] == -1)
        {
            output.coords.push_back(pt_global[i]);
            pt_index[i] = index;
            index++;
        }
    }
    for (int i = 0; i < tri_global.size(); i++)
    {
        std::array<int, 3> tri;
        for (int j = 0; j < 3; j++)
        {
            if (bcj_pre[tri_global[i][j]] != -1)
                tri[j] = pt_index[bcj_pre[tri_global[i][j]]];
            else
                tri[j] = pt_index[tri_global[i][j]];
        }
        output.tris.push_back(tri);

        //output.attribute_int.push_back(i);
    }
    output.attribute_int = surface_id;
    return;

  /*std::vector<std::array<double, 3>> pt_global;
  for (int i = 0; i < input.size(); i++)
    for (auto pt : input[i].coords) {
      pt_global.push_back(pt);
      // output.coords.push_back(pt);
    }
  int cnt = 0;
  std::vector<std::array<int, 3>> tri_global;
  for (int i = 0; i < input.size(); i++) {
    for (auto tri : input[i].tris) {
      std::array<int, 3> new_tri = {tri[0] + cnt, tri[1] + cnt, tri[2] + cnt};
      tri_global.push_back(new_tri);
      // output.tris.push_back(new_tri);
    }
    cnt += input[i].coords.size();
  }
  return;
  double eps = 1e-4;
  std::vector<int> bcj_pre(pt_global.size(),-1);
  TiGER::repair::detectDuplicatePoints(pt_global, eps, bcj_pre);
  std::unordered_map<int, int> pt_index;
  int index = 0;

  for (int i = 0; i < bcj_pre.size(); i++)
  {
      if (bcj_pre[i] == -1)
      {
        output.coords.push_back(pt_global[i]);
        pt_index[i] = index;
        index++;
      }
  }
  for (int i = 0; i < tri_global.size(); i++)
  {
    std::array<int, 3> tri;
    for (int j = 0; j < 3; j++)
    {
      if (bcj_pre[tri_global[i][j]] != -1)
        tri[j] = pt_index[bcj_pre[tri_global[i][j]]];
      else
        tri[j] = pt_index[tri_global[i][j]];
    }
    output.tris.push_back(tri);
  }*/

  double ave_eps = 0;
  int loop_count = 0;
  for (const auto& mesh : input) {
    Eigen::MatrixXi topo;
    TiGER::list_to_matrix<3>(mesh.tris, topo);
    Eigen::MatrixXi loop;
    igl::boundary_loop(topo, loop);

    for (int k = 0; k < loop.rows() - 1; k++) {
      ave_eps += (Eigen::RowVector3d(mesh.coords[loop(k, 0)].data()) -
                  Eigen::RowVector3d(mesh.coords[loop(k + 1, 0)].data()))
                     .norm();
    }
    if (loop.rows() > 1) {
      ave_eps /= loop.rows();
      loop_count += loop.rows();
    }
  }
  if (loop_count > 0) {
    ave_eps /= loop_count;
  }
  //double eps = ave_eps / 1000;

  std::map<std::array<double, 3>, int> vertexIndexMap;
  int newIndex = 0;

  for (const auto& mesh : input) {
    for (const auto& coord : mesh.coords) {
      bool found = false;
      int existingIndex = -1;

      for (const auto& item : vertexIndexMap) {
        if (std::sqrt(std::pow(item.first[0] - coord[0], 2) +
                      std::pow(item.first[1] - coord[1], 2) +
                      std::pow(item.first[2] - coord[2], 2)) < eps) {
          found = true;
          existingIndex = item.second;
          break;
        }
      }
      if (!found) {
        vertexIndexMap[coord] = newIndex;
        output.coords.push_back(coord);
        existingIndex = newIndex++;
      }
    }
  }

  // Now update the indices for tris in the output
  int faceid=0;
  for (const auto& mesh : input) {
    for (const auto& tri : mesh.tris) {
      std::array<int, 3> newTri;
      for (int i = 0; i < 3; ++i) {
        newTri[i] = vertexIndexMap[mesh.coords[tri[i]]];
      }
      output.tris.push_back(newTri);
      output.attribute_int.push_back(faceid);
    }
    faceid++;
  }
  // std::cout << "Total unique vertices after merging: " << output.coords.size()
  //           << std::endl;
}

}  // namespace TiGER