#pragma once

// #include "mesh_data_struct.h"
// #include "definitions.h"
// #include <functional>
// #include <vector>
#include "mesh_repair.h"
#include "binary_tree.hpp"
// #include "predicates.cpp"
#include <geom_func.h>

#include <unordered_map>
#include <omp.h>
#include <queue>
#include "alias.h"
#include "predicates.h"

typedef double REAL;

namespace TiGER {
namespace repair {

int find(std::vector<int>& belong, const int x) {
  if (belong[x] != x) {
    belong[x] = find(belong, belong[x]);  // 路径压缩
  }
  return x;
}
void unionSet(std::vector<int>& belong, const int a, const int b) {
  int root_a = find(belong, a);
  int root_b = find(belong, b);
  if (root_a != root_b) {
    belong[root_b] = root_a;
  }
}

/**
 * @brief
 *
 * @param points_in 输入的点
 * @param eps 去重的容差
 * @param bcj_pre 并查集
 *
 * points_in index =  0   1   2  3  4  5  6   7  /// 2，3，5，6 重复，4，7重复。
 * bcj_pre         = -1  -1  -1  2 -1  2  2   4  --> 并查集结果
 *
 * @return int 重复点的个数
 */
int Tool_detectPointDuplicates(
    const std::vector<std::array<double, 3>>& points_in,
    const double& eps, 
    std::vector<int>& bcj_pre){

        // TiGER::Mesh mesh_point;
        // mesh_point.Vertex.resize(points_in.size(),3);
        // for(int i=0;i<points_in.size();i++){
        //     for(int j=0;j<3;j++){
        //         mesh_point.Vertex(i,j) = points_in[i][j];
        //     }
        // }
        // TiGER::MESHIO::writeVTK("C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/error/dlr-f4_stl_point.vtk",mesh_point);

        bcj_pre.resize(points_in.size());

        std::array<double,6> box;
        for(size_t k=0;k<3;k++){
            box[k] = std::numeric_limits<double>::max();
            box[k+3] = -std::numeric_limits<double>::max();
        }
        for(int i=0;i<points_in.size();i++){
            for(size_t k=0;k<3;k++){
                box[k] = std::min(box[k],points_in[i][k]);
                box[k+3] = std::max(box[k+3],points_in[i][k]);
            }
        }
        TiGER::common::BinaryTree<Eigen::RowVector3d> envelop_tree_(box);
        size_t before = -1;
        std::vector<size_t> new_old_map(points_in.size());
        size_t cnt = 0;
        for(int i=0;i<points_in.size();i++){
            Eigen::RowVector3d v = Eigen::RowVector3d(points_in[i][0],
                                            points_in[i][1],
                                            points_in[i][2]);
            std::vector<size_t> result;
            envelop_tree_.query(v,eps,result,before);  // 
            // bool flag_dup = false;
            // for(int k=0;k<result.size();k++){
            //     Eigen::RowVector3d u = Eigen::RowVector3d(points_in[new_old_map[result[k]]][0],
            //                                             points_in[new_old_map[result[k]]][1],
            //                                             points_in[new_old_map[result[k]]][2]);
            //     double dist = (u-v).norm();
            //     if( dist <= eps){
            //         bcj_pre[i] = new_old_map[result[k]];
            //         flag_dup = true;
            //         break;
            //     }
            // }
            // if(flag_dup == false){
            //     envelop_tree_.insert(v,cnt);
            //     new_old_map[cnt] = i;
            //     bcj_pre[i] = -1;
            //     cnt++;
            // }

            if(result.size()){
                // if (result.size() > 1)
                //     throw std::runtime_error("too many duplicate points");
                bcj_pre[i] = new_old_map[result[0]];
            }
            else{
                envelop_tree_.insert(v,cnt);
                new_old_map[cnt] = i;
                bcj_pre[i] = -1;
                cnt++;
            }
        }
        return points_in.size() - cnt;
    }

int Tool_nativelyDetectPointDuplicates(
    const std::vector<std::array<double, 3>>& points_in,
    const double& eps, 
    std::vector<int>& bcj_pre){
        int count = 0;
        bcj_pre.resize(points_in.size());
        std::vector<int> belong(points_in.size());
        for(int i=0;i<points_in.size();i++){
            belong[i] = i; // 初始化 每个节点都各自为一组
        }
        // 暴力
        for(int i=0;i<points_in.size();i++){
            Eigen::Vector3d v1(points_in[i][0],points_in[i][1],points_in[i][2]);
            for(int j=i+1;j<points_in.size();j++){
                Eigen::Vector3d v2(points_in[j][0],points_in[j][1],points_in[j][2]);
                double dist = (v1-v2).norm();
                if(dist <= eps && belong[i] == i && belong[j] == j){  // 重复点
                    unionSet(belong,i,j);  // 合并到前一个节点的集合
                    count++;
                }
            }
        }
        for(int i=0;i<points_in.size();i++){
            if(belong[i] == i){
                bcj_pre[i] = -1;
            }
            else{
                bcj_pre[i] = belong[i];
            }
        }
        return count;
    }

void Tool_projectPointToSegment(
    const Eigen::RowVector3d& p,
    const Eigen::RowVector3d& a,
    const Eigen::RowVector3d& b,
    Eigen::RowVector3d& project_point,
    double& point_to_segment_distance){
        Eigen::RowVector3d v1 = b - a;
        Eigen::RowVector3d v2 = p - a;
        double v1_dot_v2 = v1.dot(v2);
        double v1_dot_v1 = v1.dot(v1);
        double t = v1_dot_v2 / v1_dot_v1;  // 比例
        if(t < 0.0){
            project_point = a;
        }
        else if(t > 1.0){
            project_point = b;
        }
        else{
            project_point = a + t * v1;
        }
        point_to_segment_distance = (project_point - p).norm();
    }

void Tool_projectPointToTriangle(
    const Eigen::RowVector3d& p,
    const Eigen::RowVector3d& a,
    const Eigen::RowVector3d& b,
    const Eigen::RowVector3d& c,
    Eigen::RowVector3d& project_point,
    double& point_to_tri_distance){
        // 三角片的法向
        Eigen::RowVector3d normal = (b - a).cross(c - a).normalized();
        // p 点在平面的投影点
        double distance = (p - a).dot(normal);
        Eigen::RowVector3d projection = p - distance * normal;
        // 是否在三角形内部
        Eigen::RowVector3d v1 = b - a;
        Eigen::RowVector3d v2 = c - a;
        Eigen::RowVector3d v3 = projection - a;

        double v1_dot_v1 = v1.dot(v1);
        double v1_dot_v2 = v1.dot(v2);
        double v1_dot_v3 = v1.dot(v3);
        double v2_dot_v2 = v2.dot(v2);
        double v2_dot_v3 = v2.dot(v3);

        double invDenom = 1.0 / (v1_dot_v1 * v2_dot_v2 - v1_dot_v2 * v1_dot_v2);
        double u = (v2_dot_v2 * v1_dot_v3 - v1_dot_v2 * v2_dot_v3) * invDenom;
        double v = (v1_dot_v1 * v2_dot_v3 - v1_dot_v2 * v1_dot_v3) * invDenom;

        // 投影点在三角形内部
        if((u >= 0) && (v >= 0) && (u + v <= 1)){
            project_point = projection;
            point_to_tri_distance = (p - project_point).norm();
            return;
        }

        // 不在三角形内部
        Tool_projectPointToSegment(projection,a,b,project_point,point_to_tri_distance);
        
        Eigen::RowVector3d point;
        double dist;
        Tool_projectPointToSegment(projection,b,c,point,dist);

        if(dist < point_to_tri_distance){
            point_to_tri_distance = dist;
            project_point = point;
        }

        Tool_projectPointToSegment(projection,c,a,point,dist);
        if(dist < point_to_tri_distance){
            point_to_tri_distance = dist;
            project_point = point;
        }
        point_to_tri_distance = (p - project_point).norm();
    }

void project_SegmentOntoTriangle(
    const Eigen::RowVector3d& e1,
    const Eigen::RowVector3d& e2,
    const Eigen::RowVector3d& a,
    const Eigen::RowVector3d& b,
    const Eigen::RowVector3d& c,
    std::vector<Eigen::RowVector3d>& project_edge){
        // 三角片的法向
        Eigen::RowVector3d normal = (b - a).cross(c - a).normalized();
        // e1 点在平面的投影点
        double distance1 = (e1 - a).dot(normal);
        Eigen::RowVector3d p1 = e1 - distance1 * normal;

        double distance2 = (e2 - a).dot(normal);
        Eigen::RowVector3d p2 = e2 - distance2 * normal;

        // std::vector<Eigen::RowVector3d> intersection_points;
        Tool_intersectLineWithTriangle(
            e1,e2,a,b,c,
            project_edge
        );
}

double project_SegmentOntoSegment(
    const Eigen::RowVector3d& a,
    const Eigen::RowVector3d& b,
    const Eigen::RowVector3d& c,
    const Eigen::RowVector3d& d,
    Eigen::RowVector3d& closestPoint1,
    Eigen::RowVector3d& closestPoint2){  // cd 上的最近点
    
        Eigen::RowVector3d e1 = b - a;
        Eigen::RowVector3d e2 = d - c;
        Eigen::RowVector3d r = a - c;

        double e1_e1 = e1.dot(e1);
        double e1_e2 = e1.dot(e2);
        double e2_e2 = e2.dot(e2);
        double e1_r = e1.dot(r);
        double e2_r = e2.dot(r);

        double denom = e1_e1 * e2_e2 - e1_e2 * e1_e2;
        double t,s;
        if(denom != 0){
            t = (e1_e2 * e2_r - e2_e2 * e1_r) / denom;
            s = (e1_e1 * e2_r - e1_e2 * e1_r) / denom;
        }
        else{
            // 平行 
            t = 0.0;
            s = e2_r / e2_e2;
            
        }

        t = std::max(0.0,std::min(1.0,t));
        s = std::max(0.0,std::min(1.0,s));

        closestPoint1 = a + t * e1;
        closestPoint2 = c + s * e2;
        
        return (closestPoint1 - closestPoint2).norm();  // 最短距离
}

// 判断投影点是否在三角形内部
bool isPointInTriangle(const Eigen::RowVector3d& p, const Eigen::RowVector3d& a,
                       const Eigen::RowVector3d& b, const Eigen::RowVector3d& c) {

    Eigen::RowVector3d v0 = c - a;
    Eigen::RowVector3d v1 = b - a;
    Eigen::RowVector3d v2 = p - a;

    double d00 = v0.dot(v0);
    double d01 = v0.dot(v1);
    double d11 = v1.dot(v1);
    double d20 = v2.dot(v0);
    double d21 = v2.dot(v1);

    double denom = d00 * d11 - d01 * d01;

    double v = (d11 * d20 - d01 * d21) / denom;
    double w = (d00 * d21 - d01 * d20) / denom;
    double u = 1.0 - v - w;

    // return (u > 0) && (v > 0) && (w > 0);  // 内部
    return (u >= 0) && (v >= 0) && (w >= 0);  // 包括边缘
}

void project_SegmentOntoSurface(
    const Eigen::RowVector3d& e1,
    const Eigen::RowVector3d& e2,
    const SurfaceMesh& surf_in,
    const double& eps,
    std::vector<Eigen::RowVector3d>& pro_points){

    double dis1 = -1; 
    double dis2 = -1; 
    Eigen::RowVector3d e1_;
    Eigen::RowVector3d e2_;
    int tri_id1 = -1;
    int tri_id2 = -1;

    for(int i=0;i<surf_in.tris.size();i++){
        Eigen::RowVector3d a = Eigen::RowVector3d(surf_in.coords[surf_in.tris[i][0]][0],
                                surf_in.coords[surf_in.tris[i][0]][1],
                                surf_in.coords[surf_in.tris[i][0]][2]);
        Eigen::RowVector3d b = Eigen::RowVector3d(surf_in.coords[surf_in.tris[i][1]][0],
                                surf_in.coords[surf_in.tris[i][1]][1],
                                surf_in.coords[surf_in.tris[i][1]][2]);
        Eigen::RowVector3d c = Eigen::RowVector3d(surf_in.coords[surf_in.tris[i][2]][0],
                                surf_in.coords[surf_in.tris[i][2]][1],
                                surf_in.coords[surf_in.tris[i][2]][2]);

        Eigen::RowVector3d p1;
        double distance1;
        TiGER::repair::project_PointOntoTriangle(
            e1,a,b,c,p1,distance1);
        
        if(fabs(distance1) <= eps && (dis1 == -1 || fabs(distance1) < dis1)){
            dis1 = distance1;
            e1_ = p1;
            tri_id1 = i;
        }

        Eigen::RowVector3d p2;
        double distance2;

        TiGER::repair::project_PointOntoTriangle(
            e2,a,b,c,p2,distance2);
        if(fabs(distance2) <= eps && (dis2 == -1 || fabs(distance2) < dis2)){
            dis2 = distance2;
            e2_ = p2;
            tri_id2 = i;
        }
    }

    // e1 投影点存在
    if(dis1 != -1 && dis1 <= eps){
        pro_points.push_back(e1_);
    }

    std::map<std::array<int,2>,std::vector<int>> edge_tris;  // 边的邻接三角形
    for(int i=0;i<surf_in.tris.size();i++){
        int p1 = surf_in.tris[i][0];
        int p2 = surf_in.tris[i][1];
        int p3 = surf_in.tris[i][2];
        std::array<int,2> edge1 = {p1,p2};
        std::array<int,2> edge2 = {p1,p3};
        std::array<int,2> edge3 = {p2,p3};
        std::sort(std::begin(edge1), std::end(edge1));
        std::sort(std::begin(edge2), std::end(edge2));
        std::sort(std::begin(edge3), std::end(edge3));
        edge_tris[edge1].push_back(i);
        edge_tris[edge2].push_back(i);
        edge_tris[edge3].push_back(i);
    }

    // 找出两段之间的三角形
    std::queue<int> q;
    q.push(tri_id1);
    while(!q.empty()){
        int index = q.front();
        q.pop();

        if(index == tri_id2)
            break;

        int p1 = surf_in.tris[index][0];
        int p2 = surf_in.tris[index][1];
        int p3 = surf_in.tris[index][2];
        Eigen::RowVector3d x = Eigen::RowVector3d(
                (surf_in.coords[p1][0] + surf_in.coords[p2][0] + surf_in.coords[p3][0]) / 3,
                (surf_in.coords[p1][1] + surf_in.coords[p2][1] + surf_in.coords[p3][1]) / 3,
                (surf_in.coords[p1][2] + surf_in.coords[p2][2] + surf_in.coords[p3][2]) / 3);
        std::array<int,2> edge1 = {p1,p2};
        std::array<int,2> edge2 = {p1,p3};
        std::array<int,2> edge3 = {p2,p3};
        std::sort(std::begin(edge1), std::end(edge1));
        std::sort(std::begin(edge2), std::end(edge2));
        std::sort(std::begin(edge3), std::end(edge3));
        double va = -1;
        int next_id = -1;
        int edge_id = -1;
        // 当前三角形与邻接三角形重心方向 点乘 两投影点方向 e1_e2_ 为正
        for(int i=0;i<edge_tris[edge1].size();i++){
            if(edge_tris[edge1][i] != index){
                int id = edge_tris[edge1][i];
                int q1 = surf_in.tris[id][0];
                int q2 = surf_in.tris[id][1];
                int q3 = surf_in.tris[id][2];
                Eigen::RowVector3d y = Eigen::RowVector3d(
                    (surf_in.coords[q1][0] + surf_in.coords[q2][0] + surf_in.coords[q3][0]) / 3,
                    (surf_in.coords[q1][1] + surf_in.coords[q2][1] + surf_in.coords[q3][1]) / 3,
                    (surf_in.coords[q1][2] + surf_in.coords[q2][2] + surf_in.coords[q3][2]) / 3);
                double d = (y - x).normalized().dot(e2_ - e1_);
                if(d > va){
                    va = d;
                    next_id = id;
                    edge_id = 1;
                }
            }
        }
        for(int i=0;i<edge_tris[edge2].size();i++){
            if(edge_tris[edge2][i] != index){
                int id = edge_tris[edge2][i];
                int q1 = surf_in.tris[id][0];
                int q2 = surf_in.tris[id][1];
                int q3 = surf_in.tris[id][2];
                Eigen::RowVector3d y = Eigen::RowVector3d(
                    (surf_in.coords[q1][0] + surf_in.coords[q2][0] + surf_in.coords[q3][0]) / 3,
                    (surf_in.coords[q1][1] + surf_in.coords[q2][1] + surf_in.coords[q3][1]) / 3,
                    (surf_in.coords[q1][2] + surf_in.coords[q2][2] + surf_in.coords[q3][2]) / 3);
                double d = (y - x).normalized().dot(e2_ - e1_);
                if(d > va){
                    va = d;
                    next_id = id;
                    edge_id = 2;
                }
            }
        }
        for(int i=0;i<edge_tris[edge3].size();i++){
            if(edge_tris[edge3][i] != index){
                int id = edge_tris[edge3][i];
                int q1 = surf_in.tris[id][0];
                int q2 = surf_in.tris[id][1];
                int q3 = surf_in.tris[id][2];
                Eigen::RowVector3d y = Eigen::RowVector3d(
                    (surf_in.coords[q1][0] + surf_in.coords[q2][0] + surf_in.coords[q3][0]) / 3,
                    (surf_in.coords[q1][1] + surf_in.coords[q2][1] + surf_in.coords[q3][1]) / 3,
                    (surf_in.coords[q1][2] + surf_in.coords[q2][2] + surf_in.coords[q3][2]) / 3);
                double d = (y - x).normalized().dot(e2_ - e1_);
                if(d > va){
                    va = d;
                    next_id = id;
                    edge_id = 3;
                }
            }
        }
        q.push(next_id);

        Eigen::RowVector3d a;
        Eigen::RowVector3d b;
        if(edge_id == 1){
            a = Eigen::RowVector3d(
                    surf_in.coords[p1][0],
                    surf_in.coords[p1][1],
                    surf_in.coords[p1][2]);
            b = Eigen::RowVector3d(
                    surf_in.coords[p2][0],
                    surf_in.coords[p2][1],
                    surf_in.coords[p2][2]); 
        }
        else if(edge_id == 3){
            a = Eigen::RowVector3d(
                    surf_in.coords[p2][0],
                    surf_in.coords[p2][1],
                    surf_in.coords[p2][2]);
            b = Eigen::RowVector3d(
                    surf_in.coords[p3][0],
                    surf_in.coords[p3][1],
                    surf_in.coords[p3][2]); 
        }
        else if(edge_id == 2){
            a = Eigen::RowVector3d(
                    surf_in.coords[p3][0],
                    surf_in.coords[p3][1],
                    surf_in.coords[p3][2]);
            b = Eigen::RowVector3d(
                    surf_in.coords[p1][0],
                    surf_in.coords[p1][1],
                    surf_in.coords[p1][2]); 
        }

        Eigen::RowVector3d closestPoint1;
        Eigen::RowVector3d closestPoint2;
        double dis = TiGER::repair::project_SegmentOntoSegment(e1,e2,a,b,closestPoint1,closestPoint2);
        pro_points.push_back(closestPoint2);
    }
    if(dis2 != -1 && dis2 <= eps){
        pro_points.push_back(e2_);
    }
}

void edge_repair(
    const SurfaceMesh& surf_in,
    const std::array<int,2>& edge_in,
    const std::vector<Eigen::RowVector3d>& points,
    SurfaceMesh& surf_out){

    int e1 = edge_in[0];
    int e2 = edge_in[1];

    surf_out.coords = surf_in.coords;
    surf_out.tris = surf_in.tris;

    surf_out.coords[e1] = {points[0].x(),points[0].y(),points[0].z()};
    
    surf_out.coords[e2] = {points.back().x(),points.back().y(),points.back().z()};

    if(points.size() == 2){ // 没有插点
        return;
    }

    std::vector<int> p_index(points.size());
    p_index[0] = e1;
    p_index.back() = e2;
    for(int i=1;i<points.size()-1;i++){
        // std::array<double,3> p = {}
        if(points.size() == 3)
            std::cout<<points[i].x()<<" -- "<<points[i].y()<<" -- "<<points[i].z()<<std::endl;
        surf_out.coords.push_back({points[i].x(),points[i].y(),points[i].z()});
        p_index[i] = surf_out.coords.size() - 1;
    }

    // surf_out.v_conn_points.resize(surf_out.coords.size());
    surf_out.v_conn_points = surf_in.v_conn_points;
    // surf_out.v_conn_tris.resize(surf_out.coords.size());
    surf_out.v_conn_tris = surf_in.v_conn_tris;

    for(int i=0;i<points.size()-2;i++){
        std::unordered_set<int> s;
        surf_out.v_conn_points.push_back(s);
        surf_out.v_conn_tris.push_back(s);
    }

    surf_out.v_conn_points[e1].erase(e2);
    surf_out.v_conn_points[e2].erase(e1);

    for(int u : surf_in.v_conn_tris[e1]){
        for(int v : surf_in.v_conn_tris[e2]){
            if(u == v){ // e1e2边的三角形
                int index = -1;
                for(int i=0;i<3;i++){
                    if(surf_in.tris[u][i] != e1 &&surf_in.tris[u][i] != e2){
                        index = surf_in.tris[u][i];
                    }
                }
                // surf_out.v_conn_tris[e1].erase(u);
                surf_out.v_conn_tris[e2].erase(u);
                surf_out.tris[u] = {index,e1,p_index[1]}; // 第一个就原地改
                surf_out.v_conn_tris[p_index[1]].insert(u);
                for(int i=2;i<p_index.size();i++){
                    int id = surf_out.tris.size();
                    surf_out.tris.push_back({index,p_index[i-1],p_index[i]}); // 
                    surf_out.v_conn_tris[index].insert(id);
                    surf_out.v_conn_tris[p_index[i-1]].insert(id);
                    surf_out.v_conn_tris[p_index[i]].insert(id);
                }
            }
        }
    }

}

void partition(
    const Eigen::RowVector3d& a,
    const Eigen::RowVector3d& b,
    const Eigen::RowVector3d& c,
    const std::vector<Eigen::RowVector3d>& project_poins,
    const std::array<int,2>& vertextid,
    const std::array<int,2>& edgeid,
    TiGER::SurfaceMesh& surf_out){
        surf_out.coords.push_back({a(0),a(1),a(2)});
        surf_out.coords.push_back({b(0),b(1),b(2)});
        surf_out.coords.push_back({c(0),c(1),c(2)});

        int point_size = project_poins.size();
        if(point_size == 1){
            if(vertextid[0] != -1){ // 唯一一个点交于三角形端点
                surf_out.tris.push_back({0,1,2});
                return;
            }
            else{
                surf_out.coords.push_back({project_poins[0](0),project_poins[0](1),project_poins[0](2)});
                if(edgeid[0] == -1){ 
                    for(int i=0;i<3;i++){
                        surf_out.tris.push_back({i,(i+1)%3,3});
                    }
                    return;
                }
                else{ // 交于第i边  0 ab 1 bc 2 ca
                    surf_out.tris.push_back({edgeid[0],(edgeid[0]+2)%3,3});
                    surf_out.tris.push_back({(edgeid[0]+1)%3,(edgeid[0]+2)%3,3});
                    return;
                }
            }
        }
        else{  // 两个交点
            if(vertextid[0] != -1 && vertextid[1] != -1){  // 共边不切
                surf_out.tris.push_back({0,1,2});
                return;
            }
            else if(vertextid[0] == -1 && vertextid[1] == -1){  //与端点不相交
                surf_out.coords.push_back({project_poins[0](0),project_poins[0](1),project_poins[0](2)});
                surf_out.coords.push_back({project_poins[1](0),project_poins[1](1),project_poins[1](2)});
                if(edgeid[0] == -1 && edgeid[1] == -1){
                    // 在内部 
                    Eigen::RowVector3d project_point1;
                    double dis1 = 0;
                    project_PointOntoSegment(project_poins[0],a,b,project_point1,dis1);
                    Eigen::RowVector3d project_point2;
                    double dis2 = 0;
                    project_PointOntoSegment(project_poins[1],a,b,project_point2,dis2);

                    if(dis1 <= dis2){
                        surf_out.tris.push_back({0,1,3});
                        Eigen::RowVector3d project_point3;
                        double dis3 = 0;
                        project_PointOntoSegment(project_poins[0],b,c,project_point3,dis3);
                        Eigen::RowVector3d project_point4;
                        double dis4 = 0;
                        project_PointOntoSegment(project_poins[1],b,c,project_point4,dis4);

                        if(dis3 <= dis4){
                            surf_out.tris.push_back({1,2,3});
                            surf_out.tris.push_back({2,0,4});
                            surf_out.tris.push_back({3,4,0});
                            surf_out.tris.push_back({3,4,2});
                        }
                        else{
                            surf_out.tris.push_back({1,2,4});

                            Eigen::RowVector3d project_point5;
                            double dis5 = 0;
                            project_PointOntoSegment(project_poins[0],c,a,project_point5,dis5);
                            Eigen::RowVector3d project_point6;
                            double dis6 = 0;
                            project_PointOntoSegment(project_poins[1],c,a,project_point6,dis6);

                            if(dis5 <= dis6){
                                surf_out.tris.push_back({2,0,3});
                                surf_out.tris.push_back({3,4,1});
                                surf_out.tris.push_back({3,4,2});
                            }
                            else{
                                surf_out.tris.push_back({2,0,4});
                                surf_out.tris.push_back({3,4,1});
                                surf_out.tris.push_back({3,4,0});
                            }
                        }
                    }
                    else{
                        surf_out.tris.push_back({0,1,4});
                        Eigen::RowVector3d project_point3;
                        double dis3 = 0;
                        project_PointOntoSegment(project_poins[0],b,c,project_point3,dis3);
                        Eigen::RowVector3d project_point4;
                        double dis4 = 0;
                        project_PointOntoSegment(project_poins[1],b,c,project_point4,dis4);

                        if(dis3 <= dis4){
                            surf_out.tris.push_back({1,2,3});
                            Eigen::RowVector3d project_point5;
                            double dis5 = 0;
                            project_PointOntoSegment(project_poins[0],c,a,project_point5,dis5);
                            Eigen::RowVector3d project_point6;
                            double dis6 = 0;
                            project_PointOntoSegment(project_poins[1],c,a,project_point6,dis6);
                            if(dis5 <= dis6){
                                surf_out.tris.push_back({2,0,3});
                                surf_out.tris.push_back({3,4,0});
                                surf_out.tris.push_back({3,4,1});
                            }
                            else{
                                surf_out.tris.push_back({2,0,4});
                                surf_out.tris.push_back({3,4,1});
                                surf_out.tris.push_back({3,4,2});
                            }
                        }
                        else{
                            surf_out.tris.push_back({1,2,4});
                            surf_out.tris.push_back({2,0,3});
                            surf_out.tris.push_back({3,4,0});
                            surf_out.tris.push_back({3,4,2});
                        }
                    }
                }
                else if(edgeid[0] != -1 && edgeid[1] != -1){
                    // 全在边上
                    if(edgeid[0] == edgeid[1]){
                        surf_out.tris.push_back({3,4,(edgeid[0]+2)%3});
                        //同边
                        if(edgeid[0] == 0){  // ab
                            Eigen::RowVector3d project_point1;
                            double dis1 = 0;
                            project_PointOntoSegment(project_poins[0],b,c,project_point1,dis1);

                            Eigen::RowVector3d project_point2;
                            double dis2 = 0;
                            project_PointOntoSegment(project_poins[1],b,c,project_point2,dis2);
                            if(dis1 <= dis2){
                                surf_out.tris.push_back({1,2,3});
                                surf_out.tris.push_back({0,2,4});
                            }
                            else{
                                surf_out.tris.push_back({1,2,4});
                                surf_out.tris.push_back({0,2,3});
                            }
                        }
                        else if(edgeid[0] == 1){  // bc
                            Eigen::RowVector3d project_point1;
                            double dis1 = 0;
                            project_PointOntoSegment(project_poins[0],c,a,project_point1,dis1);

                            Eigen::RowVector3d project_point2;
                            double dis2 = 0;
                            project_PointOntoSegment(project_poins[1],c,a,project_point2,dis2);
                            if(dis1 <= dis2){
                                surf_out.tris.push_back({2,0,3});
                                surf_out.tris.push_back({0,1,4});
                            }
                            else{
                                surf_out.tris.push_back({2,0,4});
                                surf_out.tris.push_back({0,1,3});
                            }
                        }
                        else{  // ca
                            Eigen::RowVector3d project_point1;
                            double dis1 = 0;
                            project_PointOntoSegment(project_poins[0],a,b,project_point1,dis1);

                            Eigen::RowVector3d project_point2;
                            double dis2 = 0;
                            project_PointOntoSegment(project_poins[1],a,b,project_point2,dis2);
                            if(dis1 <= dis2){
                                surf_out.tris.push_back({1,0,3});
                                surf_out.tris.push_back({1,2,4});
                            }
                            else{
                                surf_out.tris.push_back({1,0,4});
                                surf_out.tris.push_back({1,2,3});
                            }
                        }
                    }
                    else{
                        surf_out.tris.push_back({3,4,edgeid[0]});
                        surf_out.tris.push_back({3,4,(edgeid[0]+1)%3});
                        if((edgeid[0] + 1)%3 == edgeid[1]){
                            surf_out.tris.push_back({edgeid[0],(edgeid[1]+1)%3,4});
                        }
                        else{
                            surf_out.tris.push_back({(edgeid[0]+1)%3,(edgeid[0]+2)%3,4});
                        }                        
                    }
                }
                else{
                    // 一个在边
                    int p = edgeid[0] == -1 ? 1 : 0;  // 找到与边相交的点
                    surf_out.tris.push_back({(edgeid[p]+1)%3,(edgeid[p]+2)%3,1-p+3});
                    surf_out.tris.push_back({edgeid[p],(edgeid[p]+2)%3,1-p+3});

                    surf_out.tris.push_back({edgeid[p],3,4});
                    surf_out.tris.push_back({(edgeid[p]+1)%3,3,4});
                }
            }
            else{  // 与一个端点相交
                // int p = max(vertextid[0],vertextid[1]);
                int p = vertextid[0] == -1 ? 1 : 0;
                int e = std::max(edgeid[0],edgeid[1]);
                surf_out.coords.push_back({project_poins[1-p](0),project_poins[1-p](1),project_poins[1-p](2)});
                if(e == -1){  // 一点在内部
                    // surf_out.coords.push_back({project_poins[1-p](0),project_poins[1-p](1),project_poins[1-p](2)});
                    for(int i=0;i<3;i++){
                        surf_out.tris.push_back({i,(i+1)%3,3});
                    }
                    return;
                }
                else{  // 有一点在边上
                    surf_out.tris.push_back({e,(e+2)%3,3});
                    surf_out.tris.push_back({(e+1)%3,(e+2)%3,3});
                    return;
                }
            }
        }
    }

enum intersectionresult Tool_intersectLines(
    const Eigen::RowVector3d& e1,
    const Eigen::RowVector3d& e2,
    const Eigen::RowVector3d& e3,
    const Eigen::RowVector3d& e4,
    std::vector<Eigen::RowVector3d>& intersection_points){

        // 四点是否共面  两线段是否共面
        REAL* real_e1 = const_cast<REAL*>(e1.data());
        REAL* real_e2 = const_cast<REAL*>(e2.data());
        REAL* real_e3 = const_cast<REAL*>(e3.data());
        REAL* real_e4 = const_cast<REAL*>(e4.data());
        double d = TiGER_GEOM_FUNC::orient3d(real_e1,real_e2,real_e3,real_e4);
        if(std::fabs(d) < 1e-10){  // 共面
            Eigen::RowVector3d normal = (e2 - e1).cross(e4 - e3).normalized();  // 单位法向
            std::vector<Eigen::RowVector3d> axis_normal = {Eigen::RowVector3d(1,0,0),
                                            Eigen::RowVector3d(0,1,0),
                                            Eigen::RowVector3d(0,0,1)};  // 坐标轴
            int axis = 0;
            for(int i=0;i<3;i++){ // 选择与法向最近的坐标轴
                if(fabs(normal.dot(axis_normal[i])) > 0.5) {
                    axis = i;
                    break;
                }
            }
            // 确定投影坐标
            int x = (axis + 1) % 3;
            int y = (axis + 2) % 3;
            // 投影到二维平面
            Eigen::Vector2d e1_(e1[x],e1[y]);
            Eigen::Vector2d e2_(e2[x],e2[y]);
            Eigen::Vector2d e3_(e3[x],e3[y]);
            Eigen::Vector2d e4_(e4[x],e4[y]);
            std::vector<Eigen::Vector2d> points;
            enum intersectionresult result = Tool_intersectLines2D(e1_,e2_,e3_,e4_,points);
            if(result == SHAREEDGE || result == DISJOINT || result == SHAREVERTEX){
                for(auto& point : points){
                    double t = (point - e1_).norm() / (e2_ - e1_).norm();
                    intersection_points.push_back(e1 + t * (e2 - e1));
                }
                return result;
            }
            else { // 相交于一点
                double t = (points[0] - e1_).norm() / (e2_ - e1_).norm();
                intersection_points.push_back(e1 + t * (e2 - e1));
                return INTERSECT;
            }
        }
        else   // 不共面 
            return DISJOINT;
    }

enum intersectionresult Tool_intersectLines2D(
    const Eigen::Vector2d& e1,
    const Eigen::Vector2d& e2,
    const Eigen::Vector2d& e3,
    const Eigen::Vector2d& e4,
    std::vector<Eigen::Vector2d>& intersection_points){

        REAL* real_e1 = const_cast<REAL*>(e1.data());
        REAL* real_e2 = const_cast<REAL*>(e2.data());
        REAL* real_e3 = const_cast<REAL*>(e3.data());
        REAL* real_e4 = const_cast<REAL*>(e4.data());
        double d1 = TiGER_GEOM_FUNC::orient2d(real_e3, real_e4, real_e1);
        double d2 = TiGER_GEOM_FUNC::orient2d(real_e3, real_e4, real_e2);
        double d3 = TiGER_GEOM_FUNC::orient2d(real_e1, real_e2, real_e3);
        double d4 = TiGER_GEOM_FUNC::orient2d(real_e1, real_e2, real_e4);
        // 相交
        if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) &&
            ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))) {
            // e1e3e4 和 e2e3e4 面积比
            double w1 = std::fabs(d2) / (std::fabs(d1) + std::fabs(d2));
            double w2 = 1.0 - w1;
            intersection_points.push_back(w1 * e1 + w2 * e2);
            return INTERSECT;  // 相交于一点
        }

        auto onSegment = [](const Eigen::Vector2d& pi, const Eigen::Vector2d& pj,
                      const Eigen::Vector2d& pk) {
            // if ((pk - pi).norm() < 1e-10 || (pk - pj).norm() < 1e-10) return false;
            return ((std::min(pi.x(), pj.x()) <= pk.x() &&
                    pk.x() <= std::max(pi.x(), pj.x()) &&
                    std::min(pi.y(), pj.y()) <= pk.y() &&
                    pk.y() <= std::max(pi.y(), pj.y())));
        };
        
        bool e1_on = d1 == 0 && onSegment(e3,e4,e1);
        bool e2_on = d2 == 0 && onSegment(e3,e4,e2);
        bool e3_on = d3 == 0 && onSegment(e1,e2,e3);
        bool e4_on = d4 == 0 && onSegment(e1,e2,e4);

        // 线段e1e2端点与e3e4相交
        if(e1_on && e2_on){
            intersection_points.push_back(e1);
            intersection_points.push_back(e2);
            return SHAREEDGE;
        }
        else if(e1_on){
            intersection_points.push_back(e1);
            if(e3_on){
                if(e1 == e3)
                    return SHAREVERTEX;
                else{
                    intersection_points.push_back(e3);
                    return SHAREEDGE;
                }                
            }
            if(e4_on){
                if(e1 == e4)
                    return SHAREVERTEX;
                else{
                    intersection_points.push_back(e4);
                    return SHAREEDGE;
                }                
            }
            return INTERSECT;
        }
        else if(e2_on){
            intersection_points.push_back(e2);
            if(e3_on){
                if(e2 == e3)
                    return SHAREVERTEX;
                else{
                    intersection_points.push_back(e3);
                    return SHAREEDGE;
                }
            }
            if(e4_on){
                if(e2 == e4)
                    return SHAREVERTEX;
                else{
                    intersection_points.push_back(e4);
                    return SHAREEDGE;
                }
            }
            return INTERSECT;
        }
        // 线段e3e4端点与e1e2相交
        if(e3_on && e4_on){
            intersection_points.push_back(e3);
            intersection_points.push_back(e4);
            return SHAREEDGE;
        }
        else if(e3_on){
            intersection_points.push_back(e3);
            return INTERSECT;
        }
        else if(e4_on){
            intersection_points.push_back(e4);
            return INTERSECT;
        }
        
        return DISJOINT;
    }

enum intersectionresult Tool_intersectLineWithTriangle(
    const Eigen::RowVector3d& e1,
    const Eigen::RowVector3d& e2,
    const Eigen::RowVector3d& a,
    const Eigen::RowVector3d& b,
    const Eigen::RowVector3d& c,
    std::vector<Eigen::RowVector3d>& intersection_points){

        // GEOM_FUNC::exactinit();
        TiGER_GEOM_FUNC::exactinit_threadSafe();

        REAL* real_e1 = const_cast<REAL*>(e1.data());
        REAL* real_e2 = const_cast<REAL*>(e2.data());
        REAL* real_a = const_cast<REAL*>(a.data());
        REAL* real_b = const_cast<REAL*>(b.data());
        REAL* real_c = const_cast<REAL*>(c.data());
        double d1 = TiGER_GEOM_FUNC::orient3d(real_e1, real_a, real_b, real_c);
        double d2 = TiGER_GEOM_FUNC::orient3d(real_e2, real_a, real_b, real_c);

        // eps
        if (std::fabs(d1) < 1e-10) d1 = 0;
        if (std::fabs(d2) < 1e-10) d2 = 0;

        // std::cout<<"d1 "<<d1<<" "<<d2<<std::endl;

        if(d1 * d2 > 0){  // 在同一侧
            return DISJOINT;
        }
        if(d1 == 0 && d2 != 0){  // e1 与 三角形共面
            intersectionresult result = Tool_intersectPointWithTriangle(e1,a,b,c);
            if(result == DISJOINT)
                return DISJOINT;
            intersection_points.push_back(e1);
            if(result == NOD_INTERSECT_NOD)
                return LTI_INTERSECT_NOD;
            if(result == NOD_INTERSECT_EDG)
                return LTI_INTERSECT_EDG;
            if(result == NOD_INTERSECT_FAC)
                return LTI_INTERSECT_FAC;
        }
        if(d1 != 0 && d2 == 0){  // e2 与三角形共面
            intersectionresult result = Tool_intersectPointWithTriangle(e2,a,b,c);
            if(result == DISJOINT)
                return DISJOINT;
            intersection_points.push_back(e2);
            if(result == NOD_INTERSECT_NOD)
                return LTI_INTERSECT_NOD;
            if(result == NOD_INTERSECT_EDG)
                return LTI_INTERSECT_EDG;
            if(result == NOD_INTERSECT_FAC)
                return LTI_INTERSECT_FAC;
        }
        if(d1 == 0 && d2 == 0){  // e1e2 与三角形共面
            // std::cout<<"same "<<std::endl;
            // std::cout<<" "<<a(0)<<" "<<a(1)<<" "<<a(2)<<std::endl;
            // std::cout<<" "<<b(0)<<" "<<b(1)<<" "<<b(2)<<std::endl;
            // std::cout<<" "<<c(0)<<" "<<c(1)<<" "<<c(2)<<std::endl;
            // 将线和三角形投影到坐标平面
            Eigen::RowVector3d abc_normal = (b - a).cross(c - a).normalized();
            std::vector<Eigen::RowVector3d> axis_normal = {Eigen::RowVector3d(1,0,0),
                                                            Eigen::RowVector3d(0,1,0),
                                                            Eigen::RowVector3d(0,0,1)};
            int axis = 0;
            for(int i=0;i<3;i++){
                if(fabs(abc_normal.dot(axis_normal[i])) > 0.5){
                    axis = i;
                    break;
                }
            }
            int x = (axis + 1) % 3;
            int y = (axis + 2) % 3;
            Eigen::Vector2d e1_(e1[x], e1[y]);
            Eigen::Vector2d e2_(e2[x], e2[y]);
            Eigen::Vector2d a_(a[x], a[y]);
            Eigen::Vector2d b_(b[x], b[y]);
            Eigen::Vector2d c_(c[x], c[y]);

            // std::cout<<"e1- "<<e1[x]<<","<<e1[y]<<std::endl;
            // std::cout<<"e2- "<<e2[x]<<","<<e2[y]<<std::endl;
            // std::cout<<"a- "<<a[x]<<","<<a[y]<<std::endl;
            // std::cout<<"b- "<<b[x]<<","<<b[y]<<std::endl;
            // std::cout<<"c- "<<c[x]<<","<<c[y]<<std::endl;

            std::vector<Eigen::Vector2d> points;
            enum intersectionresult result = Tool_intersectLineWithTriangle2D(e1_,e2_,a_,b_,c_,points);
            if(result == DISJOINT)
                return DISJOINT;
            for(auto& point : points){
                // std::cout<<" point ("<<point.x()<<","<<point.y()<<")"<<std::endl;
                double t = (point - e1_).norm() / (e2_ - e1_).norm();
                // std::cout<<t<<" "<<"point ("<<(e1 + t * (e2 - e1)).x()<<","<<(e1 + t * (e2 - e1)).y()<<","<<(e1 + t * (e2 - e1)).z()<<")"<<std::endl;
                intersection_points.push_back(e1 + t * (e2 - e1));
            }
            return result;
        }
        // 不共面
        else{
            //  两端点不在三角面上  相交在边  端点  面内
            double ab_line = TiGER_GEOM_FUNC::orient3d(real_a, real_b, real_e1, real_e2);
            double bc_line = TiGER_GEOM_FUNC::orient3d(real_b, real_c, real_e1, real_e2);
            double ac_line = TiGER_GEOM_FUNC::orient3d(real_c, real_a, real_e1, real_e2);
            
            // 交点在点上
            if(ab_line == 0.0 && bc_line == 0.0){
                intersection_points.push_back(b);
                return LTI_INTERSECT_NOD;
            }
            else if(bc_line == 0.0 && ac_line == 0.0){
                intersection_points.push_back(c);
                return LTI_INTERSECT_NOD;
            }
            else if(ac_line == 0.0 && ab_line == 0.0){
                intersection_points.push_back(a);
                return LTI_INTERSECT_NOD;
            }
            enum intersectionresult result;
            // 交点在边上
            if(ab_line == 0.0 && bc_line*ac_line > 0.0){
                result = LTI_INTERSECT_EDG;
            }   
            else if(bc_line == 0.0 && ab_line*ac_line > 0.0){
                result = LTI_INTERSECT_EDG;
            }
            else if(ac_line == 0.0 && ab_line*bc_line > 0.0){
                result = LTI_INTERSECT_EDG;
            }
            // 交点在面内
            else if(ab_line > 0.0 && bc_line > 0.0 && ac_line > 0.0){
                result = LTI_INTERSECT_FAC;
            }
            else if(ab_line < 0.0 && bc_line < 0.0 && ac_line < 0.0){
                result = LTI_INTERSECT_FAC;
            }

            // 求交点
            Eigen::RowVector3d point;
            point(0) = TiGER_GEOM_FUNC::fixedSplitPoint(fabs(d1), fabs(d2), e1(0), e2(0));
            point(1) = TiGER_GEOM_FUNC::fixedSplitPoint(fabs(d1), fabs(d2), e1(1), e2(1));
            point(2) = TiGER_GEOM_FUNC::fixedSplitPoint(fabs(d1), fabs(d2), e1(2), e2(2));
            intersection_points.push_back(point);
            return result;
        }
    }

enum intersectionresult Tool_intersectLineWithTriangle2D(
    const Eigen::Vector2d& e1, 
    const Eigen::Vector2d& e2,
    const Eigen::Vector2d& a, 
    const Eigen::Vector2d& b,
    const Eigen::Vector2d& c,
    std::vector<Eigen::Vector2d>& intersection_points){

        REAL* real_e1 = const_cast<REAL*>(e1.data());
        REAL* real_e2 = const_cast<REAL*>(e2.data());
        REAL* real_a = const_cast<REAL*>(a.data());
        REAL* real_b = const_cast<REAL*>(b.data());
        REAL* real_c = const_cast<REAL*>(c.data());

        double d = TiGER_GEOM_FUNC::orient2d(real_c,real_a,real_b);
        double d1 = TiGER_GEOM_FUNC::orient2d(real_e1, real_a, real_b);
        double d2 = TiGER_GEOM_FUNC::orient2d(real_e1, real_b, real_c);
        double d3 = TiGER_GEOM_FUNC::orient2d(real_e1, real_c, real_a);
        double d4 = TiGER_GEOM_FUNC::orient2d(real_e2, real_a, real_b);
        double d5 = TiGER_GEOM_FUNC::orient2d(real_e2, real_b, real_c);
        double d6 = TiGER_GEOM_FUNC::orient2d(real_e2, real_c, real_a);

        // eps
        if (std::fabs(d1) < 1e-10) d1 = 0;
        if (std::fabs(d2) < 1e-10) d2 = 0;
        if (std::fabs(d3) < 1e-10) d3 = 0;
        if (std::fabs(d4) < 1e-10) d4 = 0;
        if (std::fabs(d5) < 1e-10) d5 = 0;
        if (std::fabs(d6) < 1e-10) d6 = 0;
        // std::cout<<d<<" "<<d1<<" "<<d2<<" "<<d3<<" "<<d4<<" "<<d5<<" "<<d6<<std::endl;

        bool isInner1 = d*d1 >= 0.0 && d*d2 >= 0.0 && d*d3 >= 0.0;
        bool isInner2 = d*d4 >= 0.0 && d*d5 >= 0.0 && d*d6 >= 0.0;

        if(isInner1 && isInner2){  // 直线两个点都在平面的内部或边界上时，这条直线被包含在平面内部
            // std::cout<<"two "<<std::endl;
            intersection_points.push_back(e1);
            intersection_points.push_back(e2);
            return LTI_INTERSECT_INS;
        }
        else if(isInner1){  // 直线1个点在平面的内部或边界上时，这条直线和平面必有一个交点
            intersection_points.push_back(e1);
            enum intersectionresult result = DISJOINT;
            // 直线1个点在边界
            // 有交点在三角形端点上
            int cod = 0;
            if(d1 == 0.0 && d2 == 0.0){
                cod = 1;  // b 
                result = LTI_INTERSECT_NOD;
            }
            else if(d2 == 0.0 && d3 == 0.0){
                cod = 2;  // c
                result = LTI_INTERSECT_NOD;
            }
            else if(d3 == 0.0 && d1 == 0.0){
                cod = 0;  // a
                result = LTI_INTERSECT_NOD;
            }
            
            if(result == DISJOINT){
                // 有交点在边上
                if(d1 == 0.0){
                    cod = 0;  // ab
                    result = LTI_INTERSECT_EDG;
                }
                if(d2 == 0.0){
                    cod = 1;  // bc
                    result = LTI_INTERSECT_EDG;
                }
                if(d3 == 0.0){
                    // std::cout<<" --ab-- "<<std::endl;
                    cod = 2;  //ac
                    result = LTI_INTERSECT_EDG;
                }
            }

            // 直线1个点在内部
            if(result == DISJOINT){
                double a_e1e2 = TiGER_GEOM_FUNC::orient2d(real_e1,real_e2,real_a);
                double b_e1e2 = TiGER_GEOM_FUNC::orient2d(real_e1,real_e2,real_b);
                double c_e1e2 = TiGER_GEOM_FUNC::orient2d(real_e1,real_e2,real_c);
                // 端点在直线上
                if(a_e1e2 == 0.0){
                    intersection_points.push_back(a);
                    result = LTI_INTERSECT_NOD;
                }
                else if(b_e1e2 == 0.0){
                    intersection_points.push_back(b);
                    result = LTI_INTERSECT_NOD;
                }
                else if(c_e1e2 == 0.0){
                    intersection_points.push_back(c);
                    result = LTI_INTERSECT_NOD;
                }
                // 有一条边与直线相交
                if(result == DISJOINT){
                    std::vector<Eigen::Vector2d> points;
                    if(Tool_intersectLines2D(e1,e2,a,b,points) != DISJOINT){
                        intersection_points.push_back(points[0]);
                    }
                    else if(Tool_intersectLines2D(e1,e2,b,c,points) != DISJOINT){
                        intersection_points.push_back(points[0]);
                    }
                    else if(Tool_intersectLines2D(e1,e2,c,a,points) != DISJOINT){
                        intersection_points.push_back(points[0]);
                    }
                    result = LTI_INTERSECT_EDG;
                }
            }
            else{
                // 直线与三角形有第二个交点
                if(result == LTI_INTERSECT_NOD){
                    if(cod == 0){  // a
                        if(d4*d >= 0.0 && d6*d >= 0.0){
                            if(d4 == 0.0){  // ab 共线
                                intersection_points.push_back(b);
                                return LTI_INTERSECT_EDG;
                            }
                            else if(d6 == 0){  // ac 共线
                                intersection_points.push_back(c);
                                return LTI_INTERSECT_EDG;
                            }
                            else{  //交点在对边
                                std::vector<Eigen::Vector2d> points;
                                Tool_intersectLines2D(e1,e2,b,c,points);
                                intersection_points.push_back(points[0]);
                            }
                        }
                    }
                    else if(cod == 1){
                        if(d4*d >= 0.0 && d5*d >= 0.0){
                            if(d4 == 0.0){
                                intersection_points.push_back(a);
                                return LTI_INTERSECT_NOD_EDG;
                            }
                            else if(d5 == 0.0){
                                intersection_points.push_back(c);
                                return LTI_INTERSECT_EDG;
                            }
                            else{
                                std::vector<Eigen::Vector2d> points;
                                Tool_intersectLines2D(e1,e2,a,c,points);
                                intersection_points.push_back(points[0]);
                            }
                        }
                    }
                    else if(cod == 2){
                        if(d5*d >= 0.0 && d6*d >= 0.0){
                            if(d5 == 0.0){
                                intersection_points.push_back(b);
                                return LTI_INTERSECT_EDG;
                            }
                            else if(d6 == 0.0){
                                intersection_points.push_back(a);
                                return LTI_INTERSECT_EDG;
                            }
                            else{
                                std::vector<Eigen::Vector2d> points;
                                Tool_intersectLines2D(e1,e2,a,b,points);
                                intersection_points.push_back(points[0]);
                            }
                        }
                    }
                }
                else if(result == LTI_INTERSECT_EDG){
                    /* 可能和两位两条边相交 */
                    if(cod == 0){  // ab
                    /* 如果另外1个顶点在边的正侧，存在两个交点 */
                        if(d4*d > 0.0){
                            if(d5*d < 0.0 && d6*d < 0.0){
                                double b_e1e2 = TiGER_GEOM_FUNC::orient2d(real_e1,real_e2,real_b);
                                double c_e1e2 = TiGER_GEOM_FUNC::orient2d(real_e1,real_e2,real_c);
                                if(b_e1e2 == 0.0){
                                    // result = LTI_INTERSECT_NOD;
                                    intersection_points.push_back(b);
                                    return LTI_INTERSECT_NOD;
                                }
                                else if(c_e1e2 == 0.0){
                                    // result = LTI_INTERSECT_NOD;
                                    intersection_points.push_back(c);
                                    return LTI_INTERSECT_NOD;
                                }
                                else if(b_e1e2*c_e1e2 > 0.0){
                                    // ac
                                    cod = 2;
                                }
                            }
                            else if(d5*d < 0.0){
                                cod = 1;
                            }
                            else if(d6*d < 0.0){
                                cod = 2;
                            }

                            if(cod == 1){ // bc
                                std::vector<Eigen::Vector2d> points;
                                Tool_intersectLines2D(e1,e2,b,c,points);
                                return LTI_INTERSECT_EDG;
                            }
                            else if(cod == 2){
                                std::vector<Eigen::Vector2d> points;
                                Tool_intersectLines2D(e1,e2,a,c,points);
                                return LTI_INTERSECT_EDG;
                            }
                        }
                        else if(d4 == 0.0){ /* 交点在顶点上，此时，三角形的一条边是原线段的子线段 */
                            if(d5 < 0.0){
                                intersection_points.push_back(b);
                            }
                            else if(d6 < 0.0){
                                intersection_points.push_back(a);
                            }
                            return LTI_INTERSECT_NOD;
                        }
                    }
                    else if(cod == 1){  // bc
                        if(d5*d > 0.0){
                            if(d4*d < 0.0 && d6*d < 0.0){
                                double a_e1e2 = TiGER_GEOM_FUNC::orient2d(real_e1,real_e2,real_a);
                                double c_e1e2 = TiGER_GEOM_FUNC::orient2d(real_e1,real_e2,real_c);
                                if(a_e1e2 == 0.0){
                                    intersection_points.push_back(a);
                                    return LTI_INTERSECT_NOD;
                                }
                                else if(c_e1e2 == 0.0){
                                    intersection_points.push_back(c);
                                    return LTI_INTERSECT_NOD;
                                }
                                else if(a_e1e2*c_e1e2 > 0.0){
                                    cod = 0;
                                }
                            }
                            else if(d4*d < 0.0){
                                cod = 0;
                            }
                            else if(d6*d < 0.0){
                                cod = 2;
                            }

                            if(cod == 0){  // ab
                                std::vector<Eigen::Vector2d> points;
                                Tool_intersectLines2D(e1,e2,a,b,points);
                                return LTI_INTERSECT_EDG;
                            }
                            else if(cod == 2){  // ac
                                std::vector<Eigen::Vector2d> points;
                                Tool_intersectLines2D(e1,e2,a,c,points);
                                return LTI_INTERSECT_EDG;
                            }
                        }
                        else if(d5 == 0.0){
                            if(d4 < 0.0){
                                intersection_points.push_back(b);
                            }
                            else if(d6 < 0.0){
                                intersection_points.push_back(c);
                            }
                            return LTI_INTERSECT_NOD;
                        }
                    }
                    else if(cod == 2){  // ac
                        if(d6*d > 0.0){
                            if(d4*d < 0.0 && d5*d < 0.0){
                                double b_e1e2 = TiGER_GEOM_FUNC::orient2d(real_e1,real_e2,real_b);
                                double a_e1e2 = TiGER_GEOM_FUNC::orient2d(real_e1,real_e2,real_a);
                                if(b_e1e2 == 0.0){
                                    intersection_points.push_back(b);
                                    return LTI_INTERSECT_NOD;
                                }
                                else if(a_e1e2 == 0.0){
                                    intersection_points.push_back(a);
                                    return LTI_INTERSECT_NOD;
                                }
                                else if(b_e1e2*a_e1e2 > 0.0){
                                    cod = 1;
                                }
                            }
                            else if(d4*d < 0.0){
                                cod = 0;
                            }
                            else if(d5*d < 0.0){
                                cod = 1;
                            }

                            if(cod == 0){  //ab
                                std::vector<Eigen::Vector2d> points;
                                Tool_intersectLines2D(e1,e2,a,b,points);
                                return LTI_INTERSECT_EDG;
                            }
                            else if(cod == 1){
                                std::vector<Eigen::Vector2d> points;
                                Tool_intersectLines2D(e1,e2,b,c,points);
                                return LTI_INTERSECT_EDG;
                            }
                        }
                        else if(d6 == 0.0){
                            if(d4 < 0.0){
                                intersection_points.push_back(a);
                            }
                            else if(d5 < 0.0){
                                intersection_points.push_back(c);
                            }
                            return LTI_INTERSECT_NOD;
                        }
                    }
                }
            }
        }
        else if(isInner2){  // e2 在内部
            intersection_points.push_back(e2);
            enum intersectionresult result = DISJOINT;
            // 直线1个点在边界
            // 有交点在三角形端点上
            int cod = 0;
            if(d4 == 0.0 && d5 == 0.0){
                cod = 1;
                result = LTI_INTERSECT_NOD;
            }
            else if(d5 == 0.0 && d6 == 0.0){
                cod = 2;
                result = LTI_INTERSECT_NOD;
            }
            else if(d6 == 0.0 && d4 == 0.0){
                cod = 0;
                result = LTI_INTERSECT_NOD;
            }
            
            if(result == DISJOINT){
                // 有交点在边上
                if(d4 == 0.0){
                    cod = 0;
                    result = LTI_INTERSECT_EDG;
                }
                if(d5 == 0.0){
                    cod = 1;
                    result = LTI_INTERSECT_NOD_EDG;
                }
                if(d6 == 0.0){
                    cod = 2;
                    result = LTI_INTERSECT_EDG;
                }
            }

            // 直线1个点在内部
            if(result == DISJOINT){
                double a_e1e2 = TiGER_GEOM_FUNC::orient2d(real_e1,real_e2,real_a);
                double b_e1e2 = TiGER_GEOM_FUNC::orient2d(real_e1,real_e2,real_b);
                double c_e1e2 = TiGER_GEOM_FUNC::orient2d(real_e1,real_e2,real_c);
                // 端点在直线上
                if(a_e1e2 == 0.0){
                    intersection_points.push_back(a);
                    result = LTI_INTERSECT_NOD;
                }
                else if(b_e1e2 == 0.0){
                    intersection_points.push_back(b);
                    result = LTI_INTERSECT_NOD;
                }
                else if(c_e1e2 == 0.0){
                    intersection_points.push_back(c);
                    result = LTI_INTERSECT_NOD;
                }
                // 有一条边与直线相交
                if(result == DISJOINT){
                    std::vector<Eigen::Vector2d> points;
                    if(Tool_intersectLines2D(e1,e2,a,b,points) != DISJOINT){
                        intersection_points.push_back(points[0]);
                    }
                    else if(Tool_intersectLines2D(e1,e2,b,c,points) != DISJOINT){
                        intersection_points.push_back(points[0]);
                    }
                    else if(Tool_intersectLines2D(e1,e2,c,a,points) != DISJOINT){
                        intersection_points.push_back(points[0]);
                    }
                    return LTI_INTERSECT_EDG;
                }
            }
            else{
                // 直线与三角形有第二个交点
                if(result == LTI_INTERSECT_NOD){
                    if(cod == 0){  // a
                        if(d1*d >= 0.0 && d3*d >= 0.0){
                            if(d1 == 0.0){  // ab
                                intersection_points.push_back(b);
                                return LTI_INTERSECT_EDG;
                            }
                            else if(d3 == 0.0){  // ac
                                intersection_points.push_back(c);
                                return LTI_INTERSECT_EDG;
                            }
                            else{
                                std::vector<Eigen::Vector2d> points;
                                Tool_intersectLines2D(e1,e2,b,c,points);
                                intersection_points.push_back(points[0]);
                            }
                        }                     
                    }
                    else if(cod == 1){  // b
                        if(d1*d >= 0.0 && d2*d >= 0.0){
                            if(d1 == 0.0){
                                intersection_points.push_back(a);
                                return LTI_INTERSECT_NOD_EDG;
                            }
                            else if(d2 == 0.0){
                                intersection_points.push_back(c);
                                return LTI_INTERSECT_EDG;
                            }
                            else{
                                std::vector<Eigen::Vector2d> points;
                                Tool_intersectLines2D(e1,e2,a,c,points);
                                intersection_points.push_back(points[0]);
                            }
                        }
                    }
                    else if(cod == 2){  // c
                        if(d2*d >= 0.0 && d3*d >= 0.0){
                            if(d2 == 0.0){  // bc
                                intersection_points.push_back(b);
                                return LTI_INTERSECT_EDG;
                            }
                            else if(d3 == 0.0){  // ac
                                intersection_points.push_back(a);
                                return LTI_INTERSECT_EDG;
                            }
                            else{
                                std::vector<Eigen::Vector2d> points;
                                Tool_intersectLines2D(e1,e2,a,b,points);
                                intersection_points.push_back(points[0]);
                            }
                        }
                    }
                }
                else if(result == LTI_INTERSECT_EDG){
                    /* 可能和两位两条边相交 */
                    if(cod == 0){  // ab
                        /* 如果另外1个顶点在边的正侧，存在两个交点 */
                        if(d1*d > 0.0){
                            if(d2*d < 0.0 && d3*d < 0.0){
                                double b_e1e2 = TiGER_GEOM_FUNC::orient2d(real_e1,real_e2,real_b);
                                double c_e1e2 = TiGER_GEOM_FUNC::orient2d(real_e1,real_e2,real_c);
                                if(b_e1e2 == 0.0){
                                    // result = LTI_INTERSECT_NOD;
                                    intersection_points.push_back(b);
                                    return LTI_INTERSECT_NOD;
                                }
                                else if(c_e1e2 == 0.0){
                                    // result = LTI_INTERSECT_NOD;
                                    intersection_points.push_back(c);
                                    return LTI_INTERSECT_NOD;
                                }
                                else if(b_e1e2*c_e1e2 > 0.0){
                                    // ac
                                    cod = 2;
                                }
                            }
                            else if(d2*d < 0.0){
                                cod = 1;
                            }
                            else if(d3*d < 0.0){
                                cod = 2;
                            }

                            if(cod == 1){
                                std::vector<Eigen::Vector2d> points;
                                Tool_intersectLines2D(e1,e2,b,c,points);
                                return LTI_INTERSECT_EDG;
                            }
                            else if(cod == 2){
                                std::vector<Eigen::Vector2d> points;
                                Tool_intersectLines2D(e1,e2,a,c,points);
                                return LTI_INTERSECT_EDG;
                            }
                        }
                        else if(d1*d == 0.0){
                            /* 交点在顶点上，此时，三角形的一条边是原线段的子线段 */
                            if(d2 < 0.0){
                                intersection_points.push_back(b);
                            }
                            else if(d3 < 0.0){
                                intersection_points.push_back(a);
                            }
                            return LTI_INTERSECT_NOD;
                        }
                    }
                    else if(cod == 1){
                        if(d2*d > 0.0){
                            if(d1*d < 0.0 && d3*d < 0.0){
                                double a_e1e2 = TiGER_GEOM_FUNC::orient2d(real_e1,real_e2,real_a);
                                double c_e1e2 = TiGER_GEOM_FUNC::orient2d(real_e1,real_e2,real_c);
                                if(a_e1e2 == 0.0){
                                    intersection_points.push_back(a);
                                    return LTI_INTERSECT_NOD;
                                }
                                else if(c_e1e2 == 0.0){
                                    intersection_points.push_back(c);
                                    return LTI_INTERSECT_NOD;
                                }
                                else if(a_e1e2*c_e1e2 > 0.0){
                                    cod = 0;
                                }
                            }
                            else if(d1*d < 0.0){
                                cod = 0;
                            }
                            else if(d3*d < 0.0){
                                cod = 2;
                            }

                            if(cod == 0){  // ab
                                std::vector<Eigen::Vector2d> points;
                                Tool_intersectLines2D(e1,e2,a,b,points);
                                return LTI_INTERSECT_EDG;
                            }
                            else if(cod == 2){  // ac
                                std::vector<Eigen::Vector2d> points;
                                Tool_intersectLines2D(e1,e2,a,c,points);
                                return LTI_INTERSECT_EDG;
                            }
                        }
                        else if(d2 == 0.0){
                            if(d1 < 0.0){
                                intersection_points.push_back(b);
                            }
                            else if(d3 < 0.0){
                                intersection_points.push_back(c);
                            }
                            return LTI_INTERSECT_NOD;
                        }
                    }
                    else if(cod == 2){
                        if(d3*d > 0.0){
                            if(d1*d < 0.0 && d2*d < 0.0){
                                double b_e1e2 = TiGER_GEOM_FUNC::orient2d(real_e1,real_e2,real_b);
                                double a_e1e2 = TiGER_GEOM_FUNC::orient2d(real_e1,real_e2,real_a);
                                if(b_e1e2 == 0.0){
                                    intersection_points.push_back(b);
                                    return LTI_INTERSECT_NOD;
                                }
                                else if(a_e1e2 == 0.0){
                                    intersection_points.push_back(a);
                                    return LTI_INTERSECT_NOD;
                                }
                                else if(b_e1e2*a_e1e2 > 0.0){
                                    cod = 1;
                                }
                            }
                            else if(d1*d < 0.0){
                                cod = 0;
                            }
                            else if(d2*d < 0.0){
                                cod = 1;
                            }

                            if(cod == 0){  //ab
                                std::vector<Eigen::Vector2d> points;
                                Tool_intersectLines2D(e1,e2,a,b,points);
                                return LTI_INTERSECT_EDG;
                            }
                            else if(cod == 1){
                                std::vector<Eigen::Vector2d> points;
                                Tool_intersectLines2D(e1,e2,b,c,points);
                                return LTI_INTERSECT_EDG;
                            }
                        }
                        else if(d3 == 0.0){
                            if(d1 < 0.0){
                                intersection_points.push_back(a);
                            }
                            else if(d2 < 0.0){
                                intersection_points.push_back(c);
                            }
                            return LTI_INTERSECT_NOD;
                        }
                    }
                }
            }
        }
        else{  // 都在外部
            auto onSegment = [](const Eigen::Vector2d& pi, const Eigen::Vector2d& pj,
                      const Eigen::Vector2d& pk) {
                return ((std::min(pi.x(), pj.x()) <= pk.x() &&
                        pk.x() <= std::max(pi.x(), pj.x()) &&
                        std::min(pi.y(), pj.y()) <= pk.y() &&
                        pk.y() <= std::max(pi.y(), pj.y())));
            };
            // 第一种情况  面的某条边是线的子线段
            if(d1 == 0 && d4 == 0){ // ab
                if(onSegment(e1,e2,a) && onSegment(e1,e2,b)){
                    intersection_points.push_back(a);
                    intersection_points.push_back(b);
                    return LTI_INTERSECT_EDG;
                }
                else{
                    return DISJOINT;
                }
            }
            if(d2 == 0 && d5 == 0){  // bc
                if(onSegment(e1,e2,b) && onSegment(e1,e2,c)){
                    intersection_points.push_back(b);
                    intersection_points.push_back(c);
                    return LTI_INTERSECT_EDG;
                }
                else{
                    return DISJOINT;
                }
            }
            if(d3 == 0 && d6 == 0){  // ca
                if(onSegment(e1,e2,c) && onSegment(e1,e2,a)){
                    intersection_points.push_back(c);
                    intersection_points.push_back(a);
                    return LTI_INTERSECT_EDG;
                }
                else{
                    return DISJOINT;
                }
            }

            // 第二种情况 三角形一个端点在线上 
            double a_on = TiGER_GEOM_FUNC::orient2d(real_e1,real_e2,real_a);
            double b_on = TiGER_GEOM_FUNC::orient2d(real_e1,real_e2,real_b);
            double c_on = TiGER_GEOM_FUNC::orient2d(real_e1,real_e2,real_c);
            if(a_on == 0.0){  // a 共线
                if(onSegment(e1,e2,a)){  // a 在线上
                    intersection_points.push_back(a);
                    std::vector<Eigen::Vector2d> points;
                    if(Tool_intersectLines2D(e1,e2,b,c,points) == INTERSECT){
                        intersection_points.push_back(points[0]);
                        return LTI_INTERSECT_NOD_EDG;  // 与端点和对边相交
                    }
                    return LTI_INTERSECT_NOD;
                }
                else{  // 在延长线上
                    return DISJOINT;
                }
            }
            if(b_on == 0.0){  // b 共线
                if(onSegment(e1,e2,b)){
                    intersection_points.push_back(b);
                    std::vector<Eigen::Vector2d> points;
                    if(Tool_intersectLines2D(e1,e2,a,c,points) == INTERSECT){
                        intersection_points.push_back(points[0]);
                        return LTI_INTERSECT_NOD_EDG;
                    }
                    return LTI_INTERSECT_NOD;
                }
                else{
                    return DISJOINT;
                }
            }
            if(c_on == 0.0){  // c 共线
                if(onSegment(e1,e2,c)){
                    intersection_points.push_back(c);
                    std::vector<Eigen::Vector2d> points;
                    if(Tool_intersectLines2D(e1,e2,a,b,points) == INTERSECT){
                        intersection_points.push_back(points[0]);
                        return LTI_INTERSECT_NOD_EDG;
                    }
                    return LTI_INTERSECT_NOD;
                }
                else{
                    return DISJOINT;
                }
            }

            // 第三种情况 面的2条边和直线有交点
            std::vector<Eigen::Vector2d> points;  // 与所有边的交点
            Tool_intersectLines2D(e1,e2,a,b,points);
            Tool_intersectLines2D(e1,e2,b,c,points);
            Tool_intersectLines2D(e1,e2,c,a,points);
            if(points.size() > 0){
                for(auto& point : points)
                    intersection_points.push_back(point);
                return LTI_INTERSECT_EDG;
            }
            else{
                return DISJOINT;
            }
        }
    }

enum intersectionresult Tool_intersectPointWithTriangle(
    const Eigen::RowVector3d& p,
    const Eigen::RowVector3d& a,
    const Eigen::RowVector3d& b,
    const Eigen::RowVector3d& c){
        if(p == a || p == b || p == c)
            return NOD_INTERSECT_NOD;
        REAL* real_p = const_cast<REAL*>(p.data());
        REAL* real_a = const_cast<REAL*>(a.data());
        REAL* real_b = const_cast<REAL*>(b.data());
        REAL* real_c = const_cast<REAL*>(c.data());
        double d = TiGER_GEOM_FUNC::orient3d(real_p,real_a,real_b,real_c);
        if(d != 0)  // 不共面
            return DISJOINT;
        Eigen::RowVector3d abc_normal = (b - a).cross(c - a).normalized();
        std::vector<Eigen::RowVector3d> axis_normal = {Eigen::RowVector3d(1,0,0),
                                                        Eigen::RowVector3d(0,1,0),
                                                        Eigen::RowVector3d(0,0,1)};
        int axis = 0;
        for(int i=0;i<3;i++){
            if(fabs(abc_normal.dot(axis_normal[i])) > 0.5){
                axis = i;
                break;
            }
        }
        int x = (axis + 1) % 3;
        int y = (axis + 2) % 3;
        Eigen::Vector2d p_,a_,b_,c_;
        p_(0) = p(x);
        p_(1) = p(y);
        a_(0) = a(x);
        a_(1) = a(y);
        b_(0) = b(x);
        b_(1) = b(y);
        c_(0) = c(x);
        c_(1) = c(y);
        double s1 = TiGER_GEOM_FUNC::orient2d(p_.data(), a_.data(), b_.data());
        double s2 = TiGER_GEOM_FUNC::orient2d(p_.data(), b_.data(), c_.data());
        double s3 = TiGER_GEOM_FUNC::orient2d(p_.data(), c_.data(), a_.data());
        if(s1*s2*s3 > 0) //  在内部
            return NOD_INTERSECT_FAC;
        else if(s1 == 0 || s2 == 0 || s3 == 0)
            return NOD_INTERSECT_EDG;
        return DISJOINT;
    }

} // namespace repair



} // namespace TiGER


