#pragma once

#include "mesh_data_struct.h"
#include "definitions.h"
#include <functional>
#include "Eigen/Dense"
#include "alias.h"

namespace TiGER {
namespace repair {
/**
 * @brief
 *
 * @param[in] points_in 输入的点
 * @param[in] eps 去重的容差
 * @param[out] bcj_pre 并查集
 *
 * points_in index =  0   1   2  3  4  5  6   7  /// 2，3，5，6 重复，4，7重复。
 * bcj_pre         = -1  -1  -1  2 -1  2  2   4  --> 并查集结果
 *
 * @return int 重复点的个数
 */
int Tool_detectPointDuplicates(
    const std::vector<std::array<double, 3>>& points_in,
    const double& eps, 
    std::vector<int>& bcj_pre);

/**
 * @brief 暴力法检测重复点
 *
 * @param[in] points_in 输入的点
 * @param[in] eps 去重的容差
 * @param[out] bcj_pre 并查集
 *
 * points_in index =  0   1   2  3  4  5  6   7  /// 2，3，5，6 重复，4，7重复。
 * bcj_pre         = -1  -1  -1  2 -1  2  2   4  --> 并查集结果
 *
 * @return int 重复点的个数
 */
int Tool_nativelyDetectPointDuplicates(
    const std::vector<std::array<double, 3>>& points_in,
    const double& eps, 
    std::vector<int>& bcj_pre);

/**
 * @brief 点-线段的投影点&距离
 *
 * @param[in] p 点
 * @param[in] a 线段端点
 * @param[in] b 线段端点
 * @param[out] project_point  投影点
 * @param[out] point_to_segment_distance 距离
 * 
 */
void Tool_projectPointToSegment(
    const Eigen::RowVector3d& p,
    const Eigen::RowVector3d& a,
    const Eigen::RowVector3d& b,
    Eigen::RowVector3d& project_point,
    double& point_to_segment_distance);

/**
 * @brief 点-三角形的投影点&距离
 *
 * @param[in] p 点
 * @param[in] a 三角形端点
 * @param[in] b 三角形端点
 * @param[in] c 三角形端点
 * @param[out] project_point  投影点
 * @param[out] point_to_tri_distance 距离
 * 
 */
void Tool_projectPointToTriangle(
    const Eigen::RowVector3d& p,
    const Eigen::RowVector3d& a,
    const Eigen::RowVector3d& b,
    const Eigen::RowVector3d& c,
    Eigen::RowVector3d& project_point,
    double& point_to_tri_distance);

/**
 * @brief 线-单个三角形的投影
 *
 * @param[in] e1 线段端点
 * @param[in] e2 线段端点
 * @param[in] a 三角形端点
 * @param[in] b 三角形端点
 * @param[in] c 三角形端点
 * @param[out] project_edge  投影边  可能没有
 * 
 */
void project_SegmentOntoTriangle(
    const Eigen::RowVector3d& e1,
    const Eigen::RowVector3d& e2,
    const Eigen::RowVector3d& a,
    const Eigen::RowVector3d& b,
    const Eigen::RowVector3d& c,
    std::vector<Eigen::RowVector3d>& project_edge);

/**
 * @brief 线线投影
 *
 * @param[in] a 线段端点
 * @param[in] b 线段端点
 * @param[in] c 线段端点
 * @param[in] d 线段端点
 * @param[out] closestPoint1  ab上的投影点
 * @param[out] closestPoint2  cd上的投影点
 * 
 */
double project_SegmentOntoSegment(
    const Eigen::RowVector3d& a,
    const Eigen::RowVector3d& b,
    const Eigen::RowVector3d& c,
    const Eigen::RowVector3d& d,
    Eigen::RowVector3d& closestPoint1,
    Eigen::RowVector3d& closestPoint2); // 投影点

bool isPointInTriangle(
    const Eigen::RowVector3d& p, 
    const Eigen::RowVector3d& a,
    const Eigen::RowVector3d& b, 
    const Eigen::RowVector3d& c);

/**
 * @brief 线-曲面投影
 *
 * @param[in] e1 线段端点
 * @param[in] e2 线段端点
 * @param[in] 曲面
 * @param[in] eps
 * @param[out] pro_points  投影点
 * 
 */
void project_SegmentOntoSurface(
    const Eigen::RowVector3d& e1,
    const Eigen::RowVector3d& e2,
    const SurfaceMesh& surf_in,
    const double& eps,
    std::vector<Eigen::RowVector3d>& pro_points);

/**
 * @brief 边插点
 *
 * @param[in] surf_in 输入网格
 * @param[in] edge_in 对应边
 * @param[in] points 插入点序列
 * @param[out] surf_out  输出网格
 * 
 */
void edge_repair(
    const SurfaceMesh& surf_in,
    const std::array<int,2>& edge_in,
    const std::vector<Eigen::RowVector3d>& points,
    SurfaceMesh& surf_out);

/**
 * @brief 线分割三角形
 *
 * @param[in] a 三角形端点
 * @param[in] b 三角形端点
 * @param[in] c 三角形端点
 * @param[in] project_poins 端点
 * @param[in] vertextid 默认[-1，-1]  线段端点在 三角形端点
 * @param[in] edgeid 默认[-1，-1] 线段端点在边上
 * @param[out] surf_out  分割结果
 * 
 */
void partition(
    const Eigen::RowVector3d& a,
    const Eigen::RowVector3d& b,
    const Eigen::RowVector3d& c,
    const std::vector<Eigen::RowVector3d>& project_poins,
    const std::array<int,2>& vertextid,
    const std::array<int,2>& edgeid,
    TiGER::SurfaceMesh& surf_out);

enum intersectionresult{
    // line-line
    DISJOINT = 0,  // 分离
    SHAREVERTEX,  // 共享节点
    SHAREEDGE,   // 共享边
    INTERSECT,   // 相交
    //line-tri
    // LINE_INTERSECT_NUL,
    LTI_INTERSECT_NOD,  // 相交于三角形端点
	LTI_INTERSECT_EDG,  // 相交于三角形边
	LTI_INTERSECT_FAC,  // 相交于三角形内
	LTI_INTERSECT_INS,  // 在三角形内部
    LTI_INTERSECT_EDG_FACE,  // 一个点在内部  一个在
    LTI_INTERSECT_NOD_EDG,
    LTI_INTERSECT_EDG_EDG,
    // point-tri
    NOD_INTERSECT_NOD,
    NOD_INTERSECT_EDG,
    NOD_INTERSECT_FAC
};

// enum line_tri_intersection_type{
//     LINE_INTERSECT_NUL,
//     LTI_INTERSECT_NOD,  // 相交于三角形端点
// 	LTI_INTERSECT_EDG,  // 相交于三角形边
// 	LTI_INTERSECT_FAC,  // 相交于三角形内
// 	LTI_INTERSECT_INS,  // 在三角形内部
//     LTI_INTERSECT_EDG_FACE  // 
// };

/**
 * @brief 检测3d线线相交
 *
 * @param[in] e1 线段1端点
 * @param[in] e2 线段1端点
 * @param[in] e3 线段2端点
 * @param[in] e4 线段2端点
 * @param[out] intersection_points 交点
 * 
 * @return 相交类型  DISJOINT  没有交点 intersection_points.size() = 0
 *                  SHAREVERTEX  一个交点  intersection_points.size() = 1
 *                  SHAREEDGE  两个交点 intersection_points.size() = 2
 *                  INTERSECT  一个交点 intersection_points.size() = 1
 */
enum intersectionresult Tool_intersectLines(
    const Eigen::RowVector3d& e1,
    const Eigen::RowVector3d& e2,
    const Eigen::RowVector3d& e3,
    const Eigen::RowVector3d& e4,
    std::vector<Eigen::RowVector3d>& intersection_points);

/**
 * @brief 检测2d线线相交
 *
 * @param[in] e1 线段1端点
 * @param[in] e2 线段1端点
 * @param[in] e3 线段2端点
 * @param[in] e4 线段2端点
 * @param[out] intersection_points 交点
 * 
 * @return 相交类型  DISJOINT  没有交点 intersection_points.size() = 0
 *                  SHAREVERTEX  一个交点  intersection_points.size() = 1
 *                  SHAREEDGE  两个交点 intersection_points.size() = 2
 *                  INTERSECT  一个交点 intersection_points.size() = 1
 */
enum intersectionresult Tool_intersectLines2D(
    const Eigen::Vector2d& e1,
    const Eigen::Vector2d& e2,
    const Eigen::Vector2d& e3,
    const Eigen::Vector2d& e4,
    std::vector<Eigen::Vector2d>& intersection_points);

enum intersectionresult Tool_intersectLineWithTriangle(
    const Eigen::RowVector3d& e1,
    const Eigen::RowVector3d& e2,
    const Eigen::RowVector3d& a,
    const Eigen::RowVector3d& b,
    const Eigen::RowVector3d& c,
    std::vector<Eigen::RowVector3d>& intersection_points);

enum intersectionresult Tool_intersectLineWithTriangle2D(
    const Eigen::Vector2d& e1, 
    const Eigen::Vector2d& e2,
    const Eigen::Vector2d& a, 
    const Eigen::Vector2d& b,
    const Eigen::Vector2d& c,
    std::vector<Eigen::Vector2d>& intersection_points);

enum intersectionresult Tool_intersectPointWithTriangle(
    const Eigen::RowVector3d& p,
    const Eigen::RowVector3d& a,
    const Eigen::RowVector3d& b,
    const Eigen::RowVector3d& c);

} // namespace repair
} // namespace TiGER

#include "mesh_repair.cpp"