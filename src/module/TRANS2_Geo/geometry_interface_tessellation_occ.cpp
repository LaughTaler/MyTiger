#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <Eigen/Dense>

#include "meshIO.h"
#include "list_to_matrix.hpp"
#include "../../include/alias.h"
#include "../../include/geometry_data_struct.h"
#include "./geometry_occ.h"
#include "./geometry_ferguson.h"
#include "../../include/mesh_data_struct.h"
#include "../../include/definitions.h"
#include "../../include/discretize_surface_delaunay.h"
#include "../../include/discretize_surface_delaunay_parameters.h"
#include "../../include/geometry_interface_occ.h"
// #include "dt2D_API.h"
// #include "dt2D.h"
using namespace std;

namespace TiGER
{
namespace discretecurcestl
{
void discreteCurve_uniform_test(const std::shared_ptr<TiGER::GeometryCurve>& fcurve, const int dimension,
                                const double delta_s, TiGER::LineMesh& cur_curve_mesh)
{
    double begin_u = 0.;
    double end_u = 1.;
    if( (dimension == 0 || dimension == 1 || dimension == 2) && delta_s == 0.0 )
    {
        std::array<double, 3> xyz;
        double u;
        // fcurve->arcToCoordinate(0.0, xyz, u);
        fcurve->arcToCoordinateFunction(0.0, xyz, u);
        cur_curve_mesh.coord.push_back(xyz);
        fcurve->arcToCoordinateFunction(fcurve->getTotalLenthFunction(), xyz, u);
        cur_curve_mesh.coord.push_back(xyz);
        return;
    }
    double total_lenth = fcurve->getTotalLenthFunction();
    double ave_s;
    int dim;
    if( delta_s == 0.0 )
    {
        ave_s = total_lenth / (dimension - 1);
        dim = dimension;
    }
    else
    {
        ave_s = dimension == 0 ? delta_s : min(delta_s, total_lenth / (dimension - 1));
        dim = (int)(total_lenth / ave_s);
        if( total_lenth / ave_s - (double)dim != 0.0 )
            dim++;
        ave_s = total_lenth / dim;
        dim++;
    }
    double cur_location = 0.0;
    std::array<double, 3> xyz;
    double u;
    dim--;
    // fcurve->arcToCoordinate(cur_location, xyz, u);
    fcurve->arcToCoordinateFunction(cur_location, xyz, u);
    cur_curve_mesh.coord.push_back(xyz);
    for( int i = 0; i < dim; i++ )
    {
        cur_location += ave_s;
        cur_location = fmin(cur_location, total_lenth);
        // fcurve->arcToCoordinate(cur_location, xyz, u);
        fcurve->arcToCoordinateFunction(cur_location, xyz, u);
        cur_curve_mesh.coord.push_back(xyz);
    }
    return;
}
} // namespace discretecurcestl
namespace projecthelp
{
vec2d arraytovec2d(std::array<double, 2>& aaa)
{
    vec2d a = {aaa[0], aaa[1]};
    return a;
}
vec3d arraytovec3d(std::array<double, 3>& aaa)
{
    vec3d a = {aaa[0], aaa[1], aaa[2]};
    return a;
}
std::array<double, 2> vec2dtoarray(vec2d aaa)
{
    std::array<double, 2> a = {aaa[0], aaa[1]};
    return a;
}
std::array<double, 3> vec3dtoarray(vec3d aaa)
{
    std::array<double, 3> a = {aaa[0], aaa[1], aaa[2]};
    return a;
}
int sf_project(const std::shared_ptr<TiGER::GeometrySurface>& surface, const std::array<double, 3>& xyz,
               const double eps, std::array<double, 2>& res_uv)
{
    double dis = (std::numeric_limits<double>::max)();
    int itNum = 100;
    double eps2 = 1e-10;
    vec2d iterUV = arraytovec2d(res_uv);

    vec3d P, Q0, d1u, d1v, d2u2, d2v2, d2uv;
    std::array<double, 3> Q0_, d1u_, d1v_, d2u2_, d2v2_, d2uv_;
    P.x() = xyz[0];
    P.y() = xyz[1];
    P.z() = xyz[2];
    // P = xyz;

    int i;
    for( i = 0; i < itNum; i++ )
    {
        // surface->d0(common::vec2dtoarray(iterUV), Q0_);
        Q0_ = surface->d0Function(vec2dtoarray(iterUV));
        Q0 = arraytovec3d(Q0_);
        //	//this->param_to_uv_d1(iterUV, d1u, d1v);
        // surface->d1(common::vec2dtoarray(iterUV), d1u_, d1v_);
        surface->d1Function(vec2dtoarray(iterUV), d1u_, d1v_);
        d1u = arraytovec3d(d1u_);
        d1v = arraytovec3d(d1v_);
        //	//this->param_to_uv_d2(iterUV, d2u2, d2v2, d2uv);
        // surface->d2(common::vec2dtoarray(iterUV), d2u2_, d2v2_, d2uv_);
        surface->d2Function(vec2dtoarray(iterUV), d2u2_, d2v2_, d2uv_);
        d2u2 = arraytovec3d(d2u2_);
        d2v2 = arraytovec3d(d2v2_);
        d2uv = arraytovec3d(d2uv_);
        // 待投影点与曲面上猜测点的差向量
        vec3d p_p0 = P - Q0;

        // 法向量
        vec3d normal = d1u.cross(d1v); // 叉乘
        normal.normalize();

        // 参数计算
        double E = d1u.dot(d1u);
        double F = d1u.dot(d1v);
        double G = d1v.dot(d1v);

        // double C1 = d1u.dot(normal);
        // double C2 = d1v.dot(normal);

        // double S1 = p_p0.dot(d1u);
        // double S2 = p_p0.dot(d1v);
        // double S3 = p_p0.dot(normal);

        double L = d2u2.dot(normal);
        double M = d2uv.dot(normal);
        double N = d2v2.dot(normal);

        // double denominator = C1 * C1 * G - 2 * C1 * C2 * F + C2 * C2 * E - E * G
        // + F * F;

        // double a1 = -(C1 * C2 * S2 - C1 * G * S3 - C2 * C2 * S1 + C2 * F * S3 - F
        // * S2 + G * S1) / denominator; double a2 = (C1 * C1 * S2 - C1 * C2 * S1 -
        // C1 * F * S3 + C2 * E * S3 - E * S2 + F * S1) / denominator;

        ////计算法曲率
        // double kn = (L * a1 * a1 + 2 * M * a1 * a2 + N * a2 * a2) / (E * a1 * a1
        // + 2 * F * a1 * a2 + G * a2 * a2);

        double dudu = d1u.dot(d1u);
        double dudv = d1u.dot(d1v);
        double dvdv = d1v.dot(d1v);

        double kn = (L * dudu + 2 * M * dudv + N * dvdv) / (E * dudu + 2 * F * dudv + G * dvdv);

        // 退化为切平面法
        if( std::isnan(kn) || kn < 1e-6 )
            kn = 1e-6;

        // 曲率半径
        double R = fabs(1.0 / kn);

        // 曲率中心
        vec3d m = Q0 + R * normal;

        vec3d m_P = m - P;

        double mod_m_P = m_P.squaredNorm();
        double w = 1 - R / sqrt(mod_m_P);

        vec3d q = P + w * m_P;
        vec3d q_Q0 = q - Q0;

        double C4 = q_Q0.dot(d1u);
        double C5 = q_Q0.dot(d1v);

        double deltaU = (C4 * G - C5 * F) / (E * G - F * F);
        double deltaV = -(C4 * F - C5 * E) / (E * G - F * F);

        if( fabs(deltaU) < eps2 && fabs(deltaV) < eps2 )
        {
            break;
        }

        iterUV = iterUV + vec2d(deltaU, deltaV);

        std::array<double, 4> uv_scale = surface->getUVScaleFunction();

        iterUV.x() = std::max(iterUV.x(), uv_scale[0]);
        iterUV.x() = std::min(iterUV.x(), uv_scale[1]);

        iterUV.y() = std::max(iterUV.y(), uv_scale[2]);
        iterUV.y() = std::min(iterUV.y(), uv_scale[3]);

        // surface->d0(common::vec2dtoarray(iterUV), Q0_);
        Q0_ = surface->d0Function(vec2dtoarray(iterUV));
        Q0 = arraytovec3d(Q0_);
        //	//this->param_to_uv_d1(iterUV, d1u, d1v);
        // surface->d1(common::vec2dtoarray(iterUV), d1u_, d1v_);
        // surface->d1Function(common::vec2dtoarray(iterUV), d1u_, d1v_);
        // d1u = common::arraytovec3d(d1u_); d1v = common::arraytovec3d(d1v_);
        //	//this->param_to_uv_d2(iterUV, d2u2, d2v2, d2uv);
        // surface->d2(common::vec2dtoarray(iterUV), d2u2_, d2v2_, d2uv_);
        // surface->d2Function(common::vec2dtoarray(iterUV), d2u2_, d2v2_, d2uv_);
        // d2u2 = common::arraytovec3d(d2u2_); d2v2 = common::arraytovec3d(d2v2_);
        // d2uv = common::arraytovec3d(d2uv_);

        dis = (P - Q0).squaredNorm();

        if( dis < eps ) // 应该不可能达到
        {
            // std::cout << (*dis) << std::endl;
            break;
        }
        // p_p0 = P - Q0;
        // 曲面上猜测点的切向量
        // double otho1 = p_p0.dot(d1u);
        // double otho2 = p_p0.dot(d1v);

        //// //这条判据在边界处可能失效，因此只在非边界处使用
        // if (iterUV.x() > uv_scale[0] + 1e-2 && iterUV.x() < uv_scale[1] - 1e-2 &&
        //	iterUV.y() > uv_scale[2] + 1e-2 && iterUV.y() < uv_scale[3] - 1e-2)
        //	if (fabs(otho1) < eps2 && fabs(otho2) < eps2)
        //	{
        //		break;
        //	}
    }

    if( (dis) > eps || std::isnan(dis) ) // 投影错误，使用暴力算法
    {
        // sf_project_robust(surface, xyz, res_uv, dis);
        //  baoli_pro++;
        Eigen::Vector3d xyz_coord(xyz[0], xyz[1], xyz[2]);
        std::array<double, 4> uv_scale = surface->getUVScaleFunction();
        Eigen::Vector2d uv_min = Eigen::Vector2d(uv_scale[0], uv_scale[2]);
        Eigen::Vector2d uv_max = Eigen::Vector2d(uv_scale[1], uv_scale[3]);
        Eigen::Vector2d uv_min_init = uv_min;
        Eigen::Vector2d uv_max_init = uv_max;
        int num_u = 20;
        int num_v = 20;
        int iterater_num = 10;
        double best_val;
        double pre_best_val = std::numeric_limits<double>::max();
        while( iterater_num-- )
        {
            std::pair<int, int> best_loc;
            best_val = std::numeric_limits<double>::max();
            double uv_max_min[2] = {uv_max.x() - uv_min.x(), uv_max.y() - uv_min.y()};
            double uv_min_pre[2] = {uv_min.x(), uv_min.y()};
            for( int i = 0; i < num_u; ++i )
            {
                for( int j = 0; j < num_v; ++j )
                {
                    Eigen::Vector2d current_uv_mid =
                        Eigen::Vector2d((1.0 * (i + 0.5) / num_u) * uv_max_min[0] + uv_min_pre[0],
                                        (1.0 * (j + 0.5) / num_v) * uv_max_min[1] + uv_min_pre[1]);
                    Eigen::Vector3d current_xyz_mid;
                    std::array<double, 3> xyz_array;
                    // surface->d0(common::vec2dtoarray(current_uv_mid), xyz_array);
                    xyz_array = surface->d0Function(vec2dtoarray(current_uv_mid));

                    // current_xyz_mid = common::arraytovec3d(xyz_array);
                    // double distance = (current_xyz_mid - xyz_coord).norm();
                    double distance = sqrt((xyz_array[0] - xyz[0]) * (xyz_array[0] - xyz[0]) +
                                           (xyz_array[1] - xyz[1]) * (xyz_array[1] - xyz[1]) +
                                           (xyz_array[2] - xyz[2]) * (xyz_array[2] - xyz[2]));

                    if( distance < best_val )
                    {
                        best_val = distance;
                        best_loc = std::make_pair(i, j);
                    }
                }
            }
            uv_min = Eigen::Vector2d((1.0 * (best_loc.first - 1) / num_u) * uv_max_min[0] + uv_min_pre[0],
                                     (1.0 * (best_loc.second - 1) / num_v) * uv_max_min[1] + uv_min_pre[1]);
            uv_max = Eigen::Vector2d((1.0 * (best_loc.first + 2) / num_u) * uv_max_min[0] + uv_min_pre[0],
                                     (1.0 * (best_loc.second + 2) / num_v) * uv_max_min[1] + uv_min_pre[1]);
            uv_min = uv_min.cwiseMax(uv_min_init);
            uv_max = uv_max.cwiseMin(uv_max_init);
            pre_best_val = best_val;
        }
        if( best_val < dis )
            iterUV = (uv_min + uv_max) * 0.5;
        res_uv = vec2dtoarray(iterUV);
        // spdlog::info("{}%", (dis - best_val) / dis*100);
        return 1;
    }

    res_uv = vec2dtoarray(iterUV);

    return 1;
}

} // namespace projecthelp

namespace help
{

// 定义采样点的结构
struct SamplePoint
{
    double u, v;
    static bool isEqual(double a, double b, double tolerance = 1e-5)
    {
        return std::fabs(a - b) <= tolerance;
    }

    bool operator<(const SamplePoint& other) const
    {
        if( isEqual(u, other.u) )
        {
            return v < other.v && !isEqual(v, other.v);
        }
        return u < other.u;
    }

    bool operator==(const SamplePoint& other) const
    {
        return isEqual(u, other.u) && isEqual(v, other.v);
    }

    double distance(const SamplePoint& other) const
    {
        return std::sqrt(std::pow(u - other.u, 2) + std::pow(v - other.v, 2));
    }
};

// 定义曲线上的采样点结构
struct SamplePointCurve
{
    double u; // 曲线上的参数值

    bool operator<(const SamplePointCurve& other) const
    {
        return u < other.u;
    }

    bool operator==(const SamplePointCurve& other) const
    {
        return u == other.u;
    }
};
// 定义子区域的结构，由四个采样点组成
struct SubRegion
{
    SamplePoint p1, p2, p3, p4;
};
// 定义曲线上的子区域结构
struct SubRegionCurve
{
    SamplePointCurve p1, p2; // 子区域的两个端点
};

int find_fa_id(int x, vector<int>& fa_id)
{
    if( fa_id[x] == x )
        return x;
    else
        return find_fa_id(fa_id[x], fa_id);
}
void merge_fa_id(int x, int y, vector<int>& fa_id)
{
    fa_id[find_fa_id(x, fa_id)] = find_fa_id(y, fa_id);
}
vec3d arraytovec3d(std::array<double, 3>& p)
{
    vec3d res(p[0], p[1], p[2]);
    return res;
}
void dAssemblecurve(std::vector<std::vector<std::array<double, 3>>>& curve_mesh_3d_lst)
{
    // std::cout << "use curveAssemble\n";
    std::vector<int> fa_id; // 并查集的父亲id
    int n = curve_mesh_3d_lst.size();
    fa_id.resize(2 * n + 1);
    // 曲线i的首节点id设置为2*i-1;曲线i的末节点id设置为2*i
    for( int i = 1; i <= n; i++ )
    {
        fa_id[2 * i - 1] = 2 * i - 1;
        fa_id[2 * i] = 2 * i;
    }
    // 曲线比较少用暴力合并
    if( curve_mesh_3d_lst.size() < 1000 )
    {
        for( int i = 1; i <= curve_mesh_3d_lst.size(); i++ )
        {
            double start_dis = (std::numeric_limits<double>::max)();
            double end_dis = (std::numeric_limits<double>::max)();
            vec3d start_pt = arraytovec3d(curve_mesh_3d_lst[i - 1].front());
            vec3d end_pt = arraytovec3d(curve_mesh_3d_lst[i - 1].back());
            int start_merge_id = 2 * i - 1;
            int end_merge_id = 2 * i;
            for( int j = 1; j <= curve_mesh_3d_lst.size(); j++ )
            {
                if( j == i )
                    continue;
                vec3d cur_start_pt = arraytovec3d(curve_mesh_3d_lst[j - 1].front());
                vec3d cur_end_pt = arraytovec3d(curve_mesh_3d_lst[j - 1].back());
                if( (start_pt - cur_start_pt).norm() < start_dis )
                {
                    start_merge_id = 2 * j - 1;
                    start_dis = (start_pt - cur_start_pt).norm();
                }
                if( (start_pt - cur_end_pt).norm() < start_dis )
                {
                    start_merge_id = 2 * j;
                    start_dis = (start_pt - cur_end_pt).norm();
                }
                if( (end_pt - cur_start_pt).norm() < end_dis )
                {
                    end_merge_id = 2 * j - 1;
                    end_dis = (end_pt - cur_start_pt).norm();
                }
                if( (end_pt - cur_end_pt).norm() < end_dis )
                {
                    end_merge_id = 2 * j;
                    end_dis = (end_pt - cur_end_pt).norm();
                }
            }
            merge_fa_id(2 * i - 1, start_merge_id, fa_id);
            merge_fa_id(2 * i, end_merge_id, fa_id);
        }
    }
    // 曲线比较多用kdtree合并
    /*else {
      KdTreeM kdtree(3);
      for (int i = 1; i <= n; i++) {
        kdtree.Insert3DNode(common::arraytovec3d(curve_mesh_3d_lst[i -
    1].front()), 2 * i - 1);
        kdtree.Insert3DNode(common::arraytovec3d(curve_mesh_3d_lst[i -
    1].back()), 2 * i);
      }
      for (int i = 1; i <= curve_mesh_3d_lst.size(); i++) {
        vec3d start_pt = common::arraytovec3d(curve_mesh_3d_lst[i -
    1].front()); vec3d end_pt = common::arraytovec3d(curve_mesh_3d_lst[i -
    1].back()); int start_merge_id = kdtree.nearest3DNode(start_pt); int
    end_merge_id = kdtree.nearest3DNode(end_pt); common::merge_fa_id(2 * i -
    1, start_merge_id, fa_id); common::merge_fa_id(2 * i, end_merge_id,
    fa_id);
      }
    }*/
    for( int i = 1; i <= curve_mesh_3d_lst.size(); i++ )
    {
        int start_fa_id = find_fa_id(2 * i - 1, fa_id);
        int start_fa_cv_id = (start_fa_id + 1) / 2 - 1;
        int end_fa_id = find_fa_id(2 * i, fa_id);
        int end_fa_cv_id = (end_fa_id + 1) / 2 - 1;
        if( start_fa_id % 2 == 0 )
        {
            curve_mesh_3d_lst[i - 1].front() = curve_mesh_3d_lst[start_fa_cv_id].back();
        }
        else
        {
            curve_mesh_3d_lst[i - 1].front() = curve_mesh_3d_lst[start_fa_cv_id].front();
        }
        if( end_fa_id % 2 == 0 )
        {
            curve_mesh_3d_lst[i - 1].back() = curve_mesh_3d_lst[end_fa_cv_id].back();
        }
        else
        {
            curve_mesh_3d_lst[i - 1].back() = curve_mesh_3d_lst[end_fa_cv_id].front();
        }
    }
}

double ProjectionDistancePointToSegment(const SamplePoint& p, const SamplePoint& p1, const SamplePoint& p2)
{
    Eigen::Vector2d point(p.u, p.v);
    Eigen::Vector2d segment_start(p1.u, p1.v);
    Eigen::Vector2d segment_end(p2.u, p2.v);
    Eigen::Vector2d segment = segment_end - segment_start;

    Eigen::Vector2d point_to_start = point - segment_start;
    double projection = point_to_start.dot(segment) / segment.squaredNorm();

    if( projection < 0 )
    {
        // Point projects outside the segment, closer to p1
        return (point - segment_start).norm();
    }
    else if( projection > 1 )
    {
        // Point projects outside the segment, closer to p2
        return (point - segment_end).norm();
    }
    else
    {
        // Point projects onto the segment
        Eigen::Vector2d projection_point = segment_start + projection * segment;
        return (point - projection_point).norm();
    }
}

void computeBoundingBox(const LineMesh& mesh, double& umin, double& umax, double& vmin, double& vmax)
{
    umin = mesh.coord[0][0];
    umax = mesh.coord[0][0];
    vmin = mesh.coord[0][1];
    vmax = mesh.coord[0][1];

    for( const auto& segment : mesh.segments )
    {
        // 获取每个线段的两个端点
        const auto& point1 = mesh.coord[segment[0]];
        const auto& point2 = mesh.coord[segment[1]];

        // 更新umin, umax, vmin, vmax
        umin = std::min(umin, point1[0]);
        umax = std::max(umax, point1[0]);
        umin = std::min(umin, point2[0]);
        umax = std::max(umax, point2[0]);

        vmin = std::min(vmin, point1[1]);
        vmax = std::max(vmax, point1[1]);
        vmin = std::min(vmin, point2[1]);
        vmax = std::max(vmax, point2[1]);
    }
}

double Distance(const std::array<double, 3>& a, const std::array<double, 3>& b)
{
    return sqrt(pow(a[0] - b[0], 2) + pow(a[1] - b[1], 2) + pow(a[2] - b[2], 2));
}

bool CheckDistanceCondition(std::shared_ptr<GeometrySurface>& surface, const SubRegion& region, double K0)
{
    std::array<double, 3> pnt1 = surface->d0Function({region.p1.u, region.p1.v});
    std::array<double, 3> pnt2 = surface->d0Function({region.p2.u, region.p2.v});
    std::array<double, 3> pnt3 = surface->d0Function({region.p3.u, region.p3.v});
    std::array<double, 3> pnt4 = surface->d0Function({region.p4.u, region.p4.v});

    std::array<double, 3> centerpoint_estimate{(pnt1[0] + pnt2[0] + pnt3[0] + pnt4[0]) / 4.0,
                                               (pnt1[1] + pnt2[1] + pnt3[1] + pnt4[1]) / 4.0,
                                               (pnt1[2] + pnt2[2] + pnt3[2] + pnt4[2]) / 4.0};

    double centerU = (region.p1.u + region.p4.u) / 2.0;
    double centerV = (region.p1.v + region.p4.v) / 2.0;

    std::array<double, 3> centerpoint_real = surface->d0Function({centerU, centerV});

    double distance = Distance(centerpoint_estimate, centerpoint_real);
    return distance <= K0;
}

bool CheckFirstDerivativeCondition(std::shared_ptr<GeometrySurface>& surface, const SubRegion& region, double K1)
{
    // 计算中心点的一阶导数
    double centerU = (region.p1.u + region.p4.u) / 2.0;
    double centerV = (region.p1.v + region.p4.v) / 2.0;
    std::array<double, 2> uv = {centerU, centerV};
    std::array<double, 3> d1u_center, d1v_center;
    surface->d1Function(uv, d1u_center, d1v_center);

    Eigen::Vector3d d1u_c(d1u_center.data());
    Eigen::Vector3d d1v_c(d1v_center.data());
    Eigen::Vector3d d1_c = d1u_c + d1v_c;
    // 对每个角点进行检查
    const SamplePoint* points[] = {&region.p1, &region.p2, &region.p3, &region.p4};
    for( const SamplePoint* point : points )
    {
        std::array<double, 2> uv_p = {point->u, point->v};
        std::array<double, 3> d1u_p, d1v_p;
        surface->d1Function(uv_p, d1u_p, d1v_p);
        Eigen::Vector3d d1u_p_c(d1u_p.data());
        Eigen::Vector3d d1v_p_c(d1v_p.data());
        Eigen::Vector3d d1_p = d1u_p_c + d1v_p_c;
        double dot_product = d1_c.dot(d1_p);
        double norm_product = d1_c.norm() * d1_p.norm();
        double angle = acos(dot_product / norm_product);
        double angle_degrees = angle * 180.0 / M_PI;
        if( angle_degrees > K1 )
            return false;
    }
    return true;
}
bool CheckFirstDerivativeCondition_UUU(std::shared_ptr<GeometrySurface>& surface, const SubRegion& region, double K1)
{
    // 计算中心点的一阶导数
    double centerU = (region.p1.u + region.p4.u) / 2.0;
    double centerV = (region.p1.v + region.p4.v) / 2.0;
    std::array<double, 2> uv = {centerU, centerV};
    std::array<double, 3> d1u_center, d1v_center;
    surface->d1Function(uv, d1u_center, d1v_center);

    Eigen::Vector3d d1u_c(d1u_center.data());
    Eigen::Vector3d d1_c = d1u_c;
    // 对每个角点进行检查
    const SamplePoint* points[] = {&region.p1, &region.p2, &region.p3, &region.p4};
    for( const SamplePoint* point : points )
    {
        std::array<double, 2> uv_p = {point->u, point->v};
        std::array<double, 3> d1u_p, d1v_p;
        surface->d1Function(uv_p, d1u_p, d1v_p);
        Eigen::Vector3d d1u_p_c(d1u_p.data());
        Eigen::Vector3d d1_p = d1u_p_c;
        double dot_product = d1_c.dot(d1_p);
        double norm_product = d1_c.norm() * d1_p.norm();
        double angle = acos(dot_product / norm_product);
        double angle_degrees = angle * 180.0 / M_PI;
        if( angle_degrees > K1 )
            return false;
    }
    return true;
}

bool CheckFirstDerivativeCondition_VVV(std::shared_ptr<GeometrySurface>& surface, const SubRegion& region, double K1)
{
    // 计算中心点的一阶导数
    double centerU = (region.p1.u + region.p4.u) / 2.0;
    double centerV = (region.p1.v + region.p4.v) / 2.0;
    std::array<double, 2> uv = {centerU, centerV};
    std::array<double, 3> d1u_center, d1v_center;
    surface->d1Function(uv, d1u_center, d1v_center);

    Eigen::Vector3d d1v_c(d1v_center.data());
    Eigen::Vector3d d1_c = d1v_c;
    // 对每个角点进行检查
    const SamplePoint* points[] = {&region.p1, &region.p2, &region.p3, &region.p4};
    for( const SamplePoint* point : points )
    {
        std::array<double, 2> uv_p = {point->u, point->v};
        std::array<double, 3> d1u_p, d1v_p;
        surface->d1Function(uv_p, d1u_p, d1v_p);
        Eigen::Vector3d d1v_p_c(d1v_p.data());
        Eigen::Vector3d d1_p = d1v_p_c;
        double dot_product = d1_c.dot(d1_p);
        double norm_product = d1_c.norm() * d1_p.norm();
        double angle = acos(dot_product / norm_product);
        double angle_degrees = angle * 180.0 / M_PI;
        if( angle_degrees > K1 )
            return false;
    }
    return true;
}

bool CheckCondition(std::shared_ptr<GeometrySurface>& surface, const SubRegion& region, double K0, double K1, double K2)
{
    // 检查 K0 条件
    if( !CheckDistanceCondition(surface, region, K0) )
    {
        return false;
    }

    // 检查 K1 条件
    if( !CheckFirstDerivativeCondition(surface, region, K1) )
    {
        return false;
    }

    // // 检查 K2 条件
    // if (!CheckSecondDerivativeCondition(surface, region, K2)) {
    // 	return false;
    // }

    return true;
}

void RefineSubRegion(std::shared_ptr<GeometrySurface>& surface, const SubRegion& region, double K0, double K1,
                     double K2, std::map<SamplePoint, bool>& sampledPoints, std::vector<SubRegion>& refinedRegions,
                     int depth, int maxDepth)
{

    if( depth >= maxDepth )
    {
        // 如果达到最大递归深度或满足条件，则不再细分
        refinedRegions.push_back(region);
        return;
    }

    // 计算0阶距离
    if( !CheckDistanceCondition(surface, region, K0) )
    {
        // 计算中点
        SamplePoint midP1P2 = {(region.p1.u + region.p2.u) / 2, (region.p1.v + region.p2.v) / 2};
        SamplePoint midP1P3 = {(region.p1.u + region.p3.u) / 2, (region.p1.v + region.p3.v) / 2};
        SamplePoint center = {(region.p1.u + region.p4.u) / 2, (region.p1.v + region.p4.v) / 2}; // p1和p4的中点
        SamplePoint midP2P4 = {(region.p2.u + region.p4.u) / 2, (region.p2.v + region.p4.v) / 2};
        SamplePoint midP3P4 = {(region.p3.u + region.p4.u) / 2, (region.p3.v + region.p4.v) / 2};

        // 添加新的采样点
        sampledPoints[midP1P2] = true;
        sampledPoints[midP1P3] = true;
        sampledPoints[center] = true;
        sampledPoints[midP2P4] = true;
        sampledPoints[midP3P4] = true;

        // 创建四个新的子区域
        SubRegion subRegion1 = {region.p1, midP1P2, midP1P3, center};
        SubRegion subRegion2 = {midP1P2, region.p2, center, midP2P4};
        SubRegion subRegion3 = {midP1P3, center, region.p3, midP3P4};
        SubRegion subRegion4 = {center, midP2P4, midP3P4, region.p4};

        // 递归调用细化函数
        RefineSubRegion(surface, subRegion1, K0, K1, K2, sampledPoints, refinedRegions, depth + 1, maxDepth);
        RefineSubRegion(surface, subRegion2, K0, K1, K2, sampledPoints, refinedRegions, depth + 1, maxDepth);
        RefineSubRegion(surface, subRegion3, K0, K1, K2, sampledPoints, refinedRegions, depth + 1, maxDepth);
        RefineSubRegion(surface, subRegion4, K0, K1, K2, sampledPoints, refinedRegions, depth + 1, maxDepth);
        return;
    }

    // 计算1阶距离u方向的
    if( !CheckFirstDerivativeCondition_UUU(surface, region, K1) )
    {
        // 计算中点
        SamplePoint midP1P2 = {(region.p1.u + region.p2.u) / 2, (region.p1.v + region.p2.v) / 2};
        SamplePoint midP3P4 = {(region.p3.u + region.p4.u) / 2, (region.p3.v + region.p4.v) / 2};

        // 添加新的采样点
        sampledPoints[midP1P2] = true;
        sampledPoints[midP3P4] = true;

        // 创建四个新的子区域
        SubRegion subRegion1 = {region.p1, midP1P2, region.p3, midP3P4};
        SubRegion subRegion2 = {midP1P2, region.p2, midP3P4, region.p4};

        // 递归调用细化函数
        RefineSubRegion(surface, subRegion1, K0, K1, K2, sampledPoints, refinedRegions, depth + 1, maxDepth);
        RefineSubRegion(surface, subRegion2, K0, K1, K2, sampledPoints, refinedRegions, depth + 1, maxDepth);
        return;
    }

    if( !CheckFirstDerivativeCondition_VVV(surface, region, K0) )
    {
        // 计算中点
        SamplePoint midP1P3 = {(region.p1.u + region.p3.u) / 2, (region.p1.v + region.p3.v) / 2};
        SamplePoint midP2P4 = {(region.p2.u + region.p4.u) / 2, (region.p2.v + region.p4.v) / 2};

        // 添加新的采样点
        sampledPoints[midP1P3] = true;
        sampledPoints[midP2P4] = true;

        // 创建四个新的子区域
        SubRegion subRegion1 = {region.p1, region.p2, midP1P3, midP2P4};
        SubRegion subRegion2 = {midP1P3, midP2P4, region.p3, region.p4};

        // 递归调用细化函数
        RefineSubRegion(surface, subRegion1, K0, K1, K2, sampledPoints, refinedRegions, depth + 1, maxDepth);
        RefineSubRegion(surface, subRegion2, K0, K1, K2, sampledPoints, refinedRegions, depth + 1, maxDepth);
        return;
    }
    // 计算中点
    SamplePoint midP1P2 = {(region.p1.u + region.p2.u) / 2, (region.p1.v + region.p2.v) / 2};
    SamplePoint midP1P3 = {(region.p1.u + region.p3.u) / 2, (region.p1.v + region.p3.v) / 2};
    SamplePoint center = {(region.p1.u + region.p4.u) / 2, (region.p1.v + region.p4.v) / 2}; // p1和p4的中点
    SamplePoint midP2P4 = {(region.p2.u + region.p4.u) / 2, (region.p2.v + region.p4.v) / 2};
    SamplePoint midP3P4 = {(region.p3.u + region.p4.u) / 2, (region.p3.v + region.p4.v) / 2};

    // 添加新的采样点
    sampledPoints[midP1P2] = true;
    sampledPoints[midP1P3] = true;
    sampledPoints[center] = true;
    sampledPoints[midP2P4] = true;
    sampledPoints[midP3P4] = true;

    return;
}

bool CheckDistanceConditionCurve(std::shared_ptr<TiGER::GeometryCurve> curve, const SubRegionCurve& region, double K0)
{
    std::array<double, 3> a = curve->d0Function(region.p1.u);
    std::array<double, 3> b = curve->d0Function(region.p2.u);
    std::array<double, 3> centerpoint_estimate = {(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2};

    // 计算实际的中心点（曲线上对应的点）
    double midU = (region.p1.u + region.p2.u) / 2.0;
    std::array<double, 3> centerpoint_real = curve->d0Function(midU);

    // 计算两点间的欧几里得距离
    double distance = Distance(centerpoint_estimate, centerpoint_real);
    // std::cout << "distance d0:" <<distance << std::endl;
    //  判断距离是否满足条件
    return distance <= K0;
}

bool CheckFirstDerivativeConditionCurve(std::shared_ptr<TiGER::GeometryCurve> curve, const SubRegionCurve& region,
                                        double K1)
{
    // 计算中心点的一阶导数
    double midU = (region.p1.u + region.p2.u) / 2.0;
    std::array<double, 3> centerPnt;
    curve->d1Function(midU, centerPnt);
    Eigen::Vector3d vec2(centerPnt.data());

    std::vector<SamplePointCurve> points;
    // 添加端点
    points.push_back(region.p1);
    points.push_back(region.p2);

    for( const SamplePointCurve point : points )
    {
        std::array<double, 3> temp;
        curve->d1Function(point.u, temp);
        Eigen::Vector3d vec1(temp.data());
        // 计算两个向量的夹角
        double dot_product = vec1.dot(vec2);
        double norm_product = vec1.norm() * vec2.norm();
        double angle = acos(dot_product / norm_product);
        double angle_degrees = angle * 180.0 / M_PI;
        if( angle_degrees > K1 )
            return false;
    }
    return true;
}

bool CheckConditionCurve(std::shared_ptr<TiGER::GeometryCurve> curve, const SubRegionCurve& region, double K0,
                         double K1, double K2)
{
    // 检查 K0 条件
    if( !CheckDistanceConditionCurve(curve, region, K0) )
    {
        return false;
    }

    // 检查 K1 条件
    if( !CheckFirstDerivativeConditionCurve(curve, region, K1) )
    {
        return false;
    }

    return true;
}

void RefineSubRegionCurve(std::shared_ptr<TiGER::GeometryCurve> curve, const SubRegionCurve& region, double K0,
                          double K1, double K2, std::map<double, bool>& sampledPoints, int depth, int maxDepth)
{
    if( depth >= maxDepth || CheckConditionCurve(curve, region, K0, K1, K2) )
    {
        return;
    }
    // 不满足条件，继续细化
    double midU = (region.p1.u + region.p2.u) / 2.0;
    SamplePointCurve midPoint = {midU};

    // 创建两个新的子区域
    SubRegionCurve subRegion1 = {region.p1, midPoint};
    SubRegionCurve subRegion2 = {midPoint, region.p2};

    // 将新的中点加入采样点集
    sampledPoints[midU] = true;

    // 递归细化新的子区域
    RefineSubRegionCurve(curve, subRegion1, K0, K1, K2, sampledPoints, depth + 1, maxDepth);
    RefineSubRegionCurve(curve, subRegion2, K0, K1, K2, sampledPoints, depth + 1, maxDepth);

    return;
}

void CreateSubRegionsCurve(const std::vector<double>& u_values, std::vector<SubRegionCurve>& subregions)
{
    subregions.clear(); // 清除现有内容

    for( size_t i = 0; i < u_values.size() - 1; ++i )
    {
        SubRegionCurve region;
        region.p1 = {u_values[i]};
        region.p2 = {u_values[i + 1]};
        subregions.push_back(region);
    }
}

} // namespace help


void discreteCurve_Binary(std::shared_ptr<TiGER::GeometryCurve> curve, const TiGER::CurveBParameters& args,
                          TiGER::LineMesh& segments)
{
    std::map<double, bool> sampledPoints;
    // 均匀布点，构建初始的二叉树区间，可以改成曲率自适应
    std::array<double, 2> curveUV = curve->getUVScaleFunction();
    std::vector<double> u_values;
    std::vector<TiGER::help::SubRegionCurve> subregions;
    int n = 10;
    double delta = (curveUV[1] - curveUV[0]) / n;
    for( int i = 0; i <= n; i++ )
    {
        double para = curveUV[0] + delta * i;
        if( i == 0 )
        {
            para += delta / 10;
        }
        else if( i == n )
        {
            para -= delta / 10;
        }
        u_values.emplace_back(para);
    }
    // 二叉树细分
    TiGER::help::CreateSubRegionsCurve(u_values, subregions);
    for( auto u : u_values )
    {
        sampledPoints[u] = true;
    }
    for( const auto region : subregions )
    {
        TiGER::help::RefineSubRegionCurve(curve, region, args.getDeflection(), args.getAngle(), 5, sampledPoints, 0,
                                          10);
    }
    for( auto p : sampledPoints )
    {
        double u = p.first;
        std::array<double, 3> temp_p = curve->d0Function(u);
        segments.coord.push_back(temp_p);
    }
    for( int i = 1; i < segments.coord.size(); i++ )
    {
        std::array<int, 2> temp_s = {i - 1, i};
        segments.segments.push_back(temp_s);
    }
}
// 计算二维向量的叉积
double cross(const std::array<double, 2>& v1, const std::array<double, 2>& v2)
{
    return v1[0] * v2[1] - v1[1] * v2[0];
}

// 检查两条线段是否相交
bool doSegmentsIntersect(const std::array<double, 2>& p1, const std::array<double, 2>& p2,
                         const std::array<double, 2>& q1, const std::array<double, 2>& q2)
{
    std::array<double, 2> u = {p2[0] - p1[0], p2[1] - p1[1]};
    std::array<double, 2> v = {q2[0] - q1[0], q2[1] - q1[1]};
    std::array<double, 2> w1 = {q1[0] - p1[0], q1[1] - p1[1]};
    std::array<double, 2> w2 = {q2[0] - p1[0], q2[1] - p1[1]};

    double cross1 = cross(w1, u);
    double cross2 = cross(w2, u);
    double cross3 = cross({p1[0] - q1[0], p1[1] - q1[1]}, v);
    double cross4 = cross({p2[0] - q1[0], p2[1] - q1[1]}, v);

    return cross1 * cross2 < 0 && cross3 * cross4 < 0;
}

// 检查线网是否有自相交的线段
bool hasSelfIntersections(const LineMesh& mesh)
{
    const auto& segments = mesh.segments;
    const auto& coord = mesh.coord;

    for( size_t i = 0; i < segments.size(); ++i )
    {
        const auto& seg1 = segments[i];
        std::array<double, 2> p1 = {coord[seg1[0]][0], coord[seg1[0]][1]};
        std::array<double, 2> p2 = {coord[seg1[1]][0], coord[seg1[1]][1]};

        for( size_t j = i + 1; j < segments.size(); ++j )
        {
            const auto& seg2 = segments[j];
            std::array<double, 2> q1 = {coord[seg2[0]][0], coord[seg2[0]][1]};
            std::array<double, 2> q2 = {coord[seg2[1]][0], coord[seg2[1]][1]};

            if( doSegmentsIntersect(p1, p2, q1, q2) )
            {
                return true;
            }
        }
    }
    return false;
}

bool isEar(int i, const std::vector<int>& indices, const std::vector<std::array<double, 3>>& points)
{
    int prev = indices[(i - 1 + indices.size()) % indices.size()];
    int curr = indices[i];
    int next = indices[(i + 1) % indices.size()];
    std::array<double, 3> p1 = points[prev], p2 = points[curr], p3 = points[next];

    // 计算面积，判断是否是逆时针
    double area = (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[0] - p1[0]);
    if( area <= 0 )
        return false;

    // 检查其他点是否在三角形内部
    for( int k = 0; k < points.size(); ++k )
    {
        if( k == prev || k == curr || k == next )
            continue;
        std::array<double, 3> pt = points[k];
        double u = ((pt[1] - p3[1]) * (p2[0] - p3[0]) - (pt[0] - p3[0]) * (p2[1] - p3[1])) / area;
        double v = ((pt[1] - p1[1]) * (p3[0] - p1[0]) - (pt[0] - p1[0]) * (p3[1] - p1[1])) / area;
        double w = 1 - u - v;
        if( u > 0 && v > 0 && w > 0 )
            return false;
    }
    return true;
}

void earClipping(const std::vector<std::array<double, 3>>& points, SurfaceMesh& mesh)
{
    int n = points.size();
    if( n < 3 )
        return;

    std::vector<int> indices(n);
    for( int i = 0; i < n; ++i )
        indices[i] = i;

    while( indices.size() > 3 )
    {
        bool earFound = false;
        for( int i = 0; i < indices.size(); ++i )
        {
            if( isEar(i, indices, points) )
            {
                int prev = indices[(i - 1 + indices.size()) % indices.size()];
                int curr = indices[i];
                int next = indices[(i + 1) % indices.size()];

                std::array<int, 3> tri = {prev, curr, next};
                mesh.tris.push_back(tri);
                indices.erase(indices.begin() + i);
                earFound = true;
                break;
            }
        }
        if( !earFound )
        {
            std::cerr << "No ear found. The input might not be a simple polygon." << std::endl;
            return;
        }
    }

    std::array<int, 3> lastTri = {indices[0], indices[1], indices[2]};
    mesh.tris.push_back(lastTri);
}

void discreteSurface_Binary(std::shared_ptr<GeometrySurface>& surface, const std::vector<LineMesh>& boundary_segments,
                            const TiGER::generateStlDTParameters& arg_stldt, SurfaceMesh& output_mesh)
{
    // 预处理
    std::vector<std::vector<std::array<double, 3>>> boundary_points;
    for( auto bs : boundary_segments )
        boundary_points.push_back(bs.coord);
    TiGER::help::dAssemblecurve(boundary_points);

    LineMesh Loop_in;
    std::array<double, 4> uv_scale = surface->getUVScaleFunction();
    std::array<double, 2> hint = {(uv_scale[1] - uv_scale[0]) / 2, (uv_scale[3] - uv_scale[2]) / 2};
    std::map<std::array<double, 3>, int> mp;
    for( auto bs : boundary_points )
    {
        std::vector<int> coord_idx;
        for( auto co : bs )
        {
            int n;
            if( mp.count(co) != 0 )
            {
                n = mp[co];
            }
            else
            {
                n = mp.size();
                mp[co] = n;
                std::array<double, 2> res_uv;
                double eps=1e-3;
                surface->projectFunction(co, eps, hint, res_uv);
                std::array<double, 3> pt = {res_uv[0], res_uv[1], 0};
                Loop_in.coord.emplace_back(pt);
            }
            coord_idx.push_back(n);
        }
        int nn = coord_idx.size();
        if( nn >= 2 )
        {
            for( int i = 1; i < nn; i++ )
            {
                std::array<int, 2> a = {coord_idx[i - 1], coord_idx[i]};
                Loop_in.segments.push_back(a);
            }
        }
    }
    if( hasSelfIntersections(Loop_in) )
    {
        cout << "!!!!!!self intersection!!!!!" << endl;
        mp.clear();
        surface->SmoothFace();
        Loop_in.coord.clear();
        Loop_in.segments.clear();
        for( auto bs : boundary_points )
        {
            std::array<double, 2> last_uv = hint;
            std::vector<int> coord_idx;
            for( auto co : bs )
            {
                int n;
                if( mp.count(co) != 0 )
                {
                    n = mp[co];
                }
                else
                {
                    double eps;
                    std::array<double, 2> res_uv = last_uv;
                    surface->projectFunction(co, eps, hint, res_uv);
                    // TiGER::projecthelp::sf_project(surface, co, 1e-6, res_uv);
                    if( std::isnan(res_uv[0]) || std::isnan(res_uv[1]) )
                    {
                        surface->projectFunction(co, eps, hint, res_uv);
                        // continue;
                    }
                    n = mp.size();
                    mp[co] = n;
                    std::array<double, 3> pt = {res_uv[0], res_uv[1], 0};
                    Loop_in.coord.emplace_back(pt);
                    // last_uv = res_uv;
                }
                coord_idx.push_back(n);
            }
            int nn = coord_idx.size();
            if( nn >= 2 )
            {
                for( int i = 1; i < nn; i++ )
                {
                    if( coord_idx[i - 1] == coord_idx[i] )
                        continue;
                    std::array<int, 2> a = {coord_idx[i - 1], coord_idx[i]};
                    Loop_in.segments.push_back(a);
                }
            }
        }
    }

    // 四叉树布点
    double umin, umax, vmin, vmax;
    TiGER::help::computeBoundingBox(Loop_in, umin, umax, vmin, vmax);
    TiGER::help::SamplePoint p1 = {umin, vmin};
    TiGER::help::SamplePoint p2 = {umax, vmin};
    TiGER::help::SamplePoint p3 = {umin, vmax};
    TiGER::help::SamplePoint p4 = {umax, vmax};
    TiGER::help::SubRegion initRegion = {p1, p2, p3, p4};
    std::map<TiGER::help::SamplePoint, bool> sampledPoints;
    std::vector<TiGER::help::SubRegion> refinedRegions;

    RefineSubRegion(surface, initRegion, arg_stldt.getDeflection(), arg_stldt.getAngle(), 0.05, sampledPoints,
                    refinedRegions, 0, arg_stldt.getLayer());

    double threshold = 1e-6;
    cout << "111111111111" << endl;
    cout << sampledPoints.size() << endl;
    cout << "111111111111" << endl;
    for( auto p : sampledPoints )
    {
        auto point_to_add = p.first;
        std::array<double, 3> temp = {point_to_add.u, point_to_add.v, 0};
        int flag = 1;
        for( const auto& curve : Loop_in.segments )
        {
            TiGER::help::SamplePoint p = {point_to_add.u, point_to_add.v};
            TiGER::help::SamplePoint p1 = {Loop_in.coord[curve[0]][0], Loop_in.coord[curve[0]][1]};
            TiGER::help::SamplePoint p2 = {Loop_in.coord[curve[1]][0], Loop_in.coord[curve[1]][1]};

            double d_temp = TiGER::help::ProjectionDistancePointToSegment(p, p1, p2);
            if( d_temp <= threshold )
            {
                flag = 0;
                break;
            }
        }
        if( flag )
        {
            Loop_in.coord.push_back(temp);
        }
    }

    // dt
    TriangleDealunay args_dt;
    args_dt.setsizealpha(arg_stldt.getsizealpha());
    std::function<double(const std::array<double, 3>)> fun3 = [](const std::array<double, 3> a) { return 5.0; };
    discretize_surface_delaunay::TMeshSurf_processTriangleLoop(Loop_in, args_dt, fun3, output_mesh);

    for( auto& pt : output_mesh.coords )
    {
        std::array<double, 3> pt_3d;
        std::array<double, 2> pt_2d = {pt[0], pt[1]};
        pt_3d = surface->d0Function(pt_2d);
        pt = pt_3d;
    }
}

bool saveShapeToSTEP(const std::string& name, const TopoDS_Shape& shape)
{
    STEPControl_Writer writer;
    IFSelect_ReturnStatus status = writer.Transfer(shape, STEPControl_AsIs);

    if( status != IFSelect_RetDone )
    {
        std::cerr << "Error transferring shape to STEP format." << std::endl;
        return false;
    }

    if( writer.Write(name.c_str()) != IFSelect_RetDone )
    {
        std::cerr << "Error writing shape to file: " << name << std::endl;
        return false;
    }

    return true;
}


void GeometryRepairOCC(TiGER::Geometry& geometry,bool divideclosed_flag,int heal,bool bool_flag){
    
}


void CAD_Tessellation(TiGER::Geometry& geometry, TiGER::SurfaceMesh& surfaceOut, TiGER::generateStlParameters& arg_stl)
{
    int currentPointCount = 0;
    std::vector<int> faceids;
    int faceid = 0;
    for( auto tempf : geometry.surfaces_ )
    {
        std::vector<double> bnd_pts_t;
        std::vector<int> bnd_triangles_t;
        tempf->discreteFaceFunction(arg_stl.getDeflection(), arg_stl.getAngle(), arg_stl.getRelative(),
                                    arg_stl.getInParallel(), arg_stl.getMinSize(), arg_stl.getInternalVerticesMode(),
                                    arg_stl.getControlSurfaceDeflection());
        tempf->getBndFunction(bnd_pts_t, bnd_triangles_t);
        if( bnd_triangles_t.size() == 0 )
        {
            cout << "occ discrete error" << endl;
        }

        std::vector<std::array<double, 3>> localPoints;
        for( size_t i = 0; i < bnd_pts_t.size(); i += 3 )
        {
            localPoints.push_back({bnd_pts_t[i], bnd_pts_t[i + 1], bnd_pts_t[i + 2]});
        }
        surfaceOut.coords.insert(surfaceOut.coords.end(), localPoints.begin(), localPoints.end());

        // 处理三角形索引的转换
        for( size_t i = 0; i < bnd_triangles_t.size(); i += 3 )
        {
            surfaceOut.tris.push_back({bnd_triangles_t[i] + currentPointCount,
                                       bnd_triangles_t[i + 1] + currentPointCount,
                                       bnd_triangles_t[i + 2] + currentPointCount});
            surfaceOut.regions.push_back(faceid);
        }
        currentPointCount += localPoints.size();
        faceid++;
    }
}

void VCAD_Tessellation(TiGER::Geometry& geometry, TiGER::SurfaceMesh& surfaceOut,
                        TiGER::generateStlDTParameters& arg_stldt)
{
    // 将所有的线都离散了,存储为linemesh
    std::vector<LineMesh> segements;
    CurveBParameters args;
    args.setDeflection(arg_stldt.getDeflection());
    args.setAngle(arg_stldt.getAngle());
    args.setLayer(arg_stldt.getLayer());

    int dimension = arg_stldt.getdemension();
    double delta = arg_stldt.getdelta();
    for( int i = 0; i < geometry.curves_.size(); i++ )
    {
        auto curve = geometry.curves_[i];
        LineMesh segment;
        // discreteCurve_Binary(curve, args, segment);
        TiGER::discretecurcestl::discreteCurve_uniform_test(curve, dimension, delta, segment);
        segements.push_back(segment);
    }

    // 对每一个面根据拓扑分配离散后的线网格
    std::vector<std::vector<LineMesh>> sf_boundary;
    sf_boundary.resize(geometry.surfaces_.size());
    for( int i = 0; i < geometry.topo_surf_to_curves_.size(); i++ )
    {
        for( int j = 0; j < geometry.topo_surf_to_curves_[i].size(); j++ )
        {
            sf_boundary[i].push_back(segements[geometry.topo_surf_to_curves_[i][j]]);
        }
    }

    // 离散每一个面
    int currentPointCount = 0;
    for( int faceid = 0; faceid < geometry.surfaces_.size(); faceid++ )
    // for( int faceid = 47; faceid <  48; faceid++ )
    {
        cout << "faceid:" << faceid << endl;
        SurfaceMesh sf_mesh;
        discreteSurface_Binary(geometry.surfaces_[faceid], sf_boundary[faceid], arg_stldt, sf_mesh);

        surfaceOut.coords.insert(surfaceOut.coords.end(), sf_mesh.coords.begin(), sf_mesh.coords.end());
        for( int i = 0; i < sf_mesh.tris.size(); i++ )
        {
            surfaceOut.tris.push_back({
                sf_mesh.tris[i][0] + currentPointCount,
                sf_mesh.tris[i][1] + currentPointCount,
                sf_mesh.tris[i][2] + currentPointCount,
            });
            surfaceOut.regions.push_back(faceid);
        }
        currentPointCount += sf_mesh.coords.size();
    }
}

void mapReal2Normal_Curve(std::shared_ptr<GeometryCurve> curve, double& real, double& normal)
{
    std::array<double, 2> para = curve->getUVScaleFunction();
    double umin = para[0];
    double umax = para[1];
    normal = (real - umin) / (umax - umin);
}

void mapNormal2Real_Curve(std::shared_ptr<GeometryCurve> curve, double& real, double& normal)
{
    std::array<double, 2> para = curve->getUVScaleFunction();
    double umin = para[0];
    double umax = para[1];
    real = umin + normal * (umax - umin);
}

std::shared_ptr<GeometryCurve> trimmedCurve(std::shared_ptr<GeometryCurve> curve, double& u1, double& u2)
{
    // 采样点的数量
    int sampleCount = 23;
    std::vector<gp_Pnt> sampledPoints;

    // 采样从 u1 到 u2 之间的点
    for( int i = 0; i < sampleCount; ++i )
    {
        double param = u1 + i * (u2 - u1) / (sampleCount - 1);

        std::array<double, 3> aaa = curve->d0Function(param);
        gp_Pnt sampledPoint{aaa[0], aaa[1], aaa[2]};
        sampledPoints.push_back(sampledPoint);
    }

    // 将采样点转换为 TColgp_Array1OfPnt
    TColgp_Array1OfPnt pointsArray(1, sampleCount);
    for( int i = 0; i < sampleCount; ++i )
    {
        pointsArray.SetValue(i + 1, sampledPoints[i]);
    }

    // 使用 GeomAPI_PointsToBSpline 拟合新的 B 样条曲线
    Handle(Geom_BSplineCurve) fittedCurve = GeomAPI_PointsToBSpline(pointsArray).Curve();
    TopoDS_Edge topo_edge = BRepBuilderAPI_MakeEdge(fittedCurve);
    BRepAdaptor_Curve curveout(topo_edge);
    std::shared_ptr<GeometryCurve_OCC> curve_temp = std::make_shared<GeometryCurve_OCC>(topo_edge, curveout);
    return curve_temp;
}

enum CurveRelationship
{
    Unrelated,
    Related
};

CurveRelationship checkCurveRelationship(int first_edge_id, int second_edge_id,
                                         std::vector<std::shared_ptr<GeometryCurve>>&curves_)
{
    std::shared_ptr<GeometryCurve> c1 = curves_[first_edge_id];
    std::shared_ptr<GeometryCurve> c2 = curves_[second_edge_id];
    
    int numSamples = 10;
    std::array<double, 2> uv1 = c1->getUVScaleFunction();
    std::array<double, 2> uv2 = c2->getUVScaleFunction();

    std::array<double, 3> tpoint = c1->d0Function(uv1[0]);
    double first_curve_min_x = tpoint[0];
    double first_curve_max_x = tpoint[0];
    double first_curve_min_y = tpoint[1];
    double first_curve_max_y = tpoint[1];
    double first_curve_min_z = tpoint[2];
    double first_curve_max_z = tpoint[2];
    double step = (uv1[1] - uv1[0]) / numSamples;
    for( int i = 0; i <= numSamples; i++ )
    {
        double para = uv1[0] + step * i;
        std::array<double, 3> point = c1->d0Function(para);
        if( point[0] < first_curve_min_x )
            first_curve_min_x = point[0];
        if( point[0] > first_curve_max_x )
            first_curve_max_x = point[0];
        if( point[1] < first_curve_min_y )
            first_curve_min_y = point[1];
        if( point[1] > first_curve_max_y )
            first_curve_max_y = point[1];
        if( point[2] < first_curve_min_z )
            first_curve_min_z = point[2];
        if( point[2] > first_curve_max_z )
            first_curve_max_z = point[2];
    }
    const double distance_threshold = 1e-6;
    int fflag = 1;
    step = (uv2[1] - uv2[0]) / numSamples;
    for( int i = 0; i <= numSamples; i++ )
    {
        double para = uv1[0] + step * i;
        std::array<double, 3> point = c2->d0Function(para);
        if( point[0] < first_curve_min_x || point[0] > first_curve_max_x || point[1] < first_curve_min_y ||
            point[1] > first_curve_max_y || point[2] < first_curve_min_z || point[2] > first_curve_max_z )
        {
        }
        else
        {
            fflag = 0;
            break;
        }
    }
    if( fflag )
    {
        return Unrelated;
    }

    fflag = 1;
    for( int i = 0; i <= numSamples; i++ )
    {
        double para = uv1[0] + step * i;
        std::array<double, 3> point = c2->d0Function(para);
        double tempe;
        double hint;
        double u;
        c1->projectFunction_dis(point, 1e-3, hint, u,tempe);
        if( tempe > distance_threshold )
        {
            fflag = 0;
            break;
        }
    }
    if( fflag )
    {
        return Related;
    }
    return Unrelated;
}

void replaceCurveIds(std::vector<std::vector<int>>& topo_surf_to_curves, int old_id, const std::vector<int>& new_ids)
{
    for( auto& curves : topo_surf_to_curves )
    {
        for( auto it = curves.begin(); it != curves.end(); )
        {
            if( *it == old_id )
            {
                it = curves.erase(it); // 移除旧的曲线ID
                // 插入新的曲线ID
                it = curves.insert(it, new_ids.begin(), new_ids.end());
                // 更新迭代器位置
                it += new_ids.size();
            }
            else
            {
                ++it; // 移动到下一个元素
            }
        }
    }
}

vector<vector<int>> deal_complex_edge(int first_edge_id, int second_edge_id, vector<vector<int>>& topo_surf_to_curves,
                                      std::vector<std::shared_ptr<GeometryCurve>>& curves_)
{
    vector<int> newedge;
    vector<int> deledge;
    newedge.clear();
    deledge.clear();
    vector<vector<int>> result;
    std::shared_ptr<GeometryCurve> c1 = curves_[first_edge_id];
    std::shared_ptr<GeometryCurve> c2 = curves_[second_edge_id];
    std::array<double, 2> uv1 = c1->getUVScaleFunction();
    std::array<double, 2> uv2 = c2->getUVScaleFunction();

    std::array<double, 3> startpoint1 = c1->d0Function(uv1[0]);
    std::array<double, 3> endpoint1 = c1->d0Function(uv1[1]);

    double num0, num1, num2, num3;
    double eps, hint;
    c2->projectFunction(startpoint1, eps, hint, num0);
    c2->projectFunction(endpoint1, eps, hint, num1);
    num2 = uv2[0];
    num3 = uv2[1];
    double therehold = 1e-6;
    // 考虑重复面无关曲线情况
    if( fabs(num0 - num1) < therehold )
    {
        result.push_back(newedge);
        result.push_back(deledge);
        return result;
    }

    // 考虑重复面无关曲线情况
    if( fabs(num2 - num3) < therehold )
    {
        result.push_back(newedge);
        result.push_back(deledge);
        return result;
    }
    if( num0>num3 )
    {
        result.push_back(newedge);
        result.push_back(deledge);
        return result;
    }
    
    // 确保曲线的一致性
    int swap01 = 0;
    if( num0 > num1 )
    {
        swap01 = 1;
        swap(num0, num1);
    }
    if( num2 > num3 )
    {
        swap(num2, num3);
    }

    //去掉不合理情况
    if( num0 - num3 >-1e-6 )
    {
        result.push_back(newedge);
        result.push_back(deledge);
        return result;
    }
    if( num2 - num1 >-1e-6 )
    {
        result.push_back(newedge);
        result.push_back(deledge);
        return result;
    }

    //正式处理
    if( fabs(num0 - num2) < therehold )
    {
        if( fabs(num1 - num3) < therehold )
        {
            replaceCurveIds(topo_surf_to_curves, second_edge_id, {first_edge_id});
            deledge.push_back(second_edge_id);
            result.push_back(newedge);
            result.push_back(deledge);
            return result;
        }
        else if( num1 > num3 )
        {
            std::array<double, 3> sppoint = c2->d0Function(num3);
            double u;

            c1->projectFunction(sppoint, eps, hint, u);
            std::shared_ptr<GeometryCurve> nc2 = trimmedCurve(c2, uv1[0], u);
            int nc2_id = curves_.size();
            curves_.push_back(nc2);

            std::shared_ptr<GeometryCurve> nc3 = trimmedCurve(c2, u, uv1[1]);
            int nc3_id = curves_.size();
            curves_.push_back(nc3);
            if( swap01 )
            {
                replaceCurveIds(topo_surf_to_curves, first_edge_id, {nc2_id, nc3_id});
                replaceCurveIds(topo_surf_to_curves, second_edge_id, {nc3_id});
            }
            else
            {
                replaceCurveIds(topo_surf_to_curves, first_edge_id, {nc2_id, nc3_id});
                replaceCurveIds(topo_surf_to_curves, second_edge_id, {nc2_id});
            }
            newedge.push_back(nc2_id);
            newedge.push_back(nc3_id);
            deledge.push_back(first_edge_id);
            deledge.push_back(second_edge_id);
            result.push_back(newedge);
            result.push_back(deledge);
            return result;
        }
        else
        {
            std::shared_ptr<GeometryCurve> nc3 = trimmedCurve(c2, num1, num3);
            int nc3_id = curves_.size();
            curves_.push_back(nc3);
            replaceCurveIds(topo_surf_to_curves, second_edge_id, {first_edge_id, nc3_id});
            newedge.push_back(nc3_id);
            deledge.push_back(second_edge_id);
            result.push_back(newedge);
            result.push_back(deledge);
            return result;
        }
    }
    else if( num0 > num2 )
    {
        if( fabs(num1 - num3) < therehold )
        {
            std::shared_ptr<GeometryCurve> nc1 = trimmedCurve(c2, num2, num0);
            int nc1_id = curves_.size();
            curves_.push_back(nc1);
            replaceCurveIds(topo_surf_to_curves, second_edge_id, {nc1_id, first_edge_id});
            newedge.push_back(nc1_id);
            deledge.push_back(second_edge_id);
            result.push_back(newedge);
            result.push_back(deledge);
            return result;
        }
        else if( num1 > num3 )
        {
            std::shared_ptr<GeometryCurve> nc1 = trimmedCurve(c2, num2, num0);
            int nc1_id = curves_.size();
            curves_.push_back(nc1);

            std::array<double, 3> sppoint = c2->d0Function(num3);
            double u;
            
            c1->projectFunction(sppoint, eps, hint, u);
            std::shared_ptr<GeometryCurve> nc2 = trimmedCurve(c2, uv1[0], u);
            int nc2_id = curves_.size();
            curves_.push_back(nc2);

            std::shared_ptr<GeometryCurve> nc3 = trimmedCurve(c2, u, uv1[1]);
            int nc3_id = curves_.size();
            curves_.push_back(nc3);

            if( swap01 )
            {
                replaceCurveIds(topo_surf_to_curves, first_edge_id, {nc2_id, nc3_id});
                replaceCurveIds(topo_surf_to_curves, second_edge_id, {nc1_id, nc3_id});
            }
            else
            {
                replaceCurveIds(topo_surf_to_curves, first_edge_id, {nc2_id, nc3_id});
                replaceCurveIds(topo_surf_to_curves, second_edge_id, {nc1_id, nc2_id});
            }
            newedge.push_back(nc1_id);
            newedge.push_back(nc2_id);
            newedge.push_back(nc3_id);
            deledge.push_back(first_edge_id);
            deledge.push_back(second_edge_id);
            result.push_back(newedge);
            result.push_back(deledge);
            return result;
        }
        else
        {
            std::shared_ptr<GeometryCurve> nc1 = trimmedCurve(c2, num2, num0);
            std::shared_ptr<GeometryCurve> nc3 = trimmedCurve(c2, num1, num3);
            int nc1_id = curves_.size();
            curves_.push_back(nc1);
            int nc3_id = curves_.size();
            curves_.push_back(nc3);
            replaceCurveIds(topo_surf_to_curves, second_edge_id, {nc1_id, first_edge_id, nc3_id});
            newedge.push_back(nc1_id);
            newedge.push_back(nc3_id);
            deledge.push_back(second_edge_id);
            result.push_back(newedge);
            result.push_back(deledge);
            return result;
        }
    }
    else
    {
        return deal_complex_edge(second_edge_id, first_edge_id, topo_surf_to_curves, curves_);
        /*if( fabs(num1 - num3) < therehold )
        {
            std::array<double, 3> sppoint = c2->d0Function(num2);
            double u;

            c1->projectFunction(sppoint, eps, hint, u);
            std::shared_ptr<GeometryCurve> nc2 = trimmedCurve(c2, uv1[0], u);
            int nc2_id = curves_.size();
            curves_.push_back(nc2);

            std::shared_ptr<GeometryCurve> nc3 = trimmedCurve(c2, u, uv1[1]);
            int nc3_id = curves_.size();
            curves_.push_back(nc3);

            if( swap01 )
            {
                replaceCurveIds(topo_surf_to_curves, first_edge_id, {nc2_id, nc3_id});
                replaceCurveIds(topo_surf_to_curves, second_edge_id, {nc2_id});
            }
            else
            {
                replaceCurveIds(topo_surf_to_curves, first_edge_id, {nc2_id, nc3_id});
                replaceCurveIds(topo_surf_to_curves, second_edge_id, {nc3_id});
            }
        }
        else if( num1 > num3 )
        {
            
        
        }*/
    }

    return result;
}


void dealcomplexcurve(TiGER::Geometry& Geo)
{
    vector<int> pending;
    for( int i = 0; i < Geo.curves_.size(); i++ )
    {
        pending.push_back(i);
    }

    while( pending.size() >= 1 )
    {
        int first_curve = pending[0];
        int overlap_flag = 0;
        for( int i = 1; i < pending.size(); i++ )
        {
            int second_curve = pending[i];
            auto relationship = checkCurveRelationship(first_curve, second_curve, Geo.curves_);

            if( relationship == Related )
            {
                cout << "first_curve: " << first_curve << endl;
                cout << "second_curve: " << second_curve << endl;
                vector<vector<int>> result;
                result = deal_complex_edge(first_curve, second_curve,Geo.topo_surf_to_curves_,Geo.curves_);
                
                vector<int> newedge = result[0];
                vector<int> deledge = result[1];

                // 考虑到无关垂直的曲线
                if( deledge.size() == 0 )
                {
                    continue;
                }
                overlap_flag = 1;
                for( auto curve_idd : newedge )
                {
                    if( find(pending.begin(), pending.end(), curve_idd) == pending.end() )
                    {
                        pending.push_back(curve_idd);
                    }
                }
                for( auto curve_idd : deledge )
                {
                    auto it = std::find(pending.begin(), pending.end(), curve_idd);
                    pending.erase(it);
                }
                break;
            }
        }
        if( overlap_flag == 0 )
        {
            auto it = pending.begin();
            pending.erase(it);
        }
    }

}


} // namespace TiGER
