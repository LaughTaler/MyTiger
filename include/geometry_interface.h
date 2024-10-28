#pragma once
#include "geometry_curve_interface.h"
#include "geometry_surface_interface.h"
#include <vector>
#include <memory>
namespace SmeshGen

{
    class Geometry {
    public:
        /** 几何信息 */
        std::vector<std::shared_ptr<GeometryCurve>> curves_;
        std::vector<std::shared_ptr<GeometrySurface>> surfaces_;

        /** 拓扑信息 */
        std::vector<std::vector<int>> topo_surf_to_curves_; // 面到曲线的拓扑信息
        std::vector<std::vector<int>> topo_body_to_surfs_;  // 体到面的拓扑信息

        /**
         * @brief 读取几何到内存中
         *
         * @param[in] filename 读入文件名，支持 stp, igs, step, iges, sat, brep, x_t ...
         */
        virtual void readModel(const std::string& filename) = 0;

        /**
        * @brief 生成STL（离散背景网格）
        *
        * @param[out] tris 生成的STL三角形拓扑
        * @param[out] points 生成的STL点坐标
        */
        virtual void generateStl(
            std::vector<std::array<int, 3>>& tris,
            std::vector<std::array<double, 3>>& points) = 0;

    };

} // namespace Tiger
