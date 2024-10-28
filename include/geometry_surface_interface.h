#pragma once

#include <iostream>
#include <set>
#include <map>
#include <array>


/**
 * @brief 几何曲面基类
 *
 */
	class GeometrySurface
	{
	public:
		virtual std::array<double, 4> get_uv_scale()=0; // return umin, umax, vmin, vmax;

		virtual std::array<double, 6> get_xyz_scale() = 0; // return xmin, ymin, zmin, xmax, ymax, zmax

		/**
		 * @brief given a coord from uv space, return the cooresponding coord in 3d space
		 *
		 * @param [in] uv {u,v}
		 * @return std::array<double, 3> {x,y,z}
		 */
		virtual void d0(std::array<double, 2> uv, std::array<double, 3>& xyz) = 0;

		/**
		 * @brief
		 *
		 * @param[in] uv uv参数平面坐标
		 * @param[out] du u方向导数
		 * @param[out] dv V方向导数
		 */
		virtual void d1(const std::array<double, 2>& uv, std::array<double, 3>& du, std::array<double, 3>& dv) = 0;

		/**
		 * @brief
		 *
		 * @param[in] uv uv参数平面坐标
		 * @param[out] duu u方向二阶导数
		 * @param[out] dvv v方向二阶导数
		 * @param[out] duv uv方向二阶导数
		 */
		virtual void d2(const std::array<double, 2>& uv,
			std::array<double, 3>& duu,
			std::array<double, 3>& dvv,
			std::array<double, 3>& duv) = 0;
		/**
		 * @brief
		 *
		 * @param xyz 三维空间坐标
		 * @param eps 投影误差阈值
		 * @param hint 提示值
		 * @param res_uv 投影结果(u, v)坐标
		 */
		virtual void project(const std::array<double, 3>& xyz,
			const double& eps,
			const std::array<double, 2>& hint, std::array<double, 2>& res_uv) = 0;
	};
