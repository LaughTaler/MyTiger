#pragma once
#include <iostream>
#include <set>
#include <map>
#include <array>

	/**
	 * @brief 读取连续几何生成离散几何（离散三角形网格）
	 *
	 */
	class GeometryCurve
	{

	public:
		/**
		 * @brief given a coord from uv space, return the cooresponding coord in 3d space
		 *
		 * @param [in] uv {u,v}
		 * @return std::array<double, 3> {x,y,z}
		 */
		virtual std::array<double, 3> d0(double& u) = 0;

		/**
		 * @brief
		 *
		 * @param[in] uv uv参数平面坐标
		 * @param[out] du u方向导数
		 */
		virtual void d1(const double& u, std::array<double, 3>& du) = 0;

		/**
		 * @brief
		 *
		 * @param[in] u u参数坐标
		 * @param[out] ddu 二阶导数
		 */
		virtual void d2(const double& u, std::array<double, 3>& d1u, std::array<double, 3>& d2u) = 0;

		virtual double get_length() = 0;


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
			const double& hint,
			double& res_u) = 0;

		virtual void arcToCoordinate(double arc_length, std::array<double, 3>& xyz,double &res_u) = 0;
	};