#pragma once

#include <iostream>
#include <set>
#include <map>
#include <array>


/**
 * @brief �����������
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
		 * @param[in] uv uv����ƽ������
		 * @param[out] du u������
		 * @param[out] dv V������
		 */
		virtual void d1(const std::array<double, 2>& uv, std::array<double, 3>& du, std::array<double, 3>& dv) = 0;

		/**
		 * @brief
		 *
		 * @param[in] uv uv����ƽ������
		 * @param[out] duu u������׵���
		 * @param[out] dvv v������׵���
		 * @param[out] duv uv������׵���
		 */
		virtual void d2(const std::array<double, 2>& uv,
			std::array<double, 3>& duu,
			std::array<double, 3>& dvv,
			std::array<double, 3>& duv) = 0;
		/**
		 * @brief
		 *
		 * @param xyz ��ά�ռ�����
		 * @param eps ͶӰ�����ֵ
		 * @param hint ��ʾֵ
		 * @param res_uv ͶӰ���(u, v)����
		 */
		virtual void project(const std::array<double, 3>& xyz,
			const double& eps,
			const std::array<double, 2>& hint, std::array<double, 2>& res_uv) = 0;
	};
