#pragma once
#include "geometry_curve_interface.h"
#include "geometry_surface_interface.h"
#include <vector>
#include <memory>
namespace SmeshGen

{
    class Geometry {
    public:
        /** ������Ϣ */
        std::vector<std::shared_ptr<GeometryCurve>> curves_;
        std::vector<std::shared_ptr<GeometrySurface>> surfaces_;

        /** ������Ϣ */
        std::vector<std::vector<int>> topo_surf_to_curves_; // �浽���ߵ�������Ϣ
        std::vector<std::vector<int>> topo_body_to_surfs_;  // �嵽���������Ϣ

        /**
         * @brief ��ȡ���ε��ڴ���
         *
         * @param[in] filename �����ļ�����֧�� stp, igs, step, iges, sat, brep, x_t ...
         */
        virtual void readModel(const std::string& filename) = 0;

        /**
        * @brief ����STL����ɢ��������
        *
        * @param[out] tris ���ɵ�STL����������
        * @param[out] points ���ɵ�STL������
        */
        virtual void generateStl(
            std::vector<std::array<int, 3>>& tris,
            std::vector<std::array<double, 3>>& points) = 0;

    };

} // namespace Tiger
