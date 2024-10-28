#include <functional>
#include <iostream>

#include "./cad.h"
#include "../geometry_interface.h"
#include "../geometry_curve_interface.h"
namespace TiGER
{
    namespace API
    {
        class APIGeometryCurve : public GeometryCurve
        {
            virtual std::array<double, 3> d0(const double &u)
            {
                return d0_func(u);
            }
            virtual void d1(const double &u, std::array<double, 3> &du)
            {
                du = d1_func(u);
            }
            void set_d0_func(std::function<std::array<double, 3>(const double &)> &func)
            {
                d0_func = func;
            }
            void set_d1_func(std::function<std::array<double, 3>(const double &)> &func)
            {
                d1_func = func;
            }

        private:
            std::function<std::array<double, 3>(const double &)> d0_func;
            std::function<std::array<double, 3>(const double &)> d1_func;

        }

        struct Curve
        {
            APIGeometryCurve geo_curves;

            /* data */
        };
        TiGER set_d0_function(
            Curve &curve,
            std::function(std::array<double, 3>(const double &)) func)
        {
            curve.geo_curves.set_d0_func(func);
        }
        TiGER set_d1_function(
            Curve &curve,
            std::function(std::array<double, 3>(const double &)) func)
        {
            curve.geo_curves.set_d1_func(func);
        }
        struct Cad
        {
        };
    }

}