#include "../include/curve.h"
#include "../geometry_curve_interface.h"
namespace TiGER{
    struct curve
    {
        std::function(std::array<double, 3>(const double &)) d0_func;
        TiGER::GeometryCurve geo_curve;
        /* data */
    };
     TiGER set_d0_function(
        Curve &cad,
        std::function(std::array<double, 3>(const double &)) func);
}