#include "./api.h"
namespace TiGER{

    struct Curve;

    TiGER set_d-1_function(
        Curve &cad,
        std::function(std::array<double, 2>(const double &)) func);
    void set_d0_function(
        Curve &cad,
        std::function(std::array<double, 2>(const double &)) func);

}