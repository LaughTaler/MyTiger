#include <functional>

#include "./api.h"
namespace TiGER
{
    namespace API
    {
        struct Cad;
        struct Curve;
        struct Surface;

        
        TiGER::STATE set_d1_function(
                          Curve &curve,
                          std::function(std::array<double, 2>(const double &)) func);
        TiGER::STATE set_d0_function(
            Curve &curve,
            std::function(std::array<double, 2>(const double &)) func);
    }

}
