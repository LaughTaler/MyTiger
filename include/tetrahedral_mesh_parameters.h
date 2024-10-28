#pragma once
#ifndef _TETRAHEDRAL_PARAMETERS_H_
#define _TETRAHEDRAL_PARAMETERS_H_
#include "./parameters.h"
namespace TiGER
{
class TetrahedraParameters : public APIPara
{
    DEFINE_SETTER(int, infolevel, TetrahedraParameters, 1);
    DEFINE_SETTER(int, debug, TetrahedraParameters, 0);
    DEFINE_SETTER(int, constrain, TetrahedraParameters, 1);
    DEFINE_SETTER(int, refine, TetrahedraParameters, 1);
    DEFINE_SETTER(int, optlevel, TetrahedraParameters, 3);
    DEFINE_SETTER(int, optloop, TetrahedraParameters, 3);
    DEFINE_SETTER(int, nthread, TetrahedraParameters, -1);
    DEFINE_SETTER(double, size, TetrahedraParameters, 1);
    DEFINE_SETTER(double, decay, TetrahedraParameters, -1);
    DEFINE_SETTER(double, minEdge, TetrahedraParameters, -1);
    DEFINE_SETTER(double, maxEdge, TetrahedraParameters, -1);
    DEFINE_SETTER(double, minSize, TetrahedraParameters, -1);
    DEFINE_SETTER(double, maxSize, TetrahedraParameters, -1);
    DEFINE_SETTER(double, backgroundSpace, TetrahedraParameters, -1);
    DEFINE_SETTER(double, growsize, TetrahedraParameters, -1);
    DEFINE_SETTER(double, optangle, TetrahedraParameters, 15);
    DEFINE_SETTER(double, optratio, TetrahedraParameters, 0.2);
    DEFINE_SETTER(std::vector<int>, layer, TetrahedraParameters, {});
    DEFINE_SETTER(std::vector<double>, hole, TetrahedraParameters, {});
};
} // namespace TiGER

#endif // _EXTRUDE_PARAMETERS_H_