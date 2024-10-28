#pragma once

#pragma once
#ifndef DISCRETIZE_SURFACE_PARAMETERS_H_
#define DISCRETIZE_SURFACE_PARAMETERS_H_
#include "./parameters.h"
namespace TiGER
{
struct IsotropicSurfaceParametersTri : public APIPara
{
public:
    IsotropicSurfaceParametersTri()
        : mUserSpecifiedMaxEdgeLen(0)
        , mUserSpecifiedMinEdgeLen(0)
        , mMaxAngle(0)
        , mMaxDeviation(0)
        , mSwapCellWithNoInteriorPoints(false)
        , OutputInformation(false)
        , Savepath("")
        , ExcessRate(1.25)
        , iteration(5)
    {
    }

    DEFINE_SETTER(double, mUserSpecifiedMaxEdgeLen, IsotropicSurfaceParametersTri, 0);
    DEFINE_SETTER(double, mUserSpecifiedMinEdgeLen, IsotropicSurfaceParametersTri, 0);
    DEFINE_SETTER(double, mMaxAngle, IsotropicSurfaceParametersTri, 0);
    DEFINE_SETTER(double, mMaxDeviation, IsotropicSurfaceParametersTri, 0);
    DEFINE_SETTER(bool, mSwapCellWithNoInteriorPoints, IsotropicSurfaceParametersTri, false);
    DEFINE_SETTER(bool, OutputInformation, IsotropicSurfaceParametersTri, false);
    DEFINE_SETTER(std::string, Savepath, IsotropicSurfaceParametersTri, "");
    DEFINE_SETTER(double, ExcessRate, IsotropicSurfaceParametersTri, 1.25);
    DEFINE_SETTER(int, iteration, IsotropicSurfaceParametersTri, 5);
    
};

struct AnisotropicSurfaceParametersTri : public APIPara
{
public:
    AnisotropicSurfaceParametersTri(){};
    DEFINE_SETTER(int, MaxLayer, AnisotropicSurfaceParametersTri, 0);
    DEFINE_SETTER(double, growth_ratio, AnisotropicSurfaceParametersTri, 1.2);
    DEFINE_SETTER(bool, Topology_priority, AnisotropicSurfaceParametersTri, true);
    DEFINE_SETTER(double, Angle_threshold, AnisotropicSurfaceParametersTri, 5.0);
};
} // namespace TiGER

#endif