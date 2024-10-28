#pragma once
#ifndef _GEOMETRY_INTERFACE_PLANE_H_
#define _GEOMETRY_INTERFACE_PLANE_H_

#include "geometry_data_struct.h"
#include "alias.h"

namespace TiGER
{
    void Grid1D_createMeshLine(std::array<double, 3> p0,
        std::array<double, 3> P1,
        TiGER::GeometryCurve& line);
    void Grid2D_createMeshPlane(std::array<double, 3> p0,
        std::array<double, 3> basis1,
        std::array<double, 3> basis2,
        GeometrySurface& plane);
}

#include "geometry_interface_plane.cpp"
#endif // !_GEOMETRY_INTERFACE_PLANE_H_