#ifndef PYRAMID_SURFACE
#define PYRAMID_SURFACE
#include "./mesh_data_struct.h"
#include "./Eigen/Dense"

namespace TiGER
{
void pyramid_surface(const SurfaceMesh& surfacemeshin, VolumeMesh& volumemesh, SurfaceMesh& surfacemeshout);
}

#include "Pyramid_surface.cpp"
#endif //!_PYRAMID_SURFACE
