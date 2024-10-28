#pragma once
#include "mesh_data_struct.h"
namespace TiGER
{
using ConstrainDelaunayFunc = void (*)(const SurfaceMesh&, VolumeMesh&);
} // namespace TiGER
