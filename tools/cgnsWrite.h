#pragma once
#ifndef CGNSWRITE_H
#define CGNSWRITE_H

#include "alias.h"
#include "mesh_data_struct.h"
#include <string>

using namespace TiGER;

int Tool_writeMeshToCGNS(const char *filename, SurfaceMesh &sm, VolumeMesh &vm);
int Tool_writeMeshToCGNS(const char* filename, SurfaceMesh& sm);

#include "cgnsWrite.cpp"
#endif // !CGNEWRITE_H
