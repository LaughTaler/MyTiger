#include "cgnsWrite.h"
#include "matrix_to_list.hpp"
#include "meshIO.h"
#include <cstdio>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include "alias.h"

#define BUFFER_LENGTH 256
using namespace MESHIO;

int main(int argc, char* argv[])
{

    const char* input_v_file = argv[1];
    const char* output_file = argv[2];

	TiGER::SurfaceMesh sm;
	TiGER::VolumeMesh vm;
    Mesh sfmesh;
    Mesh vlmesh;

   // MESHIO::readVTK(input_m_file, sfmesh, "surface_id");
    MESHIO::readVTK(input_v_file, vlmesh, "surface_id");
 //   matrix_to_list<3>(sfmesh.Vertex, sm.coords);
//    matrix_to_list<3>(sfmesh.Topo, sm.tris);
  //  matrix_to_list<1>(sfmesh.Masks, sm.attribute_int);
    matrix_to_list<3>(vlmesh.Vertex, vm.coords);
    matrix_to_list<6>(vlmesh.Topo, vm.prisms);
    for (int i = 0; i < 10; i++)
        sm.bc[i]= TiGER::BoundaryCondition::WALL;
	const char *filename;

    filename = output_file;
	Tool_writeMeshToCGNS(filename, sm, vm);
	return 0;
}