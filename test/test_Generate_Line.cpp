#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "geometry_data_struct.h"
#include "discretize_curve.h"
#include "geometry_interface_plane.h"
#include "meshIO.h"
#include <cstdio>
#include "alias.h"

TEST_CASE("tiger 1.5 base", "test_extrude")
{
    std::array<double, 3> p0 = {0, 0, 0};
    std::array<double, 3> p1 = {100, 100, 100};
    TiGER::GeometryCurve line;
    TiGER::Grid1D_createMeshLine(p0, p1, line);
    std::shared_ptr<TiGER::GeometryCurve> line_ptr = std::make_shared<TiGER::GeometryCurve>(std::move(line));
    std::vector<std::shared_ptr<GeometryCurve>> lines;
    std::vector<TiGER::LineMesh> segements;
    lines.push_back(line_ptr);
    TiGER::CurveParametersTanh args;
    args.setdimension(25);
    args.setbegin_s(3);
    args.setend_s(10);
    TiGER::discretize_curve::Grid1D_Tanh(lines, args, segements);
    std::string filename = "D:/mesh/TiGER/test/linemesh.vtk";
    MESHIO::writeVTK(filename, segements[0], " ");
}