#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do
                           // this in one cpp file
#include "catch.hpp"
#include "geometry_interface_occ.h"
#include "alias.h"
#include"./discretize_surface.h"
using namespace std;

TEST_CASE("tiger 1.5 base", "test_tran2fli") {
  /*TiGER::Geometry_OCC geo;
  geo.readModel("../cylinder.stp");
  std::vector<std::array<int, 3>> tris;
  std::vector<std::array<double, 3>> points;
  geo.generateStl(tris, points);
  cout<<geo.curves_.size()<<endl;
  cout<<geo.surfaces_.size()<<endl;*/

  std::vector<std::array<int, 3>> tris;
  std::vector<std::array<double, 3>> points;
  TiGER::Geometry Geo;
  TiGER::Read_Model("../cylinder.stp", Geo);
  TiGER::CAD_Tessellation(Geo, tris, points);




}