#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do
                           // this in one cpp file
#include "catch.hpp"
#include "../../include/geometry_interface_occ.h"
#include "meshIO.h"
#include "../../include/alias.h"
#include"../../include/discretize_surface.h"
#include"list_to_matrix.hpp"
using namespace std;




TEST_CASE("tiger 1.5 base", "test_tran2fli") {
  /*TiGER::Geometry_OCC geo;
  geo.readModel("../cylinder.stp");[�ļ���ImportStep.cpp]
  std::vector<std::array<int, 3>> tris;
  std::vector<std::array<double, 3>> points;
  geo.generateStl(tris, points);
  cout<<geo.curves_.size()<<endl;
  cout<<geo.surfaces_.size()<<endl;*/

  TiGER::SurfaceMesh surfaceOut;
  TiGER::Geometry Geo;
  TiGER::readModelParameters arg;
  //arg.setdivideclosed_flag(1);
  //arg.setgetname_flag(0);
  //arg.setheal(7); 
  //arg.setbool_flag(1);
  //arg.setbool_para(0);
  //TiGER::Read_Model("../vki.stp", Geo, arg);
  TiGER::Read_Model("../F22.stp", Geo, arg);
  //TiGER::Read_Model("../cylinder.stp", Geo,arg);
  cout << "readModel success" << endl;
  cout << Geo.curves_.size() << endl;
  std::array<double, 3> point;
  double u;
  Geo.curves_[716]->arcToCoordinateFunction(0,point,u);
  cout<<u<<endl;



  //TiGER::dealcomplexcurve(Geo);
  //cout << Geo.curves_.size() << endl;
  ////use guide
  //std::shared_ptr<TiGER::GeometryCurve> c = Geo.curves_[0];
  //double n1=0.2, n2=0.8;
  //double r1, r2;
  //TiGER::mapNormal2Real_Curve(c, n1, r1);
  //TiGER::mapNormal2Real_Curve(c, n2, r2);
  //std::shared_ptr<TiGER::GeometryCurve> c2= TiGER::trimmedCurve(c, r1, r2);


  TiGER::generateStlParameters arg_stl;
  arg_stl.setDeflection(0.001);
  arg_stl.setAngle(5);
  arg_stl.setRelative(1);
  arg_stl.setInParallel(0);
  arg_stl.setMinSize(1e-12);
  arg_stl.setInternalVerticesMode(1);
  arg_stl.setControlSurfaceDeflection(0);
  TiGER::CAD_Tessellation(Geo, surfaceOut, arg_stl);


  //TiGER::generateStlDTParameters arg_stldt;
  //arg_stldt.setDeflection(0.1);
  //arg_stldt.setAngle(5);
  //arg_stldt.setLayer(10);
  //arg_stldt.setdemension(20);
  //arg_stldt.setdelta(0);
  //arg_stldt.setsizealpha(0.7);
  //TiGER::VCAD_Tessellation(Geo, surfaceOut, arg_stldt);

  if( surfaceOut.coords.size() == 0 )
      return;
  TiGER::Mesh m;
  TiGER::list_to_matrix<3>(surfaceOut.coords, m.Vertex);
  TiGER::list_to_matrix<3>(surfaceOut.tris, m.Topo);
  TiGER::list_to_matrix<1>(surfaceOut.attribute_int, m.Masks);
  TiGER::MESHIO::writeVTK("../111.vtk", m);
}