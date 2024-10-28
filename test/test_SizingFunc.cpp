#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "matrix_to_list.hpp"
#include "meshIO.h"
#include "sizing_field_interface.h"
#include <iostream>
#include <array>
#include "alias.h"

TEST_CASE("tiger 1.5 base", "test_sizingfunction"){
    TiGER::Mesh m;
    TiGER::MESHIO::readVTK("E://model//J20_stl.vtk", m);
    // TiGER::MESHIO::readVTK("/home/zxy/test/hummer_debug_autogrid.vtk", m);
    
    TiGER::SurfaceMesh surf;
    TiGER::BackGroundParameters args;
    TiGER::SizingManager sf;

    TiGER::matrix_to_list<3>(m.Vertex, surf.coords);
    TiGER::matrix_to_list<3>(m.Topo, surf.tris);
    TiGER::matrix_to_list<1>(m.Masks, surf.regions);
    TiGER::size_field::setSizeField(surf, args, sf);

    std::cout << "Num of source:" << sf.sf.size() << std::endl;

    TiGER::PointSourceParameters args1;
    TiGER::CylinderSourceParameters args2;
    TiGER::CubicSourceParameters arg3;

    
    // std::array<double,3> xyz0={0,0,0};
    // TiGER::size_field::addPointSource(xyz0, args1, sf);

    // std::array<std::array<double, 3>, 2> xyz1 = { 
    //     std::array<double, 3>{0.0, 0.0, 0.0}, 
    //     std::array<double, 3>{0.0, 0.0, 5.0} 
    // };
    // std::array<double,1> r={1};
    // TiGER::size_field::addCylinderSource(xyz1, r, args2, sf);

    std::array<std::array<double, 3>, 2> xyz2 = {0,0,0,0,0,10};
    std::array<double, 2> r2 = {1, 1};
    TiGER::size_field::addCylinderSource(xyz2,r2 , args2, sf);
    std::array<double, 3> q{ 0,0,0 };
    for (int i = 0; i < 10; i++) {
        q[0] += i * 0.2;
        double  ans = sf.sf[0]->getsize(q[0], q[1], q[2]);
    }
    std::cout << "success" << std::endl;

    // std::array<std::array<double, 3>, 8> xyz2 = { 
    //     std::array<double, 3>{0.0, 0.0, 0.0}, 
    //     std::array<double, 3>{1.0, 0.0, 0.0},
    //     std::array<double, 3>{0.0, 1.0, 0.0},
    //     std::array<double, 3>{1.0, 1.0, 0.0},
    //     std::array<double, 3>{0.0, 0.0, 2.0}, 
    //     std::array<double, 3>{1.0, 0.0, 2.0},
    //     std::array<double, 3>{0.0, 1.0, 2.0},
    //     std::array<double, 3>{1.0, 1.0, 2.0}
    // };
    // TiGER::size_field::addCubicSource(xyz2, arg3, sf);

    std::cout << "Num of source:" << sf.sf.size() << std::endl;
    
    //std::cout<< "Start space:" << static_cast<TiGER::cylinderSource*>(sf.sf[0].get())->args.getbeginSpacing() << std::endl;
    std::cout<< "size1:"<< static_cast<TiGER::cylinderSource*>(sf.sf[0].get())->getsize(0, 0, 3) << std::endl;
    static_cast<TiGER::cylinderSource*>(sf.sf[0].get())->args.setsourcetype(1);
    std::cout<< "size2:"<< static_cast<TiGER::cylinderSource*>(sf.sf[0].get())->getsize(0, 0, 3) << std::endl;

    

    // static_cast<TiGER::backgroundsizefield*>(sf.sf[0].get())->exportSizeFieldToVTK("E://model//J20_size_func.vtk");
}


// TEST_CASE("tiger 1.5 base", "test_sizingfunction"){
//     std::cout << "Start load" << std::endl;
//     TiGER::Mesh m;
//     TiGER::MESHIO::readVTK("E:/model/yuanpan_debug_autogrid.vtk", m, "surface_id");

//     std::cout << "Load success" << std::endl;
    
//     TiGER::SurfaceMesh surf;
//     TiGER::BackGroundParameters args;
//     TiGER::SizingManager sf;

//     std::cout << "Num of source:" << sf.sf.size() << std::endl;

//     TiGER::PointSourceParameters args1;
//     TiGER::CylinderSourceParameters args2;
//     TiGER::CubicSourceParameters arg3;

//     std::cout << "sf'size:" << sf.sf.size() << std::endl;
//     std::cout << "coords'size:" << surf.coords.size() << std::endl;
//     std::cout << "tris'size:" << surf.tris.size() << std::endl;
//     std::cout << "regions'size:" << surf.regions.size() << std::endl;
//     std::cout << "regions's id:" << surf.regions[0] << surf.regions[1] << std::endl;
    
//     TiGER::size_field::setSizeField(surf, args, sf);
//     std::cout << "sf'size:" << sf.sf.size() << std::endl;

//     // std::array<double,3> xyz0={0,0,0};
//     // TiGER::size_field::addPointSource(xyz0, args1, sf);

//     // std::array<std::array<double, 3>, 2> xyz1 = { 
//     //     std::array<double, 3>{0.0, 0.0, 0.0}, 
//     //     std::array<double, 3>{0.0, 0.0, 5.0} 
//     // };
//     // std::array<double,1> r={1};
//     // TiGER::size_field::addCylinderSource(xyz1, r, args2, sf);

//     std::array<std::array<double, 3>, 8> xyz2 = { 
//         std::array<double, 3>{0.0, 0.0, 0.0}, 
//         std::array<double, 3>{1.0, 0.0, 0.0},
//         std::array<double, 3>{0.0, 1.0, 0.0},
//         std::array<double, 3>{1.0, 1.0, 0.0},
//         std::array<double, 3>{0.0, 0.0, 2.0}, 
//         std::array<double, 3>{1.0, 0.0, 2.0},
//         std::array<double, 3>{0.0, 1.0, 2.0},
//         std::array<double, 3>{1.0, 1.0, 2.0}
//     };
//     TiGER::size_field::addCubicSource(xyz2, arg3, sf);

//     std::cout << "Num of source:" << sf.sf.size() << std::endl;
    
//     // std::cout<< "Start space:" << static_cast<TiGER::cylinderSource*>(sf.sf[1].get())->args.getbeginSpacing() << std::endl;
//     std::cout<< "formal size:"<< static_cast<TiGER::cylinderSource*>(sf.sf[0].get())->getsize(0.5, 0.5, 1) << std::endl;
//     std::cout<< "formal size:"<< static_cast<TiGER::cylinderSource*>(sf.sf[0].get())->getsize(0.5, 0.5, 3) << std::endl;

//     static_cast<TiGER::cylinderSource*>(sf.sf[0].get())->args.setsourcetype(1);
//     std::cout<< "constant size:"<< static_cast<TiGER::cylinderSource*>(sf.sf[0].get())->getsize(0.5, 0.5, 1) << std::endl;
//     std::cout<< "constant size:"<< static_cast<TiGER::cylinderSource*>(sf.sf[0].get())->getsize(0.5, 0.5, 3) << std::endl;

// }


// TEST_CASE("tiger 1.5 base", "test_sizingfunction"){
//     std::cout << "Start load" << std::endl;
//     TiGER::Mesh m;
//     TiGER::MESHIO::readVTK("E:/model/solidassem3_debug_autogrid.vtk", m, "surface_id");

//     std::cout << "Load success" << std::endl;
    
//     TiGER::SurfaceMesh surf;
//     TiGER::BackGroundParameters args;
//     TiGER::SizingManager sf;

//     std::cout << "Num of source:" << sf.sf.size() << std::endl;

//     TiGER::PointSourceParameters args1;
//     TiGER::CylinderSourceParameters args2;
//     TiGER::CubicSourceParameters arg3;

//     TiGER::matrix_to_list<3>(m.Vertex, surf.coords);
//     TiGER::matrix_to_list<3>(m.Topo, surf.tris);
//     TiGER::matrix_to_list<1>(m.Masks, surf.regions);

//     std::cout << "sf'size:" << sf.sf.size() << std::endl;
//     std::cout << "coords'size:" << surf.coords.size() << std::endl;
//     std::cout << "tris'size:" << surf.tris.size() << std::endl;
//     std::cout << "regions'size:" << surf.regions.size() << std::endl;

//     //std::cout << "regions's id:" << surf.regions[0] << surf.regions[1] << std::endl;
    
//     TiGER::size_field::setSizeField(surf, args, sf);
//     std::cout << "sf'size:" << sf.sf.size() << std::endl;

//     std::array<double,3> xyz0={0,0,0};
//     std::array<double,3> xyz1={1,0,0};
//     std::array<double,3> xyz2={3,4,0};
//     std::array<double,3> xyz3={1,2,3};
    
//     std::cout << "Generate success" << std::endl;
//     std::cout << "Size:" << TiGER::size_field::getSizeValue(xyz0[0],xyz0[1],xyz0[2],sf) << std::endl;
//     REQUIRE(TiGER::size_field::getSizeValue(xyz0[0],xyz0[1],xyz0[2],sf)>0);


//     std::vector<int> point_index = {4681};
//     std::vector<int> line_index = {8, 14};
//     std::vector<int> face_index = {7};
    
//     std::cout << static_cast<TiGER::backgroundsizefield*>(sf.sf[0].get())->backgroundmesh.point_attribute_double[4681] <<std::endl;
//     std::cout << static_cast<TiGER::backgroundsizefield*>(sf.sf[0].get())->backgroundmesh.point_attribute_double[8908] <<std::endl;
//     std::cout << static_cast<TiGER::backgroundsizefield*>(sf.sf[0].get())->backgroundmesh.point_attribute_double[3837] <<std::endl;
//     std::cout << static_cast<TiGER::backgroundsizefield*>(sf.sf[0].get())->backgroundmesh.point_attribute_double[4274] <<std::endl;
//     std::cout << static_cast<TiGER::backgroundsizefield*>(sf.sf[0].get())->backgroundmesh.point_attribute_double[4271] <<std::endl;


//     std::cout<<"Begin modify"<<std::endl;
//     TiGER::size_field::setAttributeSize("point", point_index, 1.0, *static_cast<TiGER::backgroundsizefield*>(sf.sf[0].get()));
//     std::cout<<"Modify 1"<<std::endl;
//     std::cout << static_cast<TiGER::backgroundsizefield*>(sf.sf[0].get())->backgroundmesh.point_attribute_double[4681] <<std::endl;
//     TiGER::size_field::setAttributeSize("line", line_index, 1.0, *static_cast<TiGER::backgroundsizefield*>(sf.sf[0].get()));
//     std::cout<<"Modify 2"<<std::endl;
//     std::cout << static_cast<TiGER::backgroundsizefield*>(sf.sf[0].get())->backgroundmesh.point_attribute_double[4274] <<std::endl;
//     std::cout << static_cast<TiGER::backgroundsizefield*>(sf.sf[0].get())->backgroundmesh.point_attribute_double[4271] <<std::endl;
//     TiGER::size_field::setAttributeSize("face", face_index, 1.0, *static_cast<TiGER::backgroundsizefield*>(sf.sf[0].get()));
//     std::cout << static_cast<TiGER::backgroundsizefield*>(sf.sf[0].get())->backgroundmesh.point_attribute_double[3837] <<std::endl;
//     std::cout << static_cast<TiGER::backgroundsizefield*>(sf.sf[0].get())->backgroundmesh.point_attribute_double[3790] <<std::endl;
//     std::cout << static_cast<TiGER::backgroundsizefield*>(sf.sf[0].get())->backgroundmesh.point_attribute_double[3903] <<std::endl;
//     std::cout<<"Modify 3"<<std::endl;

// }




