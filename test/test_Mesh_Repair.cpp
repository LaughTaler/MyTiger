#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include <chrono>
// #include"geometry_data_struct.h"
// #include"geometry_interface_occ.h"

// #include "../tools/mesh_repair.h"
#include "alias.h"
#include "meshIO.h"
#include "../tools/matrix_to_list.hpp"
#include "../tools/list_to_matrix.hpp"
#include "../tools/mesh_repair.h"
#include "../tools/merge_volume_mesh.hpp"
#include "../tools/merge_surface_mesh.hpp"
#include "discrete_geometry_repair.h"
#include <queue>
// using namespace std;


// TEST_CASE("tiger 1.5 Tool_nativelyDetectPointDuplicates", "test_nativedetectDuplicatePoints") {
//     TiGER::Mesh mesh;    
//     TiGER::MESHIO::readVTK("C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/cylinder/test_AFT_Tri_output.vtk",mesh);
//     std::vector<std::array<double, 3>> points_in;
//     TiGER::matrix_to_list<3>(mesh.Vertex, points_in);
//     double eps = 1e-4;
//     std::vector<int> bcj_pre;
//     int num = TiGER::repair::Tool_nativelyDetectPointDuplicates(points_in,eps,bcj_pre);
//     cout<<"duplicate num : "<<num<<endl;
    
//     std::string filename = "C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/cylinder/1.txt";
//     FILE *f = fopen(filename.c_str(), "w");

//     for(int i=0;i<bcj_pre.size();i++){
//         fprintf(f, "%d ", bcj_pre[i]);
//     }
// }

// TEST_CASE("tiger 1.5 Tool_detectPointDuplicates", "test_detectDuplicatePoints") {
//     // std::vector<std::array<double,3>> points_in(3);
//     // points_in[0] = {-1,0,0};
//     // points_in[1] = {1,0,0};
//     // points_in[2] = {0,0,0};
//     // double eps = 1;
//     // std::vector<int> bcj_pre;
//     // int num = TiGER::repair::Tool_detectPointDuplicates(points_in,eps,bcj_pre);
//     // std::cout<<"duplicate num : "<<num<<std::endl;

//     TiGER::Mesh mesh;    
//     TiGER::MESHIO::readVTK("C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/dirty_geometry/aircraft-14_stl.vtk",mesh);
//     std::vector<std::array<double, 3>> points_in;
//     std::vector<std::array<int, 3>> tris;
//     TiGER::matrix_to_list<3>(mesh.Vertex, points_in);
//     TiGER::matrix_to_list<3>(mesh.Topo,tris);
//     double eps = 1e-5;
//     std::vector<int> bcj_pre;

//     auto start = std::chrono::high_resolution_clock::now();
//     int num = TiGER::repair::Tool_detectPointDuplicates(points_in,eps,bcj_pre);
//     auto finish = std::chrono::high_resolution_clock::now();
//     // 计算持续时间
//     std::chrono::duration<double> elapsed = finish - start;
 
//     // 输出运行时间
//     std::cout << "test_detectPunctureSurface time: " << elapsed.count() << " sec" << std::endl;

//     cout<<"duplicate num : "<<num<<endl;
    
//     std::vector<std::array<double, 3>> points_out;
//     std::unordered_map<int, int> pt_index;
//     int index = 0;
//     for (int i = 0; i < bcj_pre.size(); i++)  // 删除多余 点
//     {
//         if (bcj_pre[i] == -1)
//         {
//             points_out.push_back(points_in[i]);
//             pt_index[i] = index;
//             index++;
//         }
//     }

//     for (int i = 0; i < tris.size(); i++)
//     {
//         for (int j = 0; j < 3; j++)
//         {
//             if (bcj_pre[tris[i][j]] != -1)
//                 tris[i][j] = pt_index[bcj_pre[tris[i][j]]];
//             else
//                 tris[i][j] = pt_index[tris[i][j]];
//         }
//     }

//     TiGER::Mesh mesh_tri;
//     TiGER::list_to_matrix<3>(points_out,mesh_tri.Vertex);
//     TiGER::list_to_matrix<3>(tris,mesh_tri.Topo);

//     TiGER::MESHIO::writeVTK("C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/myexample/UAV-0309_surface_mesh.o.vtk",mesh_tri);
// }

// TEST_CASE("tiger 1.5 mergemesh", "test_mergemesh") {
//     vector<TiGER::SurfaceMesh> suf_in;
//     TiGER::SurfaceMesh sur_out;

//     for(int i=0;i<65;i++){
//         TiGER::SurfaceMesh piece;
//         TiGER::Mesh mesh;  
//         std::vector<std::array<double, 3>> points_in;
//         std::vector<std::array<int, 3>> tris;  
//         std::string filename = "face" + std::to_string(i) + ".vtk";
//         TiGER::MESHIO::readVTK("C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/face/"+filename,mesh);
//         TiGER::matrix_to_list<3>(mesh.Vertex, points_in);
//         TiGER::matrix_to_list<3>(mesh.Topo,tris);
//         std::cout<<"mesh"<<i<<" point size "<<points_in.size()<<std::endl;
//         piece.coords = points_in;
//         piece.tris = tris;
//         suf_in.push_back(piece);
//     }

//     // TiGER::SurfaceMesh piece;
//     // TiGER::Mesh mesh;  
//     // std::vector<std::array<double, 3>> points_in;
//     // std::vector<std::array<int, 3>> tris;  
//     // TiGER::MESHIO::readVTK("C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/merge/mergemesh4.vtk",mesh);
//     // TiGER::matrix_to_list<3>(mesh.Vertex, points_in);
//     // TiGER::matrix_to_list<3>(mesh.Topo,tris);
//     // std::cout<<"mesh1 point size "<<points_in.size()<<std::endl;
//     // piece.coords = points_in;
//     // piece.tris = tris;
//     // suf_in.push_back(piece);

//     // TiGER::MESHIO::readVTK("C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/merge/mergemesh1.vtk",mesh);
//     // TiGER::matrix_to_list<3>(mesh.Vertex, points_in);
//     // TiGER::matrix_to_list<3>(mesh.Topo,tris);
//     // std::cout<<"mesh2 point size "<<points_in.size()<<std::endl;
//     // piece.coords = points_in;
//     // piece.tris = tris;
//     // suf_in.push_back(piece);

//     // TiGER::MESHIO::readVTK("C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/merge/mergemesh0.vtk",mesh);
//     // TiGER::matrix_to_list<3>(mesh.Vertex, points_in);
//     // TiGER::matrix_to_list<3>(mesh.Topo,tris);
//     // std::cout<<"mesh3 point size "<<points_in.size()<<std::endl;
//     // piece.coords = points_in;
//     // piece.tris = tris;
//     // suf_in.push_back(piece);

//     // TiGER::MESHIO::readVTK("C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/merge/mergemesh2.vtk",mesh);
//     // TiGER::matrix_to_list<3>(mesh.Vertex, points_in);
//     // TiGER::matrix_to_list<3>(mesh.Topo,tris);
//     // std::cout<<"mesh4 point size "<<points_in.size()<<std::endl;
//     // piece.coords = points_in;
//     // piece.tris = tris;
//     // suf_in.push_back(piece);

//     // TiGER::MESHIO::readVTK("C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/merge/mergemesh3.vtk",mesh);
//     // TiGER::matrix_to_list<3>(mesh.Vertex, points_in);
//     // TiGER::matrix_to_list<3>(mesh.Topo,tris);
//     // std::cout<<"mesh5 point size "<<points_in.size()<<std::endl;
//     // piece.coords = points_in;
//     // piece.tris = tris;
//     // suf_in.push_back(piece);
    
//     auto start = std::chrono::high_resolution_clock::now();
//     TiGER::merge_surface_mesh(suf_in,sur_out);
//     auto finish = std::chrono::high_resolution_clock::now();
//     // 计算持续时间
//     std::chrono::duration<double> elapsed = finish - start;
 
//     // 输出运行时间
//     std::cout << "time: " << elapsed.count() << " sec" << std::endl;

//     std::cout<<"after merge point size "<<sur_out.coords.size()<<std::endl; 

//     TiGER::Mesh mesh_tri;
//     TiGER::list_to_matrix(sur_out.coords,mesh_tri.Vertex);
//     TiGER::list_to_matrix(sur_out.tris,mesh_tri.Topo);
//     TiGER::MESHIO::writeVTK("C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/face/result0.5.vtk",mesh_tri);
// } 


// TEST_CASE("tiger 1.5 mergemesh", "test_mergemesh") {
//     vector<TiGER::VolumeMesh> volume_in;
//     TiGER::VolumeMesh volume_out;

//     TiGER::VolumeMesh piece;
//     TiGER::MESHIO::readVTK_toVolumemesh("C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/merge/volumne_mesh.vtk",piece);
//     std::cout<<"mesh1 point size "<<piece.coords.size()<<std::endl;
//     volume_in.push_back(piece);

//     TiGER::VolumeMesh piece1;
//     TiGER::MESHIO::readVTK_toVolumemesh("C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/merge/vol.vtk",piece1);
//     std::cout<<"mesh2 point size "<<piece1.coords.size()<<std::endl;
//     volume_in.push_back(piece1);
    
//     auto start = std::chrono::high_resolution_clock::now();
//     TiGER::merge_volume_mesh(volume_in,volume_out);
//     auto finish = std::chrono::high_resolution_clock::now();
//     // 计算持续时间
//     std::chrono::duration<double> elapsed = finish - start;
 
//     // 输出运行时间
//     std::cout << "time: " << elapsed.count() << " sec" << std::endl;

//     std::cout<<"after merge point size "<<volume_out.coords.size()<<std::endl; 

//     TiGER::MESHIO::writeVTK("C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/merge/volume_result.vtk",volume_out);
// } 

// TEST_CASE("tiger 1.5 TMeshSurf_detectSurfacePunctures", "test_detectPunctureSurface") {
//     TiGER::Mesh mesh;    
//     TiGER::MESHIO::readVTK("C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/dirty_geometry/aircraft-14_stl.vtk",mesh);
    
//     std::vector<std::array<double, 3>> coords;
//     std::vector<std::array<int, 3>> tris;
//     TiGER::matrix_to_list<3>(mesh.Vertex, coords);
//     TiGER::matrix_to_list<3>(mesh.Topo, tris);

//     std::cout<<"point size "<<coords.size()<<std::endl;
//     std::cout<<"tri size "<<tris.size()<<std::endl; 

//     TiGER::SurfaceMesh surf_in;
//     surf_in.coords = coords;
//     surf_in.tris = tris;

//     std::vector<TiGER::discrete_geometry_repair::EdgeTriPair> edge_tri;

//     auto start = std::chrono::high_resolution_clock::now();
//     TiGER::discrete_geometry_repair::TMeshSurf_detectSurfacePunctures(
//         surf_in,
//         edge_tri
//     );
//     auto finish = std::chrono::high_resolution_clock::now();
//     // 计算持续时间
//     std::chrono::duration<double> elapsed = finish - start;
 
//     // 输出运行时间
//     std::cout << "test_detectPunctureSurface time: " << elapsed.count() << " sec" << std::endl;
//     std::cout<<"puncture pair size "<<edge_tri.size()<<std::endl;

//     TiGER::Mesh mesh_edge;  
//     mesh_edge.Vertex = mesh.Vertex;

//     std::map<std::array<int,2>,int> edge_list;
//     for(int i=0;i<edge_tri.size();i++){
//         edge_list[edge_tri[i].edge]++;        
//     }
//     std::cout<<" puncture edge size "<<edge_list.size()<<std::endl;
//     mesh_edge.Topo.resize(edge_list.size(),2);
//     int x = 0;
//     for(auto& edge_count : edge_list){
//         for(int j=0;j<2;j++){
//             mesh_edge.Topo(x,j) = edge_count.first[j];
//         }
//         x++;
//     }

//     TiGER::MESHIO::writeVTK("C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/dirty_geometry/aircraft-14_stl_punc_edge.vtk",mesh_edge);

//     TiGER::Mesh mesh_tri;
//     mesh_tri.Vertex = mesh.Vertex;

//     // std::unordered_set<int> tri_list;
//     std::map<int,int> tri_list;
//     for(int i=0;i<edge_tri.size();i++){
//         tri_list[edge_tri[i].tri]++;
//     }
//     std::cout<<" puncture tri size "<<tri_list.size()<<std::endl;
//     mesh_tri.Topo.resize(tri_list.size(),3);
//     int index = 0;
//     for(auto tri_count : tri_list){
//         for(int j=0;j<3;j++){
//             mesh_tri.Topo(index,j) = tris[tri_count.first][j];
//         }
//         index++;
//     }
//     TiGER::MESHIO::writeVTK("C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/dirty_geometry/aircraft-14_stl_punc_tri.vtk",mesh_tri);
// }

// TEST_CASE("tiger 1.5 detectSufaceQuality", "test_detectSufaceQuality"){
//     TiGER::Mesh mesh;    
//     TiGER::MESHIO::readVTK("C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/cylinder/cylinder_stl.vtk",mesh);
    
//     std::vector<std::array<double, 3>> coords;
//     std::vector<std::array<int, 3>> tris;
//     TiGER::matrix_to_list<3>(mesh.Vertex, coords);
//     TiGER::matrix_to_list<3>(mesh.Topo, tris);

//     TiGER::TiGER::DiscreteGeometry surf_in;
//     std::vector<TiGER::TiGER::DiscreteGeometry::TriangleIndex> tris_out;

//     surf_in.coords = coords;
//     std::shared_ptr<TiGER::DiscreteGeometrySurface> surface_ptr = std::make_shared<TiGER::DiscreteGeometrySurface>();
//     surface_ptr->tris = tris;
//     surf_in.surfaces.push_back(surface_ptr);

//     TiGER::RepairParameters args;
//     args.setoperator_option(0);
//     std::vector<double> tolerance(2);
//     tolerance[0] = 100;
//     args.settolerance(tolerance);

//     TiGER::discrete_geometry_repair::detectSufaceQuality(
//         surf_in,
//         args,
//         tris_out
//     );

//     std::cout<<"tri size "<<tris_out.size()<<std::endl; 
//     for(int i=0;i<tris_out.size();i++){
//         std::cout<<tris_out[i][1]<<" ";
//     }
//     std::cout<<std::endl;
// }

// TEST_CASE("tiger 1.5 TMeshSurf_identifyMeshAdjacency", "test_detectAdjacence") {
//     TiGER::Mesh mesh;    
//     TiGER::MESHIO::readVTK("C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/boeing/awacs20new_stl.vtk",mesh);
    
//     std::vector<std::array<double, 3>> coords;
//     std::vector<std::array<int, 3>> tris;
//     TiGER::matrix_to_list<3>(mesh.Vertex, coords);
//     TiGER::matrix_to_list<3>(mesh.Topo, tris);

//     TiGER::TiGER::DiscreteGeometry surf_in;
//     std::vector<TiGER::TiGER::DiscreteGeometry::Triangle_Triangle> tri_out;

//     surf_in.coords = coords;
//     std::shared_ptr<TiGER::DiscreteGeometrySurface> surface_ptr = std::make_shared<TiGER::DiscreteGeometrySurface>();
//     surface_ptr->tris = tris;
//     surf_in.surfaces.push_back(surface_ptr);

//     TiGER::RepairParameters args;
//     args.setoperator_option(1);
//     std::vector<double> tolerance(2);
//     tolerance[0] = 0.01;
//     args.settolerance(tolerance);

//     TiGER::discrete_geometry_repair::TMeshSurf_identifyMeshAdjacency(
//         surf_in,
//         args,
//         tri_out  // 两个面id 一对
//     );

//     std::cout<<"size "<<tri_out.size()<<std::endl;

//     // TiGER::Mesh mesh_tri;
//     // mesh_tri.Vertex = mesh.Vertex;
//     // mesh_tri.Topo.resize(2*tri_out.size(),3);
//     // for(int i=0;i<tri_out.size();i++){
//     //     for(int k=0;k<2;k++){
//     //         for(int j=0;j<3;j++){
//     //             mesh_tri.Topo(2*i+k,j) = tri_out[i][k][j];
//     //         }
//     //     }
//     // }
//     // TiGER::MESHIO::writeVTK("C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/boeing/tri.vtk",mesh_tri);
// }


// TEST_CASE("tiger 1.5 TMeshSurf_identifyFreeEdges", "test_detectFreeEdge") {
//     TiGER::Mesh mesh;    
//     TiGER::MESHIO::readVTK("C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/dirty_geometry/aircraft-14_stl.vtk",mesh);
    
//     std::vector<std::array<double, 3>> coords;
//     std::vector<std::array<int, 3>> tris;
//     TiGER::matrix_to_list<3>(mesh.Vertex, coords);
//     TiGER::matrix_to_list<3>(mesh.Topo, tris);

//     TiGER::SurfaceMesh surf_in;
//     surf_in.coords = coords;
//     surf_in.tris = tris;
//     surf_in.v_conn_points.resize(coords.size());
//     surf_in.v_conn_tris.resize(coords.size());

//     for(int i=0;i<tris.size();i++){
//         int p1 = tris[i][0];
//         int p2 = tris[i][1];
//         int p3 = tris[i][2];
//         surf_in.v_conn_tris[p1].insert(i);
//         surf_in.v_conn_tris[p2].insert(i);
//         surf_in.v_conn_tris[p3].insert(i);

//         surf_in.v_conn_points[p1].insert(p2);
//         surf_in.v_conn_points[p1].insert(p3);

//         surf_in.v_conn_points[p2].insert(p1);
//         surf_in.v_conn_points[p2].insert(p3);

//         surf_in.v_conn_points[p3].insert(p1);
//         surf_in.v_conn_points[p3].insert(p2);
//     }

//     std::vector<TiGER::DiscreteGeometry::Edge> edges_out;

//     TiGER::discrete_geometry_repair::Tool_initializeData(surf_in);
//     auto start = std::chrono::high_resolution_clock::now();
//     TiGER::discrete_geometry_repair::TMeshSurf_identifyFreeEdges(
//         surf_in,
//         edges_out  // 两个点id定一条边
//     );
//     auto finish = std::chrono::high_resolution_clock::now();
//     // 计算持续时间
//     std::chrono::duration<double> elapsed = finish - start;
//     // 输出运行时间
//     std::cout << "test_detectFreeEdge time: " << elapsed.count() << " sec" << std::endl;
//     std::cout<<"freeedge size "<<edges_out.size()<<std::endl;

//     TiGER::Mesh mesh_edge;  
//     mesh_edge.Vertex = mesh.Vertex;
//     mesh_edge.Topo.resize(edges_out.size(),2);
//     for(int i=0;i<edges_out.size();i++){
//         for(int j=0;j<2;j++){
//             mesh_edge.Topo(i,j) = edges_out[i][j];
//         }
//     }
//     TiGER::MESHIO::writeVTK("C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/dirty_geometry/aircraft-14_stl_freeedge.vtk",mesh_edge);
// }

// TEST_CASE("tiger 1.5 TMeshSurf_identifyMultipleEdges", "test_detectMultiedge") {

//     TiGER::Mesh mesh;    
//     TiGER::MESHIO::readVTK("C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/dirty_geometry/aircraft-14_stl.vtk",mesh);
    
//     std::vector<std::array<double, 3>> coords;
//     std::vector<std::array<int, 3>> tris;
//     TiGER::matrix_to_list<3>(mesh.Vertex, coords);
//     TiGER::matrix_to_list<3>(mesh.Topo, tris);

//     // std::cout<<"point size "<<coords.size()<<std::endl;
//     // std::cout<<"tri size "<<tris.size()<<std::endl; 

//     TiGER::SurfaceMesh surf_in;
//     surf_in.coords = coords;
//     surf_in.tris = tris;
//     surf_in.v_conn_points.resize(coords.size());
//     surf_in.v_conn_tris.resize(coords.size());

//     for(int i=0;i<tris.size();i++){
//         int p1 = tris[i][0];
//         int p2 = tris[i][1];
//         int p3 = tris[i][2];
//         surf_in.v_conn_tris[p1].insert(i);
//         surf_in.v_conn_tris[p2].insert(i);
//         surf_in.v_conn_tris[p3].insert(i);

//         surf_in.v_conn_points[p1].insert(p2);
//         surf_in.v_conn_points[p1].insert(p3);

//         surf_in.v_conn_points[p2].insert(p1);
//         surf_in.v_conn_points[p2].insert(p3);

//         surf_in.v_conn_points[p3].insert(p1);
//         surf_in.v_conn_points[p3].insert(p2);
//     }

//     std::vector<TiGER::DiscreteGeometry::Edge> edges_out;

//     TiGER::discrete_geometry_repair::Tool_initializeData(surf_in);
//     auto start = std::chrono::high_resolution_clock::now();
//     TiGER::discrete_geometry_repair::TMeshSurf_identifyMultipleEdges(
//         surf_in,
//         edges_out  // 两个点id定一条边
//     );
//     auto finish = std::chrono::high_resolution_clock::now();
//     // 计算持续时间
//     std::chrono::duration<double> elapsed = finish - start;
//     // 输出运行时间
//     std::cout << "test_detectMultiedge time: " << elapsed.count() << " sec" << std::endl;
//     std::cout<<"multiedge size "<<edges_out.size()<<std::endl;

//     TiGER::Mesh mesh_edge;  
//     mesh_edge.Vertex = mesh.Vertex;
//     mesh_edge.Topo.resize(edges_out.size(),2);
//     for(int i=0;i<edges_out.size();i++){
//         for(int j=0;j<2;j++){
//             mesh_edge.Topo(i,j) = edges_out[i][j];
//         }
//     }
//     TiGER::MESHIO::writeVTK("C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/dirty_geometry/aircraft-14_stl_multiedge.vtk",mesh_edge);
// }

// TEST_CASE("tiger 1.5 TMeshSurf_identifyMultiplePoints", "test_detectMultipoint") {

//     TiGER::Mesh mesh;    
//     TiGER::MESHIO::readVTK("C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/dirty_geometry/aircraft-14_stl.vtk",mesh);
    
//     std::vector<std::array<double, 3>> coords;
//     std::vector<std::array<int, 3>> tris;
//     TiGER::matrix_to_list<3>(mesh.Vertex, coords);
//     TiGER::matrix_to_list<3>(mesh.Topo, tris);

//     TiGER::SurfaceMesh surf_in;
//     surf_in.coords = coords;
//     surf_in.tris = tris;
//     surf_in.v_conn_points.resize(coords.size());
//     surf_in.v_conn_tris.resize(coords.size());

//     for(int i=0;i<tris.size();i++){
//         int p1 = tris[i][0];
//         int p2 = tris[i][1];
//         int p3 = tris[i][2];
//         surf_in.v_conn_tris[p1].insert(i);
//         surf_in.v_conn_tris[p2].insert(i);
//         surf_in.v_conn_tris[p3].insert(i);

//         surf_in.v_conn_points[p1].insert(p2);
//         surf_in.v_conn_points[p1].insert(p3);

//         surf_in.v_conn_points[p2].insert(p1);
//         surf_in.v_conn_points[p2].insert(p3);

//         surf_in.v_conn_points[p3].insert(p1);
//         surf_in.v_conn_points[p3].insert(p2);
//     }

//     std::vector<TiGER::DiscreteGeometry::PointIndex> points_out;

//     TiGER::discrete_geometry_repair::Tool_initializeData(surf_in);
//     auto start = std::chrono::high_resolution_clock::now();
//     TiGER::discrete_geometry_repair::TMeshSurf_identifyMultiplePoints(
//         surf_in,
//         points_out
//     );
//     auto finish = std::chrono::high_resolution_clock::now();
//     // 计算持续时间
//     std::chrono::duration<double> elapsed = finish - start;
//     // 输出运行时间
//     std::cout << "test_detectMultipoint time: " << elapsed.count() << " sec" << std::endl;
//     std::cout<<"multipoint size "<<points_out.size()<<std::endl;

//     // for(int i=0;i<points_out.size();i++){
//     //     std::cout<<points_out[i]<<" ";
//     // }
//     // std::cout<<std::endl;

//     TiGER::Mesh mesh_point;
//     mesh_point.Vertex.resize(points_out.size(),3);
//     for(int i=0;i<points_out.size();i++){
//         for(int j=0;j<3;j++){
//             mesh_point.Vertex(i,j) = coords[points_out[i]][j];
//         }
//     }
//     TiGER::MESHIO::writeVTK("C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/dirty_geometry/aircraft-14_stl_multipoint.vtk",mesh_point);
// }

// TEST_CASE("tiger 1.5 test", "test_freeedge") {

//     TiGER::SurfaceMesh surf_in;
//     std::vector<TiGER::TiGER::DiscreteGeometry::Edge> edges_out;

//     std::array<double, 3> pt0 = {0,0,0};
//     std::array<double, 3> pt1 = {0,0,1};
//     std::array<double, 3> pt2 = {0,1,0};
//     // std::vector<std::array<int, 3>> tris;
//     std::array<int, 3> tri = {0,1,2};
//     // tris.push_back(tri);
//     surf_in.coords.push_back(pt0);
//     surf_in.coords.push_back(pt1);
//     surf_in.coords.push_back(pt2);
//     surf_in.tris.push_back(tri);
    

//     TiGER::discrete_geometry_repair::Tool_initializeData(surf_in);
//     TiGER::discrete_geometry_repair::TMeshSurf_identifyFreeEdges(
//         surf_in,
//         edges_out  // 两个点id定一条边
//     );

//     std::cout<<edges_out.size()<<std::endl;
// }

// TEST_CASE("tiger 1.5 merge_volume_mesh", "test_merge_volume_mesh") {
//     std::vector<TiGER::VolumeMesh> volume_in;
//     VolumeMesh volume_out;
//     for(int i=0;i<2;i++){
//         VolumeMesh mesh;
//         mesh.coords.resize(4);
//         mesh.coords[0] = {0.0,0.0,0.0};
//         mesh.coords[1] = {1.0,1.0,1.0};
//         mesh.coords[2] = {-1.0,2.5,0.8};
//         mesh.coords[3] = {0.0,1.5,-2.0};
//         mesh.tetras.resize(1);
//         mesh.tetras[0] = {0,1,2,3};
//         volume_in.push_back(mesh);
//     }
//     std::cout<<"int "<<std::endl;
//     std::cout<<volume_in[0].coords.size()<<std::endl;
//     std::cout<<volume_in[0].tetras.size()<<std::endl;
//     merge_volume_mesh(volume_in,volume_out);
//     std::cout<<"out "<<endl;
//     std::cout<<volume_out.coords.size()<<std::endl;
//     std::cout<<volume_out.tetras.size()<<std::endl;
// }

// TEST_CASE("tiger 1.5 Tool_projectPointToSegment", "test_project_PointOntoSegment"){
//     Eigen::RowVector3d p = {5.0,5.0,0.0};
//     Eigen::RowVector3d a = {0.0,0.0,0.0};
//     Eigen::RowVector3d b = {1.0,1.0,1.0};
//     Eigen::RowVector3d project_point;
//     double point_to_segment_distance;
//     TiGER::repair::Tool_projectPointToSegment(p,a,b,project_point,point_to_segment_distance);
//     std::cout<<"point ("<<project_point.x()<<","<<project_point.y()<<","<<project_point.z()<<" distance "<<point_to_segment_distance<<std::endl;
// }


// TEST_CASE("tiger 1.5 Tool_projectPointToTriangle", "test_project_PointOntoTriangle"){
//     Eigen::RowVector3d p = {5.0,5.0,0.0};
//     Eigen::RowVector3d a = {0.0,0.0,0.0};
//     Eigen::RowVector3d b = {1.0,1.0,1.0};
//     Eigen::RowVector3d c = {1.0,0.0,0.0};
//     Eigen::RowVector3d project_point;
//     double point_to_tri_distance;
//     TiGER::repair::Tool_projectPointToTriangle(p,a,b,c,project_point,point_to_tri_distance);
//     std::cout<<"point ("<<project_point.x()<<","<<project_point.y()<<","<<project_point.z()<<" distance "<<point_to_tri_distance<<std::endl;
// }

// TEST_CASE("tiger 1.5 Tool_intersectLines", "test_line_line_intersection"){
//     Eigen::RowVector3d e1 = {0.0,0.0,1.0};
//     Eigen::RowVector3d e2 = {1.0,1.0,0.0};
//     Eigen::RowVector3d e3 = {0.0,0.0,0.0};
//     Eigen::RowVector3d e4 = {5.0,5.0,5.0};
//     // Eigen::RowVector3d intersection_point;
//     std::vector<Eigen::RowVector3d> intersection_points;
//     TiGER::repair::intersectionresult result = TiGER::repair::Tool_intersectLines(e1,e2,e3,e4,intersection_points);
    
//     std::cout<<result<<std::endl;;
//     for(auto intersection_point : intersection_points)
//         std::cout<<" intersection_point ("<<intersection_point.x()<<","<<intersection_point.y()<<","<<intersection_point.z()<<")"<<std::endl;

//     e1 = {0.0,0.0,0.0};
//     e2 = {5.0,5.0,0.0};
//     e3 = {1.0,1.0,0.0};
//     e4 = {3.0,3.0,0.0};
//     intersection_points.clear();
//     result = TiGER::repair::Tool_intersectLines(e1,e2,e3,e4,intersection_points);
    
//     std::cout<<result<<std::endl;;
//     for(auto intersection_point : intersection_points)
//         std::cout<<" intersection_point ("<<intersection_point.x()<<","<<intersection_point.y()<<","<<intersection_point.z()<<")"<<std::endl;

//     e1 = {0.0,0.0,1.0};
//     e2 = {1.0,1.0,0.0};
//     e3 = {0.0,0.0,1.0};
//     e4 = {5.0,5.0,5.0};
//     intersection_points.clear();
//     result = TiGER::repair::Tool_intersectLines(e1,e2,e3,e4,intersection_points);
    
//     std::cout<<result<<std::endl;;
//     for(auto intersection_point : intersection_points)
//         std::cout<<" intersection_point ("<<intersection_point.x()<<","<<intersection_point.y()<<","<<intersection_point.z()<<")"<<std::endl;

//     e1 = {1.0,0.0,0.0};
//     e2 = {0.0,1.0,0.0};
//     e3 = {0.0,0.0,0.0};
//     e4 = {5.0,5.0,5.0};
//     intersection_points.clear();
//     result = TiGER::repair::Tool_intersectLines(e1,e2,e3,e4,intersection_points);
    
//     std::cout<<result<<std::endl;;
//     for(auto intersection_point : intersection_points)
//         std::cout<<" intersection_point ("<<intersection_point.x()<<","<<intersection_point.y()<<","<<intersection_point.z()<<")"<<std::endl;
    
//     e1 = {1.0,0.0,0.0};
//     e2 = {1.0,1.0,0.0};
//     e3 = {0.0,0.0,0.0};
//     e4 = {5.0,5.0,0.0};
//     intersection_points.clear();
//     result = TiGER::repair::Tool_intersectLines(e1,e2,e3,e4,intersection_points);
    
//     std::cout<<result<<std::endl;;
//     for(auto intersection_point : intersection_points)
//         std::cout<<" intersection_point ("<<intersection_point.x()<<","<<intersection_point.y()<<","<<intersection_point.z()<<")"<<std::endl;

//     e1 = {0.0,0.0,0.0};
//     e2 = {2.0,2.0,0.0};
//     e3 = {1.0,1.0,0.0};
//     e4 = {5.0,5.0,0.0};
//     intersection_points.clear();
//     result = TiGER::repair::Tool_intersectLines(e1,e2,e3,e4,intersection_points);
    
//     std::cout<<result<<std::endl;;
//     for(auto intersection_point : intersection_points)
//         std::cout<<" intersection_point ("<<intersection_point.x()<<","<<intersection_point.y()<<","<<intersection_point.z()<<")"<<std::endl;

// }

// TEST_CASE("tiger 1.5 Tool_intersectLineWithTriangle", "test_line_tri_intersection"){
//     //
//     Eigen::RowVector3d e1 = {1.0,0.0,0.0};
//     Eigen::RowVector3d e2 = {1.0,2.0,0.0};
//     Eigen::RowVector3d a = {0.0,0.0,0.0};
//     Eigen::RowVector3d b = {1.0,1.0,0.0};
//     Eigen::RowVector3d c = {2.0,0.0,0.0};

//     std::vector<Eigen::RowVector3d> intersection_points;
//     TiGER::repair::intersectionresult result = TiGER::repair::Tool_intersectLineWithTriangle(e1,e2,a,b,c,intersection_points);
//     std::cout<<result<<std::endl;;
//     for(auto intersection_point : intersection_points)
//         std::cout<<" intersection_point ("<<intersection_point.x()<<","<<intersection_point.y()<<","<<intersection_point.z()<<")"<<std::endl;
    

//     e1 = {1.0,1.0,0.0};
//     e2 = {1.0,-2.0,0.0};
//     intersection_points.clear();
//     result = TiGER::repair::Tool_intersectLineWithTriangle(e1,e2,a,b,c,intersection_points);
//     std::cout<<result<<std::endl;;
//     for(auto intersection_point : intersection_points)
//         std::cout<<" intersection_point ("<<intersection_point.x()<<","<<intersection_point.y()<<","<<intersection_point.z()<<")"<<std::endl;
    

//     e1 = {-1.0,0.0,0.0};
//     e2 = {3.0,0.0,0.0};
//     intersection_points.clear();
//     result = TiGER::repair::Tool_intersectLineWithTriangle(e1,e2,a,b,c,intersection_points);
//     std::cout<<result<<std::endl;;
//     for(auto intersection_point : intersection_points)
//         std::cout<<" intersection_point ("<<intersection_point.x()<<","<<intersection_point.y()<<","<<intersection_point.z()<<")"<<std::endl;
    
//     e1 = {1.0,2.0,0.0};
//     e2 = {1.0,-1.0,0.0};
//     intersection_points.clear();
//     result = TiGER::repair::Tool_intersectLineWithTriangle(e1,e2,a,b,c,intersection_points);
//     std::cout<<result<<std::endl;;
//     for(auto intersection_point : intersection_points)
//         std::cout<<" intersection_point ("<<intersection_point.x()<<","<<intersection_point.y()<<","<<intersection_point.z()<<")"<<std::endl;
    
//     e1 = {0.0,-1.0,0.0};
//     e2 = {2.0,1.0,0.0};
//     intersection_points.clear();
//     result = TiGER::repair::Tool_intersectLineWithTriangle(e1,e2,a,b,c,intersection_points);
//     std::cout<<result<<std::endl;;
//     for(auto intersection_point : intersection_points)
//         std::cout<<" intersection_point ("<<intersection_point.x()<<","<<intersection_point.y()<<","<<intersection_point.z()<<")"<<std::endl;
    
//     e1 = {1.0,0.5,0.0};
//     e2 = {3.0,-0.5,0.0};
//     intersection_points.clear();
//     result = TiGER::repair::Tool_intersectLineWithTriangle(e1,e2,a,b,c,intersection_points);
//     std::cout<<result<<std::endl;;
//     for(auto intersection_point : intersection_points)
//         std::cout<<" intersection_point ("<<intersection_point.x()<<","<<intersection_point.y()<<","<<intersection_point.z()<<")"<<std::endl;
    
//     e1 = {1.0,0.5,0.0};
//     e2 = {1.0,-1.0,0.0};
//     intersection_points.clear();
//     result = TiGER::repair::Tool_intersectLineWithTriangle(e1,e2,a,b,c,intersection_points);
//     std::cout<<result<<std::endl;;
//     for(auto intersection_point : intersection_points)
//         std::cout<<" intersection_point ("<<intersection_point.x()<<","<<intersection_point.y()<<","<<intersection_point.z()<<")"<<std::endl;
    
//     e1 = {1.0,0.5,-1.0};
//     e2 = {1.0,0.5,1.0};
//     intersection_points.clear();
//     result = TiGER::repair::Tool_intersectLineWithTriangle(e1,e2,a,b,c,intersection_points);
//     std::cout<<result<<std::endl;;
//     for(auto intersection_point : intersection_points)
//         std::cout<<" intersection_point ("<<intersection_point.x()<<","<<intersection_point.y()<<","<<intersection_point.z()<<")"<<std::endl;
    
// }

// TEST_CASE("tiger 1.5 intersectionFaceSets", "test_intersectionFaceSets") {
//     TiGER::Mesh mesh;    
//     TiGER::MESHIO::readVTK("/mnt/c/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/myexample/in.vtk",mesh,"surface_id");
    
//     TiGER::SurfaceMesh surf_in_A;
//     TiGER::matrix_to_list<3>(mesh.Vertex, surf_in_A.coords);
//     TiGER::matrix_to_list<3>(mesh.Topo, surf_in_A.tris);
//     TiGER::matrix_to_list<1>(mesh.Masks,surf_in_A.attribute_int);
//     std::cout<<"point size "<<surf_in_A.coords.size()<<std::endl;
//     std::cout<<"tri size "<<surf_in_A.tris.size()<<" "<<surf_in_A.tris[0].size()<<std::endl; 
//     std::cout<<surf_in_A.attribute_int[0]<<" "<<surf_in_A.attribute_int.back()<<std::endl;
    
//     TiGER::Mesh mesh1;    
//     TiGER::MESHIO::readVTK("/mnt/c/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/myexample/out.vtk",mesh1,"surface_id");
//     TiGER::SurfaceMesh surf_in_B;
//     TiGER::matrix_to_list<3>(mesh1.Vertex, surf_in_B.coords);
//     TiGER::matrix_to_list<3>(mesh1.Topo, surf_in_B.tris);
//     TiGER::matrix_to_list<1>(mesh1.Masks,surf_in_B.attribute_int);
//     std::cout<<"point size "<<surf_in_B.coords.size()<<std::endl;
//     std::cout<<"tri size "<<surf_in_B.tris.size()<<std::endl; 

//     TiGER::SurfaceMesh surf_out;
//     std::vector<TiGER::discrete_geometry_repair::InterEdge> edge_out;
//     TiGER::RepairParameters args;

//     args.seteps_bool_option(true);
//     args.setbool_eps(1e-4);

//     TiGER::discrete_geometry_repair::intersectionFaceSets(
//         surf_in_A,
//         surf_in_B,
//         args,
//         surf_out,
//         edge_out  // 抽取交线
//     );

//     std::set<int> face_in;
//     std::string filename = "/mnt/c/Users/10170/Desktop/face.txt";
//     FILE *f = fopen(filename.c_str(), "w");
//     for(int i=0;i<surf_out.tris.size();i++){
//         if(face_in.find(surf_out.attribute_int[i]) == face_in.end()){
//             fprintf(f,"faceId %d \n",surf_out.attribute_int[i]);
//             for(int k=0;k<surf_out.parent_face[i].size();k++)
//                 fprintf(f, "%d ", surf_out.parent_face[i][k]);
//             fprintf(f,"\n");
//             face_in.insert(surf_out.attribute_int[i]);
//         }
//     }

//     TiGER::Mesh mesh_tri;
//     TiGER::list_to_matrix<3>(surf_out.coords,mesh_tri.Vertex);
//     TiGER::list_to_matrix<3>(surf_out.tris,mesh_tri.Topo);
//     TiGER::list_to_matrix<1>(surf_out.attribute_int,mesh_tri.Masks);
//     TiGER::MESHIO::writeVTK("/mnt/c/Users/10170/Desktop/tri.vtk",mesh_tri);

//     //交线
//     std::cout<<"edges size "<<edge_out.size()<<std::endl;
//     std::map<int,std::vector<std::array<int,2>>> edges;
//     for(auto& it : edge_out){
//         std::array<int,2> edge;
//         edge[0] = it.startid;
//         edge[1] = it.endid;
//         // std::cout<<" -- --- "<<std::endl;
//         for(int id:it.faceid){
//             // std::cout<<" !! "<<id<<std::endl;
//             edges[id].push_back(edge);
//         }
//     }

//     for(auto& e : edges){
//         TiGER::Mesh mesh_edge;  
//         TiGER::list_to_matrix<3>(surf_out.coords,mesh_edge.Vertex);
//         mesh_edge.Topo.resize(e.second.size(),2);
//         for(int i=0;i<e.second.size();i++){
//             for(int j=0;j<2;j++){
//                 mesh_edge.Topo(i,j) = e.second[i][j];
//             }
//         }
//         std::string filename = std::to_string(e.first) + ".vtk";
//         TiGER::MESHIO::writeVTK("/mnt/c/Users/10170/Desktop/"+filename,mesh_edge);
//     }
// }


// TEST_CASE("tiger 1.5 intersectionFaces", "test_intersectionFaces") {
//     TiGER::Mesh mesh;    
//     TiGER::MESHIO::readVTK("/mnt/c/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/error/big_cubic.o.vtk",mesh);
    
//     TiGER::SurfaceMesh surf_in_A;
//     TiGER::matrix_to_list<3>(mesh.Vertex, surf_in_A.coords);
//     TiGER::matrix_to_list<3>(mesh.Topo, surf_in_A.tris);
//     TiGER::matrix_to_list<1>(mesh.Masks,surf_in_A.attribute_int);
//     std::cout<<"point size "<<surf_in_A.coords.size()<<std::endl;
//     std::cout<<"tri size "<<surf_in_A.tris.size()<<" "<<surf_in_A.tris[0].size()<<std::endl; 
//     std::cout<<surf_in_A.tris[0][0]<<" "<<surf_in_A.tris[0][1]<<" "<<surf_in_A.tris[0][2]<<std::endl;
    
//     TiGER::Mesh mesh1;    
//     TiGER::MESHIO::readVTK("/mnt/c/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/error/small_cubic.o.vtk",mesh1);
//     TiGER::SurfaceMesh surf_in_B;
//     TiGER::matrix_to_list<3>(mesh1.Vertex, surf_in_B.coords);
//     TiGER::matrix_to_list<3>(mesh1.Topo, surf_in_B.tris);
//     TiGER::matrix_to_list<1>(mesh1.Masks,surf_in_B.attribute_int);
//     std::cout<<"point size "<<surf_in_B.coords.size()<<std::endl;
//     std::cout<<"tri size "<<surf_in_B.tris.size()<<std::endl; 

//     TiGER::SurfaceMesh surf_out;
//     std::vector<TiGER::DiscreteGeometry::Edge> edge_out;
//     TiGER::RepairParameters args;

//     TiGER::discrete_geometry_repair::intersectionFaces_new(surf_in_A,surf_in_B,args,surf_out,edge_out);

//     TiGER::Mesh mesh_tri;
//     TiGER::list_to_matrix<3>(surf_out.coords,mesh_tri.Vertex);
//     TiGER::list_to_matrix<3>(surf_out.tris,mesh_tri.Topo);
//     TiGER::list_to_matrix<1>(surf_out.attribute_int,mesh_tri.Masks);
//     TiGER::MESHIO::writeVTK("/mnt/c/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/error/cubic_tri.vtk",mesh_tri);

//     TiGER::Mesh mesh_edge;  
//     TiGER::list_to_matrix<3>(surf_out.coords,mesh_edge.Vertex);
//     mesh_edge.Topo.resize(edge_out.size(),2);
//     for(int i=0;i<edge_out.size();i++){
//         for(int j=0;j<2;j++){
//             mesh_edge.Topo(i,j) = edge_out[i][j];
//         }
//     }
//     TiGER::MESHIO::writeVTK("/mnt/c/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/error/cubic_edge.vtk",mesh_edge);
// }


// TEST_CASE("tiger 1.5 project_SegmentOntoTriangle", "test_project_SegmentOntoTriangle") {
//     std::cout<<"example1"<<endl;
//     Eigen::RowVector3d e1 = {2.0,0.0,0.0};
//     Eigen::RowVector3d e2 = {0.0,1.0,0.0};
//     Eigen::RowVector3d a = {0.0,0.0,0.0};
//     Eigen::RowVector3d b = {1.0,1.0,0.0};
//     Eigen::RowVector3d c = {1.0,0.0,0.0};
//     std::vector<Eigen::RowVector3d> project_edge;
//     TiGER::repair::project_SegmentOntoTriangle(e1,e2,a,b,c,project_edge);
//     for(auto point : project_edge){
//         std::cout<<" point ("<<point.x()<<","<<point.y()<<","<<point.z()<<")"<<std::endl;
//     }
//     std::cout<<"example2"<<endl;
//     e1 = {0.5,0.4,0.0};
//     e2 = {0.0,1.0,0.0};
//     a = {0.0,0.0,0.0};
//     b = {1.0,1.0,0.0};
//     c = {1.0,0.0,0.0};
//     std::vector<Eigen::RowVector3d> project_edge2;
//     TiGER::repair::project_SegmentOntoTriangle(e1,e2,a,b,c,project_edge2);
//     for(auto point : project_edge2){
//         std::cout<<" point ("<<point.x()<<","<<point.y()<<","<<point.z()<<")"<<std::endl;
//     }
//     std::cout<<"example3"<<endl;
//     e1 = {1.0,1.0,0.0};
//     e2 = {0.0,1.0,0.0};
//     a = {0.0,0.0,0.0};
//     b = {1.0,1.0,0.0};
//     c = {1.0,0.0,0.0};
//     std::vector<Eigen::RowVector3d> project_edge3;
//     TiGER::repair::project_SegmentOntoTriangle(e1,e2,a,b,c,project_edge3);
//     for(auto point : project_edge3){
//         std::cout<<" point ("<<point.x()<<","<<point.y()<<","<<point.z()<<")"<<std::endl;
//     }
//     std::cout<<"example4"<<endl;
//     e1 = {0.5,0.0,0.0};
//     e2 = {2.0,0.0,0.0};
//     a = {0.0,0.0,0.0};
//     b = {1.0,1.0,0.0};
//     c = {1.0,0.0,0.0};
//     std::vector<Eigen::RowVector3d> project_edge4;
//     TiGER::repair::project_SegmentOntoTriangle(e1,e2,a,b,c,project_edge4);
//     for(auto point : project_edge4){
//         std::cout<<" point ("<<point.x()<<","<<point.y()<<","<<point.z()<<")"<<std::endl;
//     }
//     std::cout<<"example5"<<endl;
//     e1 = {0.5,0.5,0.0};
//     e2 = {0.5,0.5,1.0};
//     a = {0.0,0.0,0.0};
//     b = {1.0,1.0,0.0};
//     c = {1.0,0.0,0.0};
//     std::vector<Eigen::RowVector3d> project_edge5;
//     TiGER::repair::project_SegmentOntoTriangle(e1,e2,a,b,c,project_edge5);
//     for(auto point : project_edge5){
//         std::cout<<" point ("<<point.x()<<","<<point.y()<<","<<point.z()<<")"<<std::endl;
//     }

// }


// TEST_CASE("tiger 1.5 partition", "test_partition") {

//     Eigen::RowVector3d b = {0.0,0.0,0.0};
//     Eigen::RowVector3d c = {2.0,0.0,0.0};
//     Eigen::RowVector3d a = {1.0,1.0,0.0};

//     std::vector<Eigen::RowVector3d> project_poins;
//     std::array<int,2> vertextid;
//     std::array<int,2> edgeid;
//     for(int i=0;i<2;i++){
//         vertextid[i] = -1;
//         edgeid[i] = -1;
//     }
//     project_poins.push_back(Eigen::RowVector3d(2.0,0.0,0.0));
//     vertextid[0] = 2;
//     TiGER::SurfaceMesh surf_out;
//     TiGER::repair::partition(
//         a,b,c,
//         project_poins,
//         vertextid,
//         edgeid,
//         surf_out
//     );
//     TiGER::Mesh mesh_tri;
//     TiGER::list_to_matrix<3>(surf_out.coords,mesh_tri.Vertex);
//     TiGER::list_to_matrix<3>(surf_out.tris,mesh_tri.Topo);
//     TiGER::MESHIO::writeVTK("C:/Users/10170/Desktop/partion.vtk",mesh_tri);

//     project_poins.clear();
//     project_poins.push_back(Eigen::RowVector3d(1,0.5,0.0));
//     vertextid[0] = -1;
//     TiGER::SurfaceMesh surf_out1;
//     TiGER::repair::partition(
//         a,b,c,
//         project_poins,
//         vertextid,
//         edgeid,
//         surf_out1
//     );
//     TiGER::Mesh mesh_tri1;
//     TiGER::list_to_matrix<3>(surf_out1.coords,mesh_tri1.Vertex);
//     TiGER::list_to_matrix<3>(surf_out1.tris,mesh_tri1.Topo);
//     TiGER::MESHIO::writeVTK("C:/Users/10170/Desktop/partion1.vtk",mesh_tri1);

//     project_poins.clear();
//     project_poins.push_back(Eigen::RowVector3d(0.5,0.5,0.0));
//     edgeid[0] = 0;
//     TiGER::SurfaceMesh surf_out2;
//     TiGER::repair::partition(
//         a,b,c,
//         project_poins,
//         vertextid,
//         edgeid,
//         surf_out2
//     );
//     TiGER::Mesh mesh_tri2;
//     TiGER::list_to_matrix<3>(surf_out2.coords,mesh_tri2.Vertex);
//     TiGER::list_to_matrix<3>(surf_out2.tris,mesh_tri2.Topo);
//     TiGER::MESHIO::writeVTK("C:/Users/10170/Desktop/partion2.vtk",mesh_tri2);

//     project_poins.clear();
//     project_poins.push_back(Eigen::RowVector3d(0.0,0.0,0.0));
//     project_poins.push_back(Eigen::RowVector3d(2.0,0.0,0.0));
//     edgeid[0] = -1;
//     vertextid[0] = 1;
//     vertextid[1] = 2;
//     TiGER::SurfaceMesh surf_out3;
//     TiGER::repair::partition(
//         a,b,c,
//         project_poins,
//         vertextid,
//         edgeid,
//         surf_out3
//     );
//     TiGER::Mesh mesh_tri3;
//     TiGER::list_to_matrix<3>(surf_out3.coords,mesh_tri3.Vertex);
//     TiGER::list_to_matrix<3>(surf_out3.tris,mesh_tri3.Topo);
//     TiGER::MESHIO::writeVTK("C:/Users/10170/Desktop/partion3.vtk",mesh_tri3);

//     project_poins.clear();
//     project_poins.push_back(Eigen::RowVector3d(0.5,0.3,0.0));
//     project_poins.push_back(Eigen::RowVector3d(1.0,0.8,0.0));
//     vertextid[0] = -1;
//     vertextid[1] = -1;
//     TiGER::SurfaceMesh surf_out4;
//     TiGER::repair::partition(
//         a,b,c,
//         project_poins,
//         vertextid,
//         edgeid,
//         surf_out4
//     );
//     TiGER::Mesh mesh_tri4;
//     TiGER::list_to_matrix<3>(surf_out4.coords,mesh_tri4.Vertex);
//     TiGER::list_to_matrix<3>(surf_out4.tris,mesh_tri4.Topo);
//     TiGER::MESHIO::writeVTK("C:/Users/10170/Desktop/partion4.vtk",mesh_tri4);

//     project_poins.clear();
//     project_poins.push_back(Eigen::RowVector3d(0.5,0.0,0.0));
//     project_poins.push_back(Eigen::RowVector3d(1.0,0.0,0.0));
//     edgeid[0] = 1;
//     edgeid[1] = 1;
//     TiGER::SurfaceMesh surf_out5;
//     TiGER::repair::partition(
//         a,b,c,
//         project_poins,
//         vertextid,
//         edgeid,
//         surf_out5
//     );
//     TiGER::Mesh mesh_tri5;
//     TiGER::list_to_matrix<3>(surf_out5.coords,mesh_tri5.Vertex);
//     TiGER::list_to_matrix<3>(surf_out5.tris,mesh_tri5.Topo);
//     TiGER::MESHIO::writeVTK("C:/Users/10170/Desktop/partion5.vtk",mesh_tri5);

//     project_poins.clear();
//     project_poins.push_back(Eigen::RowVector3d(0.5,0.5,0.0));
//     project_poins.push_back(Eigen::RowVector3d(1.8,0.2,0.0));
//     edgeid[0] = 0;
//     edgeid[1] = 2;
//     TiGER::SurfaceMesh surf_out6;
//     TiGER::repair::partition(
//         a,b,c,
//         project_poins,
//         vertextid,
//         edgeid,
//         surf_out6
//     );
//     TiGER::Mesh mesh_tri6;
//     TiGER::list_to_matrix<3>(surf_out6.coords,mesh_tri6.Vertex);
//     TiGER::list_to_matrix<3>(surf_out6.tris,mesh_tri6.Topo);
//     TiGER::MESHIO::writeVTK("C:/Users/10170/Desktop/partion6.vtk",mesh_tri6);

//     project_poins.clear();
//     project_poins.push_back(Eigen::RowVector3d(0.5,0.5,0.0));
//     project_poins.push_back(Eigen::RowVector3d(1.5,0.2,0.0));
//     edgeid[0] = 0;
//     edgeid[1] = -1;
//     TiGER::SurfaceMesh surf_out7;
//     TiGER::repair::partition(
//         a,b,c,
//         project_poins,
//         vertextid,
//         edgeid,
//         surf_out7
//     );
//     TiGER::Mesh mesh_tri7;
//     TiGER::list_to_matrix<3>(surf_out7.coords,mesh_tri7.Vertex);
//     TiGER::list_to_matrix<3>(surf_out7.tris,mesh_tri7.Topo);
//     TiGER::MESHIO::writeVTK("C:/Users/10170/Desktop/partion7.vtk",mesh_tri7);

//     project_poins.clear();
//     project_poins.push_back(Eigen::RowVector3d(0.5,0.3,0.0));
//     project_poins.push_back(Eigen::RowVector3d(2.0,0.0,0.0));
//     edgeid[0] = -1;
//     vertextid[1] = 2;
//     TiGER::SurfaceMesh surf_out8;
//     TiGER::repair::partition(
//         a,b,c,
//         project_poins,
//         vertextid,
//         edgeid,
//         surf_out8
//     );
//     TiGER::Mesh mesh_tri8;
//     TiGER::list_to_matrix<3>(surf_out8.coords,mesh_tri8.Vertex);
//     TiGER::list_to_matrix<3>(surf_out8.tris,mesh_tri8.Topo);
//     TiGER::MESHIO::writeVTK("C:/Users/10170/Desktop/partion8.vtk",mesh_tri8);
    
//     project_poins.clear();
//     project_poins.push_back(Eigen::RowVector3d(0.5,0.5,0.0));
//     project_poins.push_back(Eigen::RowVector3d(2.0,0.0,0.0));
//     edgeid[0] = 0;
//     vertextid[1] = 2;
//     TiGER::SurfaceMesh surf_out9;
//     TiGER::repair::partition(
//         a,b,c,
//         project_poins,
//         vertextid,
//         edgeid,
//         surf_out9
//     );
//     TiGER::Mesh mesh_tri9;
//     TiGER::list_to_matrix<3>(surf_out9.coords,mesh_tri9.Vertex);
//     TiGER::list_to_matrix<3>(surf_out9.tris,mesh_tri9.Topo);
//     TiGER::MESHIO::writeVTK("C:/Users/10170/Desktop/partion9.vtk",mesh_tri9);
// }


// TEST_CASE("tiger 1.5 print", "test_prinnt"){
//     Eigen::RowVector3d e1 = {2.0,0.0,1.0};
//     Eigen::RowVector3d e2 = {-2.0,0.5,1.0};

//     TiGER::Mesh mesh_edge;  
//     mesh_edge.Vertex.resize(2,3);
//     mesh_edge.Vertex(0,0) = 2.0;
//     mesh_edge.Vertex(0,1) = 0.0;
//     mesh_edge.Vertex(0,2) = 0.0;
//     mesh_edge.Vertex(1,0) = -2.0;
//     mesh_edge.Vertex(1,1) = 0.5;
//     mesh_edge.Vertex(1,2) = 1.0;
//     mesh_edge.Topo.resize(1,2);
//     mesh_edge.Topo(0,0) = 0;
//     mesh_edge.Topo(0,1) = 1;

//     TiGER::MESHIO::writeVTK("C:/Users/10170/Desktop/segment.vtk",mesh_edge);

//     // Eigen::RowVector3d a = {0.0,0.0,0.0};
//     // Eigen::RowVector3d b = {0.0,1.0,0.0};
//     // Eigen::RowVector3d c = {1.0,0.0,0.0};

//     TiGER::SurfaceMesh surf_in;
//     std::array<double, 3> pt0 = {0,0,0};
//     std::array<double, 3> pt1 = {0,1,0};
//     std::array<double, 3> pt2 = {1,0,0};
//     std::array<double, 3> pt3 = {-0.5,0.5,-0.5};
//     std::array<int, 3> tri = {0,1,2};
//     std::array<int, 3> tri1 = {0,1,3};
//     surf_in.coords.push_back(pt0);
//     surf_in.coords.push_back(pt1);
//     surf_in.coords.push_back(pt2);
//     surf_in.coords.push_back(pt3);
//     surf_in.tris.push_back(tri);
//     surf_in.tris.push_back(tri1);
    
//     std::vector<TiGER::SurfaceMesh> surf_out;
//     for(int i=0;i<surf_in.tris.size();i++){
//         // std::cout<<"("<<surf_in.coords[surf_in.tris[i][0]][0]<<","<<surf_in.coords[surf_in.tris[i][0]][1]<<","<<surf_in.coords[surf_in.tris[i][0]][2]<<")"<<std::endl;
//         Eigen::RowVector3d a = Eigen::RowVector3d(surf_in.coords[surf_in.tris[i][0]][0],
//                                 surf_in.coords[surf_in.tris[i][0]][1],
//                                 surf_in.coords[surf_in.tris[i][0]][2]);
//         // std::cout<<"("<<surf_in.coords[surf_in.tris[i][1]][0]<<","<<surf_in.coords[surf_in.tris[i][1]][1]<<","<<surf_in.coords[surf_in.tris[i][1]][2]<<")"<<std::endl;                        
//         Eigen::RowVector3d b = Eigen::RowVector3d(surf_in.coords[surf_in.tris[i][1]][0],
//                                 surf_in.coords[surf_in.tris[i][1]][1],
//                                 surf_in.coords[surf_in.tris[i][1]][2]);
//         // std::cout<<"("<<surf_in.coords[surf_in.tris[i][2]][0]<<","<<surf_in.coords[surf_in.tris[i][2]][1]<<","<<surf_in.coords[surf_in.tris[i][2]][2]<<")"<<std::endl;
//         Eigen::RowVector3d c = Eigen::RowVector3d(surf_in.coords[surf_in.tris[i][2]][0],
//                                 surf_in.coords[surf_in.tris[i][2]][1],
//                                 surf_in.coords[surf_in.tris[i][2]][2]);

//         std::vector<Eigen::RowVector3d> project_edge;
//         // 投影
//         TiGER::repair::project_SegmentOntoTriangle(
//             e1,e2,
//             a,b,c,
//             project_edge);
//         if(project_edge.size() == 0)
//             continue;
//         // std::cout<<i<<" "<<project_edge[0](0)<<" "<<project_edge[0](1)<<" "<<project_edge[0](2)<<std::endl;
//         // 维护
//         std::array<int,2> vertextid{-1,-1};
//         std::array<int,2> edgeid{-1,-1};
//         Eigen::RowVector3d project_point;
//         double dis1 = -1;
//         TiGER::repair::project_PointOntoSegment(
//             project_edge[0],
//             a,b,
//             project_point,
//             dis1);
//         double dis2 = -1;
//         TiGER::repair::project_PointOntoSegment(
//             project_edge[0],
//             b,c,
//             project_point,
//             dis2);
//         double dis3 = -1;
//         TiGER::repair::project_PointOntoSegment(
//             project_edge[0],
//             c,a,
//             project_point,
//             dis3);
//         if(dis1 == 0){  // ab
//             if(dis2 == 0){  // bc
//                 vertextid[0] = 1;
//             }
//             else if(dis3 == 0){  // ca
//                 vertextid[0] = 0;
//             }
//             else{
//                 edgeid[0] = 0;
//             }
//         }
//         else if(dis2 == 0){
//             if(dis3 == 0){
//                 vertextid[0] = 2;
//             }
//             else{
//                 edgeid[0] = 1;
//             }
//         }
//         else if(dis3 == 0){
//             edgeid[0] = 2;
//         }

//         if(project_edge.size() == 2){
//             // std::cout<<i<<" "<<project_edge[1](0)<<" "<<project_edge[1](1)<<" "<<project_edge[1](2)<<std::endl;
//             Eigen::RowVector3d project_point1;
//             double dis3 = -1;
//             TiGER::repair::project_PointOntoSegment(
//                 project_edge[1],
//                 a,b,
//                 project_point1,
//                 dis3);
//             double dis4 = -1;
//             TiGER::repair::project_PointOntoSegment(
//                 project_edge[1],
//                 b,c,
//                 project_point1,
//                 dis4);
//             double dis5 = -1;
//             TiGER::repair::project_PointOntoSegment(
//                 project_edge[1],
//                 c,a,
//                 project_point1,
//                 dis5);
//             if(dis3 == 0){  // ab
//                 if(dis4 == 0){  // bc
//                     vertextid[1] = 1;
//                 }
//                 else if(dis5 == 0){  // ca
//                     vertextid[1] = 0;
//                 }
//                 else{
//                     edgeid[1] = 0;
//                 }
//             }
//             else if(dis4 == 0){
//                 if(dis5 == 0){
//                     vertextid[1] = 2;
//                 }
//                 else{
//                     edgeid[1] = 1;
//                 }
//             }
//             else if(dis5 == 0){
//                 edgeid[1] = 2;
//             }
//         }
//         // std::cout<<i<<" "<<vertextid[0]<<" "<<vertextid[1]<<std::endl;
//         // std::cout<<i<<" "<<edgeid[0]<<" "<<edgeid[1]<<std::endl;

//         // // 分割
//         TiGER::SurfaceMesh piece;
//         TiGER::repair::partition(
//             a,b,c,
//             project_edge,
//             vertextid,
//             edgeid,
//             piece
//         );
//         // std::cout<<piece.coords.size()<<" "<<piece.tris.size()<<std::endl;
//         surf_out.push_back(piece);
//     }

//     TiGER::SurfaceMesh result;
//     TiGER::merge_surface_mesh(
//         surf_out,
//         result
//     );

//     // std::cout<<result.coords.size()<<" "<<result.tris.size()<<std::endl;

//     TiGER::Mesh mesh_tri;
//     TiGER::list_to_matrix<3>(result.coords,mesh_tri.Vertex);
//     TiGER::list_to_matrix<3>(result.tris,mesh_tri.Topo);

//     TiGER::MESHIO::writeVTK("C:/Users/10170/Desktop/print.vtk",mesh_tri);
// }


// TEST_CASE("tiger 1.5 project_SegmentOntoSegment", "project_SegmentOntoSegment"){
//     Eigen::RowVector3d a = {0.0,0.0,0.0};
//     Eigen::RowVector3d b = {0.0,1.0,0.0};
//     Eigen::RowVector3d c = {1.0,0.0,0.0};
//     Eigen::RowVector3d d = {-1.0,1.0,5.0};
//     Eigen::RowVector3d closestPoint1;
//     Eigen::RowVector3d closestPoint2;
//     double dis = TiGER::repair::project_SegmentOntoSegment(a,b,c,d,closestPoint1,closestPoint2);
//     std::cout<<dis<<std::endl;
//     std::cout<<closestPoint1.x()<<" "<<closestPoint1.y()<<" "<<closestPoint1.z()<<std::endl;
//     std::cout<<closestPoint2.x()<<" "<<closestPoint2.y()<<" "<<closestPoint2.z()<<std::endl;
// }

// TEST_CASE("tiger 1.5 print", "print"){
//     Eigen::RowVector3d a = {-1.0,0.0,0.0};
//     Eigen::RowVector3d b = {0.0,2.0,0.0};
//     Eigen::RowVector3d c = {1.0,0.0,0.0};

//     Eigen::RowVector3d x = {-2.0,1.0,5.0};
//     Eigen::RowVector3d y = {1.0,-1.0,2.0};

//     std::vector<Eigen::RowVector3d> points;
//     std::vector<double> diss;
    
//     Eigen::RowVector3d closestPoint1;
//     Eigen::RowVector3d closestPoint2;
//     double dis = TiGER::repair::project_SegmentOntoSegment(a,b,x,y,closestPoint1,closestPoint2);
//     std::cout<<" -- "<<closestPoint2.x()<<" "<<closestPoint2.y()<<" "<<closestPoint2.z()<<std::endl;
//     if(!(closestPoint2 == x || closestPoint2 == y)){
//         std::cout<<closestPoint1.x()<<" "<<closestPoint1.y()<<" "<<closestPoint1.z()<<std::endl;
//         points.push_back(closestPoint1);
//         diss.push_back(dis);
//     }

//     dis = TiGER::repair::project_SegmentOntoSegment(b,c,x,y,closestPoint1,closestPoint2);
//     std::cout<<" -- "<<closestPoint2.x()<<" "<<closestPoint2.y()<<" "<<closestPoint2.z()<<std::endl;
//     if(!(closestPoint2 == x || closestPoint2 == y)){
//         std::cout<<closestPoint1.x()<<" "<<closestPoint1.y()<<" "<<closestPoint1.z()<<std::endl;
//         points.push_back(closestPoint1);
//         diss.push_back(dis);
//     }

//     dis = TiGER::repair::project_SegmentOntoSegment(c,a,x,y,closestPoint1,closestPoint2);
//     std::cout<<" -- "<<closestPoint2.x()<<" "<<closestPoint2.y()<<" "<<closestPoint2.z()<<std::endl;
//     if(!(closestPoint2 == x || closestPoint2 == y)){
//         std::cout<<closestPoint1.x()<<" "<<closestPoint1.y()<<" "<<closestPoint1.z()<<std::endl;
//         if(diss[0] <= diss[1]){
//             if(dis < diss[0]){
//                 points[0] = closestPoint1;
//                 diss[0] = dis;
//             }
//             else if(dis > diss[0] && dis < diss[1]){
//                 points[1] = closestPoint1;
//                 diss[1] = dis;
//             }
//         }
//         else{
//             if(dis < diss[1]){
//                 points[1] = closestPoint1;
//                 diss[1] = dis;
//             }
//             else if(dis > diss[1] && dis < diss[0]){
//                 points[0] = closestPoint1;
//                 diss[0] = dis;
//             }
//         }
//     }

//     for(int i=0;i<points.size();i++){
//         std::cout<<diss[i]<<" ("<<points[i].x()<<" "<<points[i].y()<<" "<<points[i].z()<<")"<<std::endl;
//     }

// }

// TEST_CASE("tiger 1.5 project", "test_project"){
//     // Eigen::RowVector3d e1 = {2.0,0.0,1.0};
//     // Eigen::RowVector3d e2 = {-2.0,0.5,1.0};
//     Eigen::RowVector3d e1 = {2.492269992828369,-0.05197148397564888,0.4495629668235779};
//     Eigen::RowVector3d e2 = {2.492269992828369,0.9480285048484802,0.4495629668235779};

//     // TiGER::Mesh mesh_edge;  
//     // mesh_edge.Vertex.resize(2,3);
//     // mesh_edge.Vertex(0,0) = 2.0;
//     // mesh_edge.Vertex(0,1) = 0.0;
//     // mesh_edge.Vertex(0,2) = 1.0;
//     // mesh_edge.Vertex(1,0) = -2.0;
//     // mesh_edge.Vertex(1,1) = 0.5;
//     // mesh_edge.Vertex(1,2) = 1.0;
//     // mesh_edge.Topo.resize(1,2);
//     // mesh_edge.Topo(0,0) = 0;
//     // mesh_edge.Topo(0,1) = 1;
//     // TiGER::MESHIO::writeVTK("C:/Users/10170/Desktop/segment.vtk",mesh_edge);


//     TiGER::Mesh mesh;    
//     TiGER::MESHIO::readVTK("C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/tris.o.vtk",mesh);
    
//     TiGER::SurfaceMesh surf_in;
//     TiGER::matrix_to_list<3>(mesh.Vertex, surf_in.coords);
//     TiGER::matrix_to_list<3>(mesh.Topo, surf_in.tris);

//     double eps = 0.5;

//     // Eigen::RowVector3d a = {0.0,0.0,0.0};
//     // Eigen::RowVector3d b = {0.0,1.0,0.0};
//     // Eigen::RowVector3d c = {1.0,0.0,0.0};
//     // TiGER::SurfaceMesh surf_in;
//     // std::array<double, 3> pt0 = {0,0,0};
//     // std::array<double, 3> pt1 = {0,1,0};
//     // std::array<double, 3> pt2 = {1,0,0};
//     // std::array<double, 3> pt3 = {-0.5,0.5,-0.5};
//     // std::array<int, 3> tri = {0,1,2};
//     // std::array<int, 3> tri1 = {0,1,3};
//     // surf_in.coords.push_back(pt0);
//     // surf_in.coords.push_back(pt1);
//     // surf_in.coords.push_back(pt2);
//     // surf_in.coords.push_back(pt3);
//     // surf_in.tris.push_back(tri);
//     // surf_in.tris.push_back(tri1);

//     std::vector<Eigen::RowVector3d> pro_points;

//     TiGER::repair::project_SegmentOntoSurface(
//         e1,
//         e2,
//         surf_in,
//         eps,
//         pro_points);

//     // double dis1 = -1; 
//     // double dis2 = -1; 
//     // Eigen::RowVector3d e1_;
//     // Eigen::RowVector3d e2_;
//     // int tri_id1 = -1;
//     // int tri_id2 = -1;

//     // // 线段端点的投影点（最近距离）
//     // for(int i=0;i<surf_in.tris.size();i++){
//     //     // std::cout<<"("<<surf_in.coords[surf_in.tris[i][0]][0]<<","<<surf_in.coords[surf_in.tris[i][0]][1]<<","<<surf_in.coords[surf_in.tris[i][0]][2]<<")"<<std::endl;
//     //     Eigen::RowVector3d a = Eigen::RowVector3d(surf_in.coords[surf_in.tris[i][0]][0],
//     //                             surf_in.coords[surf_in.tris[i][0]][1],
//     //                             surf_in.coords[surf_in.tris[i][0]][2]);
//     //     // std::cout<<"("<<surf_in.coords[surf_in.tris[i][1]][0]<<","<<surf_in.coords[surf_in.tris[i][1]][1]<<","<<surf_in.coords[surf_in.tris[i][1]][2]<<")"<<std::endl;                        
//     //     Eigen::RowVector3d b = Eigen::RowVector3d(surf_in.coords[surf_in.tris[i][1]][0],
//     //                             surf_in.coords[surf_in.tris[i][1]][1],
//     //                             surf_in.coords[surf_in.tris[i][1]][2]);
//     //     // std::cout<<"("<<surf_in.coords[surf_in.tris[i][2]][0]<<","<<surf_in.coords[surf_in.tris[i][2]][1]<<","<<surf_in.coords[surf_in.tris[i][2]][2]<<")"<<std::endl;
//     //     Eigen::RowVector3d c = Eigen::RowVector3d(surf_in.coords[surf_in.tris[i][2]][0],
//     //                             surf_in.coords[surf_in.tris[i][2]][1],
//     //                             surf_in.coords[surf_in.tris[i][2]][2]);

//     //     Eigen::RowVector3d p1;
//     //     double distance1;
//     //     TiGER::repair::project_PointOntoTriangle(
//     //         e1,a,b,c,p1,distance1);
        
//     //     if(fabs(distance1) <= eps && (dis1 == -1 || fabs(distance1) < dis1)){
//     //         dis1 = distance1;
//     //         e1_ = p1;
//     //         tri_id1 = i;
//     //     }

//     //     Eigen::RowVector3d p2;
//     //     double distance2;

//     //     TiGER::repair::project_PointOntoTriangle(
//     //         e2,a,b,c,p2,distance2);
//     //     if(fabs(distance2) <= eps && (dis2 == -1 || fabs(distance2) < dis2)){
//     //         dis2 = distance2;
//     //         e2_ = p2;
//     //         tri_id2 = i;
//     //     }

//     //     // 投影点在三角形内部
//     //     // // 端点投影
//     //     // Eigen::RowVector3d normal = (b - a).cross(c - a).normalized();
//     //     // // e1 点在平面的投影点
//     //     // double distance1 = (e1 - a).dot(normal);
//     //     // Eigen::RowVector3d p1 = e1 - distance1 * normal;
        
//     //     // if(TiGER::repair::isPointInTriangle(p1,a,b,c) && (dis1 == -1 || distance1 < dis1) && fabs(distance1) <= eps){
//     //     //     // std::cout<<(e1 - p1).norm()<<std::endl;
//     //     //     // std::cout<<distance1<<" "<<p1.x()<<" "<<p1.y()<<" "<<p1.z()<<" -- 1 "<<i<<std::endl;
//     //     //     // std::cout<<e1.x()<<" "<<e1.y()<<" "<<e2.z()<<std::endl;
//     //     //     dis1 = distance1;
//     //     //     e1_ = p1;
//     //     //     tri_id1 = i;
//     //     // }

//     //     // double distance2 = (e2 - a).dot(normal);
//     //     // Eigen::RowVector3d p2 = e2 - distance2 * normal;
        
//     //     // if(TiGER::repair::isPointInTriangle(p2,a,b,c) && (dis2 == -1 || fabs(distance2) < dis2)){
//     //     //     // std::cout<<distance2<<" "<<p2.x()<<" "<<p2.y()<<" "<<p2.z()<<" -- 2 "<<i<<std::endl;
//     //     //     dis2 = distance2;
//     //     //     e2_ = p2;
//     //     //     tri_id2 = i;
//     //     // }
//     // }

//     // if(dis1 != -1 && dis1 <= eps){
//     //     pro_points.push_back(e1_);
//     // }
//     // // if(dis2 != -1 && dis2 <= eps){
//     // //     pro_points.push_back(e2_);
//     // // }

//     // // if(tri_id1 == tri_id2){
//     // //     if(dis2 != -1 && dis2 <= eps){
//     // //         pro_points.push_back(e2_);
//     // //     }
//     // //     // TiGER::Mesh mesh_tri;
//     // //     // mesh_tri.Vertex.resize(pro_points.size(),3);
//     // //     // for(int i=0;i<pro_points.size();i++){
//     // //     //     mesh_tri.Vertex(i,0) = pro_points[i].x();
//     // //     //     mesh_tri.Vertex(i,1) = pro_points[i].y();
//     // //     //     mesh_tri.Vertex(i,2) = pro_points[i].z();
//     // //     // }
//     // //     // TiGER::MESHIO::writeVTK("C:/Users/10170/Desktop/print.vtk",mesh_tri);
//     // // }

//     // // 构建三角形邻接关系
//     // std::map<std::array<int,2>,std::vector<int>> edge_tris;  // 边的邻接三角形
//     // for(int i=0;i<surf_in.tris.size();i++){
//     //     int p1 = surf_in.tris[i][0];
//     //     int p2 = surf_in.tris[i][1];
//     //     int p3 = surf_in.tris[i][2];
//     //     std::array<int,2> edge1 = {p1,p2};
//     //     std::array<int,2> edge2 = {p1,p3};
//     //     std::array<int,2> edge3 = {p2,p3};
//     //     std::sort(std::begin(edge1), std::end(edge1));
//     //     std::sort(std::begin(edge2), std::end(edge2));
//     //     std::sort(std::begin(edge3), std::end(edge3));
//     //     edge_tris[edge1].push_back(i);
//     //     edge_tris[edge2].push_back(i);
//     //     edge_tris[edge3].push_back(i);
//     // }

//     // // 找出两段之间的三角形
//     // std::queue<int> q;
//     // q.push(tri_id1);
//     // while(!q.empty()){
//     //     int index = q.front();
//     //     q.pop();

//     //     if(index == tri_id2)
//     //         break;

//     //     int p1 = surf_in.tris[index][0];
//     //     int p2 = surf_in.tris[index][1];
//     //     int p3 = surf_in.tris[index][2];
//     //     Eigen::RowVector3d x = Eigen::RowVector3d(
//     //             (surf_in.coords[p1][0] + surf_in.coords[p2][0] + surf_in.coords[p3][0]) / 3,
//     //             (surf_in.coords[p1][1] + surf_in.coords[p2][1] + surf_in.coords[p3][1]) / 3,
//     //             (surf_in.coords[p1][2] + surf_in.coords[p2][2] + surf_in.coords[p3][2]) / 3);
//     //     std::array<int,2> edge1 = {p1,p2};
//     //     std::array<int,2> edge2 = {p1,p3};
//     //     std::array<int,2> edge3 = {p2,p3};
//     //     std::sort(std::begin(edge1), std::end(edge1));
//     //     std::sort(std::begin(edge2), std::end(edge2));
//     //     std::sort(std::begin(edge3), std::end(edge3));
//     //     double va = -1;
//     //     int next_id = -1;
//     //     int edge_id = -1;
//     //     // 当前三角形与邻接三角形重心方向 点乘 两投影点方向 e1_e2_ 为正
//     //     for(int i=0;i<edge_tris[edge1].size();i++){
//     //         if(edge_tris[edge1][i] != index){
//     //             int id = edge_tris[edge1][i];
//     //             int q1 = surf_in.tris[id][0];
//     //             int q2 = surf_in.tris[id][1];
//     //             int q3 = surf_in.tris[id][2];
//     //             Eigen::RowVector3d y = Eigen::RowVector3d(
//     //                 (surf_in.coords[q1][0] + surf_in.coords[q2][0] + surf_in.coords[q3][0]) / 3,
//     //                 (surf_in.coords[q1][1] + surf_in.coords[q2][1] + surf_in.coords[q3][1]) / 3,
//     //                 (surf_in.coords[q1][2] + surf_in.coords[q2][2] + surf_in.coords[q3][2]) / 3);
//     //             double d = (y - x).normalized().dot(e2_ - e1_);
//     //             // std::cout<<index<<" is "<<id<<" "<<d<<" 1"<<std::endl;
//     //             if(d > va){
//     //                 va = d;
//     //                 next_id = id;
//     //                 edge_id = 1;
//     //             }
//     //         }
//     //     }
//     //     for(int i=0;i<edge_tris[edge2].size();i++){
//     //         if(edge_tris[edge2][i] != index){
//     //             int id = edge_tris[edge2][i];
//     //             int q1 = surf_in.tris[id][0];
//     //             int q2 = surf_in.tris[id][1];
//     //             int q3 = surf_in.tris[id][2];
//     //             Eigen::RowVector3d y = Eigen::RowVector3d(
//     //                 (surf_in.coords[q1][0] + surf_in.coords[q2][0] + surf_in.coords[q3][0]) / 3,
//     //                 (surf_in.coords[q1][1] + surf_in.coords[q2][1] + surf_in.coords[q3][1]) / 3,
//     //                 (surf_in.coords[q1][2] + surf_in.coords[q2][2] + surf_in.coords[q3][2]) / 3);
//     //             double d = (y - x).normalized().dot(e2_ - e1_);
//     //             // std::cout<<index<<" is "<<id<<" "<<d<<" 2"<<std::endl;
//     //             if(d > va){
//     //                 va = d;
//     //                 next_id = id;
//     //                 edge_id = 2;
//     //             }
//     //         }
//     //     }
//     //     for(int i=0;i<edge_tris[edge3].size();i++){
//     //         if(edge_tris[edge3][i] != index){
//     //             int id = edge_tris[edge3][i];
//     //             int q1 = surf_in.tris[id][0];
//     //             int q2 = surf_in.tris[id][1];
//     //             int q3 = surf_in.tris[id][2];
//     //             Eigen::RowVector3d y = Eigen::RowVector3d(
//     //                 (surf_in.coords[q1][0] + surf_in.coords[q2][0] + surf_in.coords[q3][0]) / 3,
//     //                 (surf_in.coords[q1][1] + surf_in.coords[q2][1] + surf_in.coords[q3][1]) / 3,
//     //                 (surf_in.coords[q1][2] + surf_in.coords[q2][2] + surf_in.coords[q3][2]) / 3);
//     //             double d = (y - x).normalized().dot(e2_ - e1_);
//     //             // std::cout<<index<<" is "<<id<<" "<<d<<" 3"<<std::endl;
//     //             if(d > va){
//     //                 va = d;
//     //                 next_id = id;
//     //                 edge_id = 3;
//     //             }
//     //         }
//     //     }
//     //     // std::cout<<index<<" --> "<<next_id<<" "<<edge_id<<std::endl;
//     //     q.push(next_id);

//     //     Eigen::RowVector3d a;
//     //     Eigen::RowVector3d b;
//     //     if(edge_id == 1){
//     //         // std::cout<<p1<<" ?? "<<p2<<std::endl;
//     //         a = Eigen::RowVector3d(
//     //                 surf_in.coords[p1][0],
//     //                 surf_in.coords[p1][1],
//     //                 surf_in.coords[p1][2]);
//     //         b = Eigen::RowVector3d(
//     //                 surf_in.coords[p2][0],
//     //                 surf_in.coords[p2][1],
//     //                 surf_in.coords[p2][2]); 
//     //     }
//     //     else if(edge_id == 3){
//     //         // std::cout<<p2<<" ?? "<<p3<<std::endl;
//     //         a = Eigen::RowVector3d(
//     //                 surf_in.coords[p2][0],
//     //                 surf_in.coords[p2][1],
//     //                 surf_in.coords[p2][2]);
//     //         b = Eigen::RowVector3d(
//     //                 surf_in.coords[p3][0],
//     //                 surf_in.coords[p3][1],
//     //                 surf_in.coords[p3][2]); 
//     //     }
//     //     else if(edge_id == 2){
//     //         // std::cout<<p3<<" ?? "<<p1<<std::endl;
//     //         a = Eigen::RowVector3d(
//     //                 surf_in.coords[p3][0],
//     //                 surf_in.coords[p3][1],
//     //                 surf_in.coords[p3][2]);
//     //         b = Eigen::RowVector3d(
//     //                 surf_in.coords[p1][0],
//     //                 surf_in.coords[p1][1],
//     //                 surf_in.coords[p1][2]); 
//     //     }

//     //     Eigen::RowVector3d closestPoint1;
//     //     Eigen::RowVector3d closestPoint2;
//     //     double dis = TiGER::repair::project_SegmentOntoSegment(e1,e2,a,b,closestPoint1,closestPoint2);
//     //     pro_points.push_back(closestPoint2);
//     // }
//     // pro_points.push_back(e2_);


//     // for(int i=0;i<surf_in.tris.size();i++){
//     //     // std::cout<<"("<<surf_in.coords[surf_in.tris[i][0]][0]<<","<<surf_in.coords[surf_in.tris[i][0]][1]<<","<<surf_in.coords[surf_in.tris[i][0]][2]<<")"<<std::endl;
//     //     Eigen::RowVector3d a = Eigen::RowVector3d(surf_in.coords[surf_in.tris[i][0]][0],
//     //                             surf_in.coords[surf_in.tris[i][0]][1],
//     //                             surf_in.coords[surf_in.tris[i][0]][2]);
//     //     // std::cout<<"("<<surf_in.coords[surf_in.tris[i][1]][0]<<","<<surf_in.coords[surf_in.tris[i][1]][1]<<","<<surf_in.coords[surf_in.tris[i][1]][2]<<")"<<std::endl;                        
//     //     Eigen::RowVector3d b = Eigen::RowVector3d(surf_in.coords[surf_in.tris[i][1]][0],
//     //                             surf_in.coords[surf_in.tris[i][1]][1],
//     //                             surf_in.coords[surf_in.tris[i][1]][2]);
//     //     // std::cout<<"("<<surf_in.coords[surf_in.tris[i][2]][0]<<","<<surf_in.coords[surf_in.tris[i][2]][1]<<","<<surf_in.coords[surf_in.tris[i][2]][2]<<")"<<std::endl;
//     //     Eigen::RowVector3d c = Eigen::RowVector3d(surf_in.coords[surf_in.tris[i][2]][0],
//     //                             surf_in.coords[surf_in.tris[i][2]][1],
//     //                             surf_in.coords[surf_in.tris[i][2]][2]);

//     //     std::vector<Eigen::RowVector3d> points;
//     //     std::vector<double> diss;

//     //     Eigen::RowVector3d closestPoint1;
//     //     Eigen::RowVector3d closestPoint2;
//     //     double dis = TiGER::repair::project_SegmentOntoSegment(e1,e2,a,b,closestPoint1,closestPoint2);
//     //     // diss.push_back(dis);
//     //     if(closestPoint1 != e1 && closestPoint1 != e2 && dis < eps){
//     //         // pro_points.push_back(closestPoint2);
//     //         if(i == 26 || i == 27)
//     //             std::cout<<i<<" ab "<<closestPoint1.x()<<" "<<closestPoint1.y()<<" "<<closestPoint1.z()<<" -- "<<closestPoint2.x()<<" "<<closestPoint2.y()<<" "<<closestPoint2.z()<<std::endl;
//     //         points.push_back(closestPoint2);
//     //         diss.push_back(dis);
//     //     }
//     //     dis = TiGER::repair::project_SegmentOntoSegment(e1,e2,b,c,closestPoint1,closestPoint2);
//     //     // diss.push_back(dis);
//     //     if(closestPoint1 != e1 && closestPoint1 != e2 && dis <= eps){
//     //         // pro_points.push_back(closestPoint2);
//     //         if(i == 26 || i == 27)
//     //             std::cout<<i<<" bc "<<closestPoint1.x()<<" "<<closestPoint1.y()<<" "<<closestPoint1.z()<<" -- "<<closestPoint2.x()<<" "<<closestPoint2.y()<<" "<<closestPoint2.z()<<std::endl;
//     //         points.push_back(closestPoint2);
//     //         diss.push_back(dis);
//     //     }
//     //     dis = TiGER::repair::project_SegmentOntoSegment(e1,e2,c,a,closestPoint1,closestPoint2);
//     //     if(closestPoint1 != e1 && closestPoint1 != e2 && dis <= eps){
//     //         if(i == 26 || i == 27)
//     //             std::cout<<i<<" ca "<<closestPoint1.x()<<" "<<closestPoint1.y()<<" "<<closestPoint1.z()<<" -- "<<closestPoint2.x()<<" "<<closestPoint2.y()<<" "<<closestPoint2.z()<<std::endl;
//     //         if(points.size() < 2){
//     //             points.push_back(closestPoint2);
//     //             diss.push_back(dis);
//     //         }
//     //         else{
//     //             // pro_points.push_back(closestPoint2);
//     //             // std::cout<<i<<" ca "<<closestPoint2.x()<<" "<<closestPoint2.y()<<" "<<closestPoint2.z()<<std::endl;
//     //             if(diss[0] <= diss[1]){
//     //                 if(dis < diss[0]){
//     //                     points[0] = closestPoint2;
//     //                     diss[0] = dis;
//     //                 }
//     //                 else if(dis > diss[0] && dis < diss[1]){
//     //                     points[1] = closestPoint2;
//     //                     diss[1] = dis;
//     //                 }
//     //             }
//     //             else{
//     //                 if(dis < diss[1]){
//     //                     points[1] = closestPoint2;
//     //                     diss[1] = dis;
//     //                 }
//     //                 else if(dis > diss[1] && dis < diss[0]){
//     //                     points[0] = closestPoint2;
//     //                     diss[0] = dis;
//     //                 }
//     //             }
//     //         }
//     //     }  
//     //     // std::cout<<i<<" "<<points[0].x()<<" "<<points[0].y()<<" "<<points[0].z()<<std::endl;
//     //     // std::cout<<i<<" "<<points[1].x()<<" "<<points[1].y()<<" "<<points[1].z()<<std::endl;

//     //     if(i == tri_id1 || i == tri_id2){
//     //         if(points.size() == 2){
//     //             if(diss[0] <= diss[1]){
//     //                 pro_points.push_back(points[0]);
//     //             }
//     //             else{
//     //                 pro_points.push_back(points[1]);
//     //             }
//     //         }
//     //         else{
//     //             for(int i=0;i<points.size();i++){
//     //                 pro_points.push_back(points[i]);
//     //             }
//     //         }
//     //     }
//     //     else{
//     //         for(int i=0;i<points.size();i++){
//     //             pro_points.push_back(points[i]);
//     //         }
//     //         // pro_points.push_back(points[0]);
//     //         // pro_points.push_back(points[1]);
//     //     }
//     // }

//     TiGER::Mesh mesh_tri;
//     // TiGER::list_to_matrix<3>(result.coords,mesh_tri.Vertex);
//     // TiGER::list_to_matrix<3>(result.tris,mesh_tri.Topo);
//     mesh_tri.Vertex.resize(pro_points.size(),3);
//     for(int i=0;i<pro_points.size();i++){
//         mesh_tri.Vertex(i,0) = pro_points[i].x();
//         mesh_tri.Vertex(i,1) = pro_points[i].y();
//         mesh_tri.Vertex(i,2) = pro_points[i].z();
//     }
//     TiGER::MESHIO::writeVTK("C:/Users/10170/Desktop/print.vtk",mesh_tri);
// }


TEST_CASE("tiger 1.5 project_topu", "test_project_topu"){
    TiGER::Mesh mesh;    
    TiGER::MESHIO::readVTK("C:/Users/10170/Desktop/gridteam/code/FOLLOW/ti-ger-1.5/file/all.o.vtk",mesh);
    
    TiGER::SurfaceMesh surf_in;
    TiGER::matrix_to_list<3>(mesh.Vertex, surf_in.coords);
    TiGER::matrix_to_list<3>(mesh.Topo, surf_in.tris);

    surf_in.v_conn_points.resize(surf_in.coords.size());
    surf_in.v_conn_tris.resize(surf_in.coords.size());
    for(int i=0;i<surf_in.tris.size();i++){
        int p1 = surf_in.tris[i][0];
        int p2 = surf_in.tris[i][1];
        int p3 = surf_in.tris[i][2];
        surf_in.v_conn_tris[p1].insert(i);
        surf_in.v_conn_tris[p2].insert(i);
        surf_in.v_conn_tris[p3].insert(i);
        surf_in.v_conn_points[p1].insert(p2);
        surf_in.v_conn_points[p1].insert(p3);
        surf_in.v_conn_points[p2].insert(p1);
        surf_in.v_conn_points[p2].insert(p3);
        surf_in.v_conn_points[p3].insert(p1);
        surf_in.v_conn_points[p3].insert(p2);
    }
    std::vector<TiGER::DiscreteGeometry::Edge> edges_in;
    edges_in.push_back({5,7});
    edges_in.push_back({4,5});
    edges_in.push_back({4,6});
    edges_in.push_back({6,7});
    std::vector<int> tris_id = {10,12,14,36,38,39,41,42,44,46,47,48,50,51,53,54,56,57,59,61,62,63,65,67,68,70,71,73,74,76,77,78,79,81,82,141,144};
    TiGER::RepairParameters args;
    args.setprint_eps(0.5);
    TiGER::SurfaceMesh surf_out;
    TiGER::discrete_geometry_repair::project(
        surf_in,
        edges_in,
        tris_id,
        args,
        surf_out);
        
    TiGER::Mesh mesh_tri;
    TiGER::list_to_matrix<3>(surf_out.coords,mesh_tri.Vertex);
    TiGER::list_to_matrix<3>(surf_out.tris,mesh_tri.Topo);
    TiGER::MESHIO::writeVTK("C:/Users/10170/Desktop/print.vtk",mesh_tri);
}