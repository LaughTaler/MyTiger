#pragma once

#include <igl/readSTL.h>
#include "CLI11.hpp"
#include <array>
#include <chrono> 
#include <iostream>
#include <memory>
#include <fstream>
#include <cmath>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include "meshIO.h"
#include "list_to_matrix.hpp"
#include "matrix_to_list.hpp"
#include "geometry_interface_plane.h"
#include "tiger.h"
#include "tiger_parameters.h"
#include "igl/bfs_orient.h"
#include "facet_classfication.hpp"
#include "merge_volume_mesh.hpp"
#include "surface_from_volume.hpp"
#include "cgnsWrite.h"
#include "add_source.hpp"
#include "boundary_layer_mesh.h"
#include "spdlog/spdlog.h"
#include "alias.h"
#include "MeshOrient.h"
#include "kdtree.h"

class KdTreeM {

public:

    kdtree* kd;
    kdres* res;

    KdTreeM::KdTreeM(int dimension)
    {
        kd = kd_create(dimension);
    }

    KdTreeM::~KdTreeM()
    {
        kd_free(kd);
    }

    bool KdTreeM::Insert3DNode(Eigen::Vector3d& point, int id)
    {
        int* id_ptr = new int(1);
        *id_ptr = id;
        kd_insert3(this->kd, point.x(), point.y(), point.z(), id_ptr);
        return true;
    }

    std::vector<int> KdTreeM::Query3DNodeByDistance(Eigen::Vector3d& point, double dis)
    {
        res = kd_nearest_range3(kd, point.x(), point.y(), point.z(), dis);
        std::vector<int> vec;
        for (int i = 0; i < kd_res_size(res); i++)
        {
            int* id = static_cast<int*>(kd_res_item_data(res));
            vec.push_back(id[0]);
            if (kd_res_next(res) == 0) break;
        }
        if (res != nullptr) kd_res_free(res);
        return vec;
    }
};


void dfs(int node, std::unordered_map<int, std::vector<int>> adj_list, std::set<int>& visited,
                        std::vector<int>& order)
{
    visited.insert(node);
    order.push_back(node);
    for( int neighbor : adj_list[node] )
    {
        if( visited.find(neighbor) == visited.end() )
        {
            dfs(neighbor, adj_list, visited, order);
        }
    }
}

std::vector<int> sort_point_order(const std::vector<std::pair<int, int>>& segments, int start_node)
{
    std::unordered_map<int, std::vector<int>> adj_list;
    for( const auto& segment : segments )
    {
        adj_list[segment.first].push_back(segment.second);
        adj_list[segment.second].push_back(segment.first);
    }
    std::set<int> visited;
    std::vector<int> order;
    dfs(start_node, adj_list, visited, order);
    return order;
}

void ordered_point(
    const std::vector<int>& disordered_points_id, 
    const std::vector<std::pair<int, int>>& mesh_edge_lst,
    std::unordered_map<int, bool> point_exist,
    int start_point,
    std::vector<int>& ordered_points_id)
{
    std::vector<std::pair<int, int>> temp_edge_lst;
    for( int i = 0; i < mesh_edge_lst.size(); i++ )
    { 
        if( point_exist[mesh_edge_lst[i].first] && point_exist[mesh_edge_lst[i].second] )
        {
            temp_edge_lst.push_back(mesh_edge_lst[i]);
        }
    }
    ordered_points_id = sort_point_order(temp_edge_lst, start_point);
}

static void points_map(
    const std::vector<std::array<double, 3>>& tri_mesh_points,
    const std::vector<std::array<double, 3>>& prism_mesh_points,
    std::vector<int>& ordered_points_id)
{
    std::vector<Eigen::Vector3d> line_mesh(tri_mesh_points.size());
    std::vector<Eigen::Vector3d> bl_mesh(prism_mesh_points.size());
    for( int i = 0; i < tri_mesh_points.size(); i++ )
    {
        Eigen::Vector3d temp(tri_mesh_points[i][0], tri_mesh_points[i][1], tri_mesh_points[i][2]);
        line_mesh[i] = temp;
    }
    for( int i = 0; i < prism_mesh_points.size(); i++ )
    {
        Eigen::Vector3d temp(prism_mesh_points[i][0], prism_mesh_points[i][1], prism_mesh_points[i][2]);
        bl_mesh[i] = temp;
    }
    KdTreeM tree(3);
    for( int i = 0; i < bl_mesh.size(); i++ )
    {
        Eigen::Vector3d vec = bl_mesh[i];
        tree.Insert3DNode(vec, i);
    }
    for( int i = 0; i < line_mesh.size(); i++ )
    {
        Eigen::Vector3d vec = line_mesh[i];
        double eps = 1e-25;
        auto ans = tree.Query3DNodeByDistance(vec, eps);
        while (ans.empty()) {
            eps *= 2;
            ans = tree.Query3DNodeByDistance(vec, eps);
        }
        ordered_points_id[i] = ans[0];
    }
}

void findTopLayerPoints(
    std::vector<int>& ordered_points,
    unordered_map<int, int> bl_point_topo,
    std::set<std::pair<int, int>> pyramid_edges)
{
    std::vector<int> ans;
    std::vector<vector<int>> layer(ordered_points.size());
    for( int i = 0; i < ordered_points.size(); i++ )
    {
        int temp_point = ordered_points[i];
        layer[i].push_back(temp_point);
        while( bl_point_topo.find(temp_point) != bl_point_topo.end() )
        {
            temp_point = bl_point_topo[temp_point];
            layer[i].push_back(temp_point);
        }
       // ordered_points[i] = temp_point;
    }
    
    for (int i = 0; i < ordered_points.size(); i++)
    {
        auto temp_point = layer[i].back();
        if (!i) {
            ans.push_back(temp_point);
            continue;
        }
        if (layer[i].size() > layer[i - 1].size()) {
            if (pyramid_edges.find(std::make_pair(layer[i].back(), layer[i - 1].back())) == pyramid_edges.end()) {
                for (int j = layer[i - 1].size()-1; j < layer[i].size(); j++) {
                    ans.push_back(layer[i][j]);
                }
            }
            else {
                ans.push_back(temp_point);
            }
        }
        else if (layer[i].size() == layer[i - 1].size()) {
            ans.push_back(temp_point);
        }
        else if (layer[i].size() < layer[i - 1].size()) {
            if (pyramid_edges.find(std::make_pair(layer[i].back(), layer[i - 1].back())) == pyramid_edges.end()) {
                for (int j = int(layer[i-1].size() - 2); j >= int(layer[i].size()-1); j--) {
                    ans.push_back(layer[i - 1][j]);
                }
                ans.push_back(temp_point);
            }
            else {
                ans.push_back(temp_point);
            }
        }

    }
    ordered_points = ans;

}

double vector_differ(
    std::array<double, 3> vector_a,
    std::array<double, 3> vector_b)
{
    Eigen::Vector3d va(vector_a[0], vector_a[1], vector_a[2]);
    Eigen::Vector3d vb(vector_b[0], vector_b[1], vector_b[2]);
    double ans = (va - vb).norm();
    return ans;
}

std::array<double, 3> mat_multiply_vec(
    Eigen::Matrix3d matrix_,
    std::array<double, 3> vec)
{
    Eigen::Vector3d temp(vec[0], vec[1], vec[2]);
    temp = matrix_ * temp;
    std::array<double, 3> ans{ temp.x(),temp.y(),temp.z() };
    return ans;
}

struct Args {
    std::string filepath = "D:\\code\\turbomachinery\\tiger1.5\\example\\";
    std::string geometry_name="3_0.125_cf";
    double maxsize = 10;
    std::vector<int> symmetry_surface{1,3};
    std::vector<int> wall_surface;
    std::vector<int> far_surface{ 0,2,4 };
    //double angle = 10;
    //Eigen::Vector3d axis_direction{ 1, 0, 0 };
    //std::vector<int> source_surface;
    //std::vector<int> target_surface;
    int number_of_layer = 25;
    double step_length_layer = 0.001;
    double ratio_layer = 1.2;
};

void hybrid_mesh_generation(
    Args args,
    TiGER::Geometry geometry,
    TiGER::SizingFunction size_field,
    std::vector<TiGER::SurfaceMesh>& local_surface_mesh,
    std::vector<TiGER::LineMesh>& linemesh)
{

}

int main(int argc, char* argv[]) {
    try {
        Args args;
        // CLI::App app("test_tiger"); 
        // app.add_option("-i", args.filepath, "input file name");
        // app.add_option("-s", args.sourcename, "input source file name");
        // app.add_option("--size", args.maxsize, "global size");
        // app.add_option("-n", args.number_of_layer, "number of boundary layer mesh");
        // app.add_option("-l", args.step_length_layer, "step length of boundary layer mesh");
        // app.add_option("-r", args.ratio_layer, "growth ratio of boundary layer mesh");
        // CLI11_PARSE(app, argc, argv);

        /*------------------------------read geometry-----------------------------*/
        args.filepath = "";
        std::string input_file_name_geometry = args.filepath + args.geometry_name + ".stp";
        TiGER::Geometry geometry;
        TiGER::readModelParameters args_geometry;
        spdlog::info("the input geometry file is:{}", input_file_name_geometry);
        TiGER::Read_Model(input_file_name_geometry, geometry, args_geometry);

        TiGER::generateStlParameters args_stl;
        TiGER::SurfaceMesh stl;
        args_stl.setAngle(1);
        args_stl.setDeflection(0.001);
        TiGER::CAD_Tessellation(geometry, stl, args_stl);

        TiGER::Mesh m;
        TiGER::list_to_matrix<3>(stl.coords, m.Vertex);
        TiGER::list_to_matrix<3>(stl.tris, m.Topo);
        TiGER::MESHIO::writeVTK(args.filepath + args.geometry_name + "_stl.vtk", m);

        bool hashTable[1000] = { 0 };
        for (int i = 0; i < args.symmetry_surface.size(); i++) {
            hashTable[args.symmetry_surface[i]] = 1;
        }
        for (int i = 0; i < args.far_surface.size(); i++) {
            hashTable[args.far_surface[i]] = 1;
        }
        for (int i = 0; i < geometry.surfaces_.size(); i++) {
            if (!hashTable[i]) {
                args.wall_surface.push_back(i);
            }
        }
        /*------------------------------generate size function-----------------------------*/
        TiGER::BackGroundParameters args_sizefiled;
        TiGER::SizingManager sizingmanager;
        add_source(args.filepath + "source.txt", sizingmanager);
        TiGER::SizingFunction size_field = [&sizingmanager, &args](const double& x, const double& y, const double& z) {
            //double size = TiGER::size_field::SizingFunction_getSizeAtPoint(x, y, z, sizingmanager);
            //return std::fmin(size, args.maxsize);
            if (y > 50 || z<-550||z>60) {
                return 30;
            }
            else if (y > 40 || z < -500 || z>50) {
                return 10;
            }
            else {
                return 3;
            }

        };

        /*------------------------------calculate eps------------------------*/
        double min_length = (std::numeric_limits<double>::max)();
        for (int i = 0; i < geometry.curves_.size(); i++) {
            if (geometry.curves_[i]->getTotalLenthFunction() > 0 && geometry.curves_[i]->getTotalLenthFunction() < min_length) {
                min_length = geometry.curves_[i]->getTotalLenthFunction();
            }
        }
        double eps = min_length/10;

        /*------------------------------transforamtion matrix-----------------------------*/
        //double radians = args.angle * M_PI / 180.0;
        //Eigen::AngleAxisd r(radians, args.axis_direction);
        //Eigen::Matrix3d rotation = r.toRotationMatrix();

        /*------------------------------find periodic relationship------------------------*/
        //std::unordered_map<int, int> periodic_surface;
        //for( int i = 0; i < args.source_surface.size(); i++ )
        //{
        //    periodic_surface[args.source_surface[i]] = args.target_surface[i];
        //    periodic_surface[args.target_surface[i]] = args.source_surface[i];
        //}

        //std::unordered_map<int, int> periodic_curve;
        //for (int i = 0; i < args.source_surface.size(); i++) {
        //    int temp_source_surface = args.source_surface[i];
        //    for (int j = 0; j < geometry.topo_surf_to_curves_[temp_source_surface].size(); j++) {
        //        int source_curve_id = geometry.topo_surf_to_curves_[temp_source_surface][j];
        //        std::array<double, 3> first_source_point = geometry.curves_[source_curve_id]->d0Function(0);
        //        std::array<double, 3> second_source_point = geometry.curves_[source_curve_id]->d0Function(1);              //需要获取边界线终点坐标??????
        //        std::array<double, 3> first_compute_point = mat_multiply_vec(rotation, first_source_point);
        //        std::array<double, 3> second_compute_point = mat_multiply_vec(rotation, second_source_point);
        //        for (int k = 0; k < geometry.topo_surf_to_curves_[args.target_surface[i]].size(); k++) {
        //            int cur_curve_id = geometry.topo_surf_to_curves_[args.target_surface[i]][k];
        //            std::array<double, 3> first_cur_point = geometry.curves_[cur_curve_id]->d0Function(0);
        //            std::array<double, 3> second_cur_point = geometry.curves_[cur_curve_id]->d0Function(1);       //需要获取边界线终点坐标????????????

        //            if (vector_differ(first_cur_point, first_compute_point) < eps && vector_differ(second_cur_point, second_compute_point) < eps) {
        //                periodic_curve[cur_curve_id] = source_curve_id;
        //                periodic_curve[source_curve_id] = cur_curve_id;
        //                break;
        //            }
        //            else if (vector_differ(first_cur_point, second_compute_point) < eps  && vector_differ(second_cur_point, first_compute_point) < eps) {
        //                periodic_curve[cur_curve_id] = source_curve_id;
        //                periodic_curve[source_curve_id] = cur_curve_id;
        //                break;
        //            }
        //        }
        //    }
        //}

        /*------------------------------generate surface mesh-----------------------------*/
        TiGER::SurfaceMesh global_surface_mesh;   //global surface mesh
        std::vector<TiGER::SurfaceMesh> local_surface_mesh(geometry.surfaces_.size());  //local surface mesh
        std::vector<TiGER::LineMesh> linemesh(geometry.curves_.size());
        TiGER::CurveParameters args_linemesh;
        TiGER::IsotropicSurfaceParametersTri args_surfacemesh;
        std::unordered_map<int, bool> finished_curve_mesh;
        std::unordered_map<int, bool> finished_surface_mesh;

        ////discrete periodic curves
        //for( int i = 0; i < args.source_surface.size(); i++ )
        //{
        //    int temp_surface = args.source_surface[i];
        //    for( int j = 0;j < geometry.topo_surf_to_curves_[temp_surface].size(); j++ )
        //    {
        //        int temp_curve = geometry.topo_surf_to_curves_[temp_surface][j];
        //        int target_curve = periodic_curve[temp_curve];
        //        TiGER::discretize_curve::Grid1D_Legacy(geometry.curves_[temp_curve], size_field, args_linemesh, linemesh[temp_c*urve]);
        //        TiGER::LineMesh target_curve_mesh;        //能否只定义coord 
        //        linemesh[target_curve] = target_curve_mesh;   
        //        finished_curve_mesh[temp_curve] = true;
        //        finished_curve_mesh[target_curve] = true;
        //    }
        //}

        //discrete other curves
        for (int i = 0; i < geometry.curves_.size(); i++)
        {
            if( !finished_curve_mesh[i] )
            {
                TiGER::discretize_curve::Grid1D_Legacy(geometry.curves_[i], size_field, args_linemesh, linemesh[i]);
            }
        }

        ////discrete periodic surfaces
        //for( int i = 0; i < args.source_surface.size(); i++ )
        //{
        //    int temp_surface = args.source_surface[i];
        //    int target_surface = periodic_surface[temp_surface];
        //    std::vector<TiGER::LineMesh> surface_boundary;
        //    for( int j = 0; j < geometry.topo_surf_to_curves_[temp_surface].size(); j++ )
        //    {
        //        surface_boundary.emplace_back(linemesh[geometry.topo_surf_to_curves_[i][j]]);
        //    }
        //    TiGER::discretize_surface::CADSurf_TGrid_AFT(geometry.surfaces_[temp_surface], surface_boundary, args_surfacemesh, size_field, local_surface_mesh[temp_surface]);
        //    TiGER::SurfaceMesh target_surface_mesh;   //能否只定义coord
        //    local_surface_mesh[target_surface] = target_surface_mesh;
        //    finished_surface_mesh[temp_surface] = true;
        //    finished_surface_mesh[target_surface] = true;
        //}
        
        //args_surfacemesh.setmUserSpecifiedMinEdgeLen(2);
        //args_surfacemesh.setmUserSpecifiedMaxEdgeLen(2);
        //args_surfacemesh.setOutputInformation(1);
        //args_surfacemesh.setSavepath(args.filepath);
        
        //discrete other surfaces
        for (int i = 0; i < geometry.surfaces_.size(); i++)
        {
            if( !finished_surface_mesh[i] )
            {
            std::vector<TiGER::LineMesh> surface_boundary;
                for (int j = 0; j < geometry.topo_surf_to_curves_[i].size(); j++)
                {
                    surface_boundary.emplace_back(linemesh[geometry.topo_surf_to_curves_[i][j]]);
                }
                TiGER::discretize_surface::CADSurf_TGrid_AFT(geometry.surfaces_[i], surface_boundary, args_surfacemesh, size_field, local_surface_mesh[i]);
                //TiGER::list_to_matrix<3>(local_surface_mesh[i].coords, m.Vertex);
                //TiGER::list_to_matrix<3>(local_surface_mesh[i].tris, m.Topo);
                //std::string tempname = args.filepath+"face"+std::to_string(i)+".vtk";
                //TiGER::MESHIO::writeVTK(tempname, m);
            }
        }
        TiGER::merge_surface_mesh(local_surface_mesh, global_surface_mesh);
      //  TiGER::MESHIO::writeVTK("global_surface_mesh.vtk", global_surface_mesh);

        // /*------------------------------generate boundary layer-----------------------------*/
        TiGER::BoundaryLayerParameters args_bdry;
        args_bdry.setstep_length(args.step_length_layer);
        args_bdry.setnumber_of_layer(args.number_of_layer);
        args_bdry.setratio(args.ratio_layer);
        args_bdry.setmultiple_normal(false);
        args_bdry.setaniso_iso_blend(0);
        args_bdry.setoutputio(true);
        TiGER::TetrahedraParameters args_vmesh;
        std::vector<std::shared_ptr<TiGER::GeometrySurface>> surfaces;
        TiGER::VolumeMesh vol_mesh;
        TiGER::SurfaceMesh out_boundary_mesh;
        TiGER::SurfaceMesh tris_boundary_mesh;
        std::vector<int> tris_boundary_point_l2g;
        //reset orient
        int mesh_point_num = global_surface_mesh.coords.size();
        int mesh_tri_num = global_surface_mesh.tris.size();
        std::vector<std::vector<double>> point_list(mesh_point_num);
        std::vector<std::vector<int>> facet_list(mesh_tri_num);
        std::vector<int> blockMark;
        for (int i = 0; i < mesh_point_num; i++) {
            point_list[i].push_back(global_surface_mesh.coords[i][0]);
            point_list[i].push_back(global_surface_mesh.coords[i][1]);
            point_list[i].push_back(global_surface_mesh.coords[i][2]);
        }
        for (int i = 0; i < mesh_tri_num; i++) {
            facet_list[i].push_back(global_surface_mesh.tris[i][0]);
            facet_list[i].push_back(global_surface_mesh.tris[i][1]);
            facet_list[i].push_back(global_surface_mesh.tris[i][2]);
        }
        TiGER::resetOrientation(point_list, facet_list, blockMark);
        for (int i = 0; i < facet_list.size(); i++) {
            global_surface_mesh.tris[i][0] = facet_list[i][0];
            global_surface_mesh.tris[i][1] = facet_list[i][2];
            global_surface_mesh.tris[i][2] = facet_list[i][1];
        }

        for (int i = 0; i < geometry.surfaces_.size()+2; i++) {
            global_surface_mesh.bc[i] = TiGER::BoundaryCondition::WALL;
        }
        for (int i = 0; i < args.symmetry_surface.size(); i++) {
            int temp_id = args.symmetry_surface[i];
            global_surface_mesh.bc[temp_id] = TiGER::BoundaryCondition::SYMMETRY;
        }
        for (auto i:args.far_surface) {
            global_surface_mesh.bc[i] = TiGER::BoundaryCondition::FARFIELD;
        }

        TiGER::boundary_layer::HexGrid_createPrismaticMesh(global_surface_mesh, surfaces, args_bdry, vol_mesh, out_boundary_mesh, tris_boundary_mesh, tris_boundary_point_l2g);
        //TiGER::MESHIO::writeVTK("out_boundary_mesh.vtk", out_boundary_mesh);
        //TiGER::MESHIO::writeVTK("tris_boundary_mesh.vtk", tris_boundary_mesh);
        TiGER::MESHIO::writeVTK("vol_mesh.vtk",vol_mesh);
        TiGER::MESHIO::writeVTK("vtris_mesh.vtk", tris_boundary_mesh);

         /*------------------------------regenerate match surface mesh-----------------------------*/
        std::unordered_map<int, int> bl_curve; // 存储曲线到wall面的映射
        for( int i = 0; i < args.wall_surface.size(); i++ )
        {
            int wall_surface = args.wall_surface[i];
            for( int j = 0; j < geometry.topo_surf_to_curves_[wall_surface].size(); j++ )
            {
                bl_curve[geometry.topo_surf_to_curves_[wall_surface][j]] = wall_surface;
            }
        }

        //regenerate periodic curves
        //regenerate match curves
        for( int i = 0; i < args.symmetry_surface.size(); i++ )
        {
            int symm_surface = args.symmetry_surface[i];
            std::vector<std::vector<std::array<double, 3>>> vir_sf_boundary_mesh // 存储所有边界线，1：每条边的线网格 2：网格点组成的线网格 3. 网格点
            (geometry.topo_surf_to_curves_[symm_surface].size());
            for (int j = 0; j < geometry.topo_surf_to_curves_[symm_surface].size(); j++) {
                int temp_curve = geometry.topo_surf_to_curves_[symm_surface][j];
                vir_sf_boundary_mesh[j] = linemesh[temp_curve].coord;
            }

            //regenerate curves related to body surface
            for( int j = 0; j < geometry.topo_surf_to_curves_[symm_surface].size(); j++ )
            {
                int wall_curve = geometry.topo_surf_to_curves_[symm_surface][j];
                if( bl_curve.find(wall_curve) != bl_curve.end() )
                {
                    /// 说明是需要两边重构
                    int wall_surface = bl_curve[wall_curve]; 
                    std::vector<std::array<double, 3>> temp_line=linemesh[wall_curve].coord;
                    //TiGER::MESHIO::writeVTK("LINE", linemesh[temp_curve]);
                    std::vector<std::array<double, 3>> vol = vol_mesh.coords;
                    std::vector<int> ordered_points(temp_line.size());
                    points_map(temp_line, vol_mesh.coords, ordered_points);

                    //find top layer points
                    std::unordered_map<int, int> bl_point_topo;
                    std::set<std::pair<int, int>> pyramid_edge;
                    for (auto i : vol_mesh.pyramids) {
                        for (int j = 0; j < 5; j++) {
                            for (int k = 0; k < 5; k++) {
                                if (j != k)
                                    pyramid_edge.insert(std::make_pair(i[j],i[k]));
                            }
                        }
                    }
                    for( int i = 0; i < vol_mesh.prisms.size(); i++ )
                    {
                          for( int j = 0; j < 3; j++ )
                          {
                              bl_point_topo[vol_mesh.prisms[i][j]] = vol_mesh.prisms[i][j + 3];
                          }
                    }
                    findTopLayerPoints(ordered_points, bl_point_topo, pyramid_edge);

                    //regenerate
                    std::vector<std::array<double, 3>> temp_line_coord(ordered_points.size());
                    vir_sf_boundary_mesh[j].resize(ordered_points.size());
                    for( int k = 0; k < ordered_points.size(); k++ )
                    {
                        temp_line_coord[k]=vol_mesh.coords[ordered_points[k]];
                        vir_sf_boundary_mesh[j][k] = vol_mesh.coords[ordered_points[k]];
                    }
                    TiGER::LineMesh temp_line_mesh;
                    temp_line_mesh.coord = temp_line_coord;   
                    linemesh[wall_curve] = temp_line_mesh;   
                   // TiGER::MESHIO::writeVTK("main_line" + std::to_string(temp_surface) + std::to_string(temp_curve) + ".vtk",linemesh[temp_curve]);

                    //cut side edges
                    /// 只更新一边的网格
                    for (int k = 0; k < geometry.topo_surf_to_curves_[symm_surface].size(); k++) {
                        int cur_curve = geometry.topo_surf_to_curves_[symm_surface][k]; 
                        if (cur_curve != wall_curve&& bl_curve.find(cur_curve) == bl_curve.end()) {
                            double cur_u_max,temp_u_max;
                            double length_max = 1.0;
                            TiGER::mapNormal2Real_Curve(geometry.curves_[cur_curve], cur_u_max, length_max);
                            std::array<double,3> cur_end_coord = geometry.curves_[cur_curve]->d0Function(cur_u_max);
                            std::array<double, 3> cur_start_coord = geometry.curves_[cur_curve]->d0Function(0);
                            TiGER::mapNormal2Real_Curve(geometry.curves_[wall_curve], temp_u_max, length_max);
                            std::array<double, 3> temp_end_coord = geometry.curves_[wall_curve]->d0Function(temp_u_max);
                            std::array<double, 3> temp_start_coord = geometry.curves_[wall_curve]->d0Function(0);

                            if (vector_differ(cur_start_coord,temp_start_coord) < eps)
                            {
                                vir_sf_boundary_mesh[k][0] = vir_sf_boundary_mesh[j][0];
                            }

                            if (vector_differ(cur_start_coord, temp_end_coord) < eps)
                            {
                                vir_sf_boundary_mesh[k][0] = vir_sf_boundary_mesh[j].back();
                            }
                            if (vector_differ(cur_end_coord, temp_start_coord) < eps)
                            {
                                vir_sf_boundary_mesh[k].back() = vir_sf_boundary_mesh[j][0];
                            }
                            if (vector_differ(cur_end_coord, temp_end_coord) < eps) {
                                vir_sf_boundary_mesh[k].back() = vir_sf_boundary_mesh[j].back();
                            }
                        }
                    }
                }
            }

            //regenerate other curves
            for (int j = 0; j < geometry.topo_surf_to_curves_[symm_surface].size(); j++) {
                int temp_curve = geometry.topo_surf_to_curves_[symm_surface][j];
                if (bl_curve.find(temp_curve) == bl_curve.end()) {
                    std::array<double, 3> startPt = vir_sf_boundary_mesh[j][0];
                    std::array<double, 3> endPt = vir_sf_boundary_mesh[j].back();
                    TiGER::LineMesh temp_line_mesh;

                    //temp_line_mesh.coord.resize(2);
                    //temp_line_mesh.coord[0] = startPt;
                    //temp_line_mesh.coord[1] = endPt;
                    

   
                    //std::shared_ptr<TiGER::GeometryCurve> curves =std::make_shared<TiGER::GeometryCurve>(new TiGER::GeometryCurve);
                    //Grid1D_createMeshLine(startPt, endPt, *curves);
                    //TiGER::CurveParameters cp;
                    //TiGER::discretize_curve::Grid1D_Legacy({ curves },size_field,cp, temp_line_mesh);
                    //linemesh[temp_curve] = temp_line_mesh;
                    //TiGER::MESHIO::writeVTK("line2" + std::to_string(temp_surface) + std::to_string(temp_curve) + ".vtk", linemesh[temp_curve]);

                    std::shared_ptr<TiGER::GeometryCurve> curves = std::make_shared<TiGER::GeometryCurve>();

                    Grid1D_createMeshLine(startPt, endPt, *curves);
                    TiGER::CurveParameters cp;
                    
                    TiGER::discretize_curve::Grid1D_Legacy({ curves },size_field,cp, temp_line_mesh);
                    linemesh[temp_curve] = temp_line_mesh;
                    TiGER::MESHIO::writeVTK("vice_line" + std::to_string(symm_surface) + std::to_string(temp_curve) + ".vtk", temp_line_mesh);
                }
            }
        }

        //regenerate match surface
        for (int i = 0; i < args.symmetry_surface.size(); i++) {
            int temp_surface = args.symmetry_surface[i];
            std::vector<TiGER::LineMesh> surface_boundary;
            for (int j = 0; j < geometry.topo_surf_to_curves_[temp_surface].size(); j++)
            {
                surface_boundary.emplace_back(linemesh[geometry.topo_surf_to_curves_[temp_surface][j]]);
            }
            TiGER::discretize_surface::CADSurf_TGrid_AFT(geometry.surfaces_[temp_surface], surface_boundary, args_surfacemesh, size_field, local_surface_mesh[temp_surface]);
            TiGER::MESHIO::writeVTK("renew_surface"+to_string(temp_surface)+".vtk", local_surface_mesh[temp_surface]);

        }

       /*------------------------------generate volume mesh-----------------------------*/
        TiGER::VolumeMesh internal_volume_mesh;
        TiGER::SurfaceMesh new_surface_mesh_global;
        std::vector<TiGER::SurfaceMesh> new_surface_mesh_local;
        for (int i = 0; i < args.symmetry_surface.size(); i++) {
            int temp_surface = args.symmetry_surface[i];
            new_surface_mesh_local.push_back(local_surface_mesh[temp_surface]);
        }
        new_surface_mesh_local.push_back(tris_boundary_mesh);
        merge_surface_mesh(new_surface_mesh_local, new_surface_mesh_global);
        TiGER::MESHIO::writeVTK("final_surface_mesh.vtk", new_surface_mesh_global);
        args_vmesh.setnthread(1);
        args_vmesh.setoptloop(1);
        TiGER::tetrahedral_mesh::TetGrid_Constrained(new_surface_mesh_global, args_vmesh, internal_volume_mesh, nullptr);
        //TiGER::list_to_matrix<3>(volume_mesh.coords, m.Vertex);
        //TiGER::list_to_matrix<4>(volume_mesh.tetras, m.Topo);
        //TiGER::MESHIO::writeVTK(args.filepath + "volumne_mesh_.vtk" , m);

        //output
        TiGER::VolumeMesh volume_global;
        std::vector<TiGER::VolumeMesh> volume_local;
        volume_local.push_back(vol_mesh);
        volume_local.push_back(internal_volume_mesh);
        TiGER::merge_volume_mesh(volume_local,volume_global);

        TiGER::MESHIO::writeVTK(args.filepath+"vmesh_final.vtk",volume_global);
        //TiGER::SurfaceMesh surface_from_volume;
        //TiGER::surfaceFromVolume(volume_global, surface_from_volume);
        //std::cout<<"finish get boundry"<<std::endl;
        ////TiGER::list_to_matrix<3>(surface_from_volume.coords, m.Vertex);
        ////TiGER::list_to_matrix<3>(surface_from_volume.tris, m.Topo);
        ////TiGER::MESHIO::writeVTK(args.filepath+"boundary.vtk", m);
        //std::string cgnsfilename = args.filepath+"vmesh_final.cgns";
        //Tool_writeMeshToCGNS(cgnsfilename.c_str(), surface_from_volume,volume_global);
        //std::cout<<"tiger finish!"<<std::endl;
     }
    catch (std::exception e) {
        std::cout << e.what();
    }
}
