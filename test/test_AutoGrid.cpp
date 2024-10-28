// #define CATCH_CONFIG_MAIN // This tells Catch to provide a main() - only do
// this in one cpp file

#include <chrono>

#include <array>
#include <iostream>
#include <memory>
#include "CLI11.hpp"
#include "autogrid.h"
#include "discretize_surface_remesh.h"
#include "list_to_matrix.hpp"
#include "meshIO.h"
#include "spdlog/spdlog.h"
// #include "tiger.h"
// #include "merge_surface_mesh.hpp"
#include "definitions.h"
#include "alias.h"
#include "../tools/meshIO.h"
// #include "SmeshGen/MeshGeometry/Coons_interface.h"
#include "geometry_data_struct.h"
#include "list_to_matrix.hpp"
#include "matrix_to_list.hpp"

#include "mesh_quality.h"
#include "tetrahedral_mesh.h"
#include "tetrahedral_mesh_parameters.h"
// #include "verdict.h"
#include "discretize_curve.h"
#include "discretize_curve_parameters.h"
#include "discretize_surface.h"
#include "discretize_surface_parameters.h"
#include "geometry_interface_occ.h"
#include "sizing_field_interface.h"

#include "tiger_parameters.h"

void outputSurfaceMesh(const TiGER::SurfaceMesh& surf, const std::string& filename);

int main(int argc, char* argv[])
{
    auto beforeTime = std::chrono::steady_clock::now();
    //                                                                  //
    //                   PART 1: Parse command line arguments          //
    CLI::App app{"AutoGrid Procedure - A program to generate mesh from geometry file using top-down methods."};
    // IO options
    std::string input_geoFile;
    std::string output_vtk;
    bool output_debug = false;
    bool use_stl = false;
    bool no_remesh = false;
    // Geometry options
    bool geo_getname = false;
    bool geo_divideclosed = true;
    bool geo_no_bool = false;
    int geo_heal = 0;
    // stl options
    double stl_deflection = 0.001;
    double stl_angle = 5;
    double stl_min_size = 1e-7;
    bool stl_relative = false;
    bool stl_in_parallel = true;
    bool stl_internal_vertices_mode = true;
    bool stl_control_surface_deflection = true;
    // sizing field options
    double size_scale_factor = 1.0;
    double size_hmin = -1;
    double size_hmax = -1;
    double size_theta = 0.1;
    double size_beta = 1.2;
    double size_proximity = 2;
    // remesh options
    double remesh_epsilon = 1e-3;
    // autogrid options
    int autogrid_max_passes = 5;
    bool autogrid_use_local = false;
    double autogrid_epsilon = 1e-4;
    bool autogrid_use_round = false;

    // Define the command line options
    app.add_option("-i,--input", input_geoFile, "Input geometry file")->required();
    app.add_option("-o,--output", output_vtk, "Output vtk file");
    app.add_flag("-d,--debug", output_debug, "Enable debugging, output intermediate vtk file");
    app.add_flag("--use-stl", use_stl, "Use stl file as input");
    app.add_flag("--no-remesh", no_remesh, "Disable remeshing");

    // Geometry options
    app.add_flag("--geo-getname", geo_getname, "Geometry: get name?");
    app.add_flag("--geo-divideclosed", geo_divideclosed, "Geometry: divide closed surface?");
    app.add_flag("--geo-no-bool", geo_no_bool, "Geometry: boolean operations?");
    app.add_option("--geo-heal", geo_heal, "Geometry: heal options");

    // STL options
    app.add_option("--stl-deflection", stl_deflection, "STL: deflection");
    app.add_option("--stl-angle", stl_angle, "STL: angle");
    app.add_option("--stl-min-size", stl_min_size, "STL: min-size");
    app.add_flag("--stl-relative", stl_relative, "STL: relative");
    app.add_flag("--stl-in-parallel", stl_in_parallel, "STL: in parallel");
    app.add_flag("--stl-internal-vertices-mode", stl_internal_vertices_mode, "STL: internal vertices mode");
    app.add_flag("--stl-control-surface-deflection", stl_control_surface_deflection, "STL: control surface deflection");

    // Sizing field options
    app.add_option("--size-scale-factor", size_scale_factor, "Sizing field: scale factor");
    app.add_option("--size-hmin", size_hmin, "Sizing field: hmin");
    app.add_option("--size-hmax", size_hmax, "Sizing field: hmax");
    app.add_option("--size-theta", size_theta, "Sizing field: theta");
    app.add_option("--size-beta", size_beta, "Sizing field: beta");
    app.add_option("--size-proximity", size_proximity, "Sizing field: proximity");

    // Remesh options
    app.add_option("--remesh-epsilon", remesh_epsilon, "Remesh: tolerance");

    // Autogrid options
    app.add_option("--autogrid-epsilon", autogrid_epsilon, "Autogrid: tolerance");
    app.add_option("--autogrid-max-passes", autogrid_max_passes, "Autogrid: max iteration times");
    app.add_flag("--autogrid-use-local", autogrid_use_local, "Autogrid: use local eps?");
    app.add_flag("--autogrid-use-round", autogrid_use_round, "Autogrid: use round?");

    // Parse command line options
    CLI11_PARSE(app, argc, argv);

    // Handle default output file if not specified
    if( output_vtk.empty() )
    {
        size_t lastDot = input_geoFile.find_last_of('.');
        if( lastDot != std::string::npos )
        {
            output_vtk = input_geoFile.substr(0, lastDot) + "_mesh.vtk";
        }
        else
        {
            output_vtk = input_geoFile + "_mesh.vtk"; // No extension found
        }
    }

    // Debug output file names
    std::string autogrid_output_vtk;
    std::string remesh_output_vtk;
    std::string stl_output_vtk;
    size_t lastDot = input_geoFile.find_last_of('.');
    if( lastDot != std::string::npos )
    {
        autogrid_output_vtk = input_geoFile.substr(0, lastDot) + "_debug_autogrid.vtk";
        remesh_output_vtk = input_geoFile.substr(0, lastDot) + "_debug_remesh.vtk";
        stl_output_vtk = input_geoFile.substr(0, lastDot) + "_debug_stl.vtk";
    }
    else
    {
        autogrid_output_vtk = input_geoFile + "_debug_autogrid.vtk";
        remesh_output_vtk = input_geoFile + "_debug_remesh.vtk";
        stl_output_vtk = input_geoFile + "_debug_stl.vtk";
    }

    // Output parsed values
    spdlog::info("-------- Parse command line arguments--------- ");
    spdlog::info("IO parameters:");
    spdlog::info("\t input_geoFile: {}", input_geoFile);
    spdlog::info("\t output_vtk: {}", output_vtk);
    spdlog::info("\t output_debug: {}", output_debug);
    spdlog::info("\t use_stl: {}", use_stl);
    spdlog::info("\t no_remesh: {}", no_remesh);
    spdlog::info("geometry parameters:");
    spdlog::info("\t geo_getname: {}", geo_getname);
    spdlog::info("\t geo_divideclosed: {}", geo_divideclosed);
    spdlog::info("\t geo_no_bool: {}", geo_no_bool);
    spdlog::info("\t geo_heal: {}", geo_heal);
    spdlog::info("stl parameters:");
    spdlog::info("\t stl_deflection: {}", stl_deflection);
    spdlog::info("\t stl_angle: {}", stl_angle);
    spdlog::info("\t stl_min_size: {}", stl_min_size);
    spdlog::info("\t stl_relative: {}", stl_relative);
    spdlog::info("\t stl_in_parallel: {}", stl_in_parallel);
    spdlog::info("\t stl_internal_vertices_mode: {}", stl_internal_vertices_mode);
    spdlog::info("\t stl_control_surface_deflection: {}", stl_control_surface_deflection);
    spdlog::info("remesh parameters:");
    spdlog::info("\t remesh_epsilon: {}", remesh_epsilon);
    spdlog::info("sizing field parameters:");
    spdlog::info("\t size_scale_factor: {}", size_scale_factor);
    spdlog::info("\t size_hmin: {}", size_hmin);
    spdlog::info("\t size_hmax: {}", size_hmax);
    spdlog::info("\t size_theta: {}", size_theta);
    spdlog::info("\t size_beta: {}", size_beta);
    spdlog::info("\t size_proximity: {}", size_proximity);
    spdlog::info("autogrid parameters:");
    spdlog::info("\t autogrid_epsilon: {}", autogrid_epsilon);
    spdlog::info("\t autogrid_max_passes: {}", autogrid_max_passes);
    spdlog::info("\t autogrid_use_local: {}", autogrid_use_local);
    spdlog::info("\t autogrid_use_round: {}", autogrid_use_round);

    //                                                          //
    //                   PART 2: Read the geometry file         //
    //                       geo->stl/stl-stl                   //
    TiGER::SurfaceMesh stl_surface_out;
    TiGER::Parameters args;
    if( use_stl )
    {
        spdlog::info("---------Read the stl file     stlk->stl--------- \n");
        TiGER::Mesh m;
        TiGER::MESHIO::readVTK(input_geoFile, m, "surface_id");
        TiGER::matrix_to_list<3>(m.Vertex, stl_surface_out.coords);
        TiGER::matrix_to_list<3>(m.Topo, stl_surface_out.tris);
        TiGER::matrix_to_list<1>(m.Masks, stl_surface_out.attribute_int);
        if( output_debug )
        {
            outputSurfaceMesh(stl_surface_out, stl_output_vtk);
        }
        // return 0;
    }
    else
    {
        spdlog::info("---------Read the geometry file    geo->stl--------- \n");

        TiGER::Geometry Geo;

        TiGER::readModelParameters geo_args;
        TiGER::generateStlParameters stl_args;

        geo_args.setdivideclosed_flag(geo_divideclosed);
        geo_args.setgetname_flag(geo_getname);
        geo_args.setbool_flag(!geo_no_bool);
        geo_args.setheal(geo_heal);

        stl_args.setDeflection(stl_deflection);
        stl_args.setAngle(stl_angle);
        stl_args.setRelative(stl_relative);
        stl_args.setInParallel(stl_in_parallel);
        stl_args.setMinSize(stl_min_size);
        stl_args.setInternalVerticesMode(stl_internal_vertices_mode);
        stl_args.setControlSurfaceDeflection(stl_control_surface_deflection);

        TiGER::Read_Model(input_geoFile, Geo, geo_args);
        TiGER::CAD_Tessellation(Geo, stl_surface_out, stl_args);
        // jc 把面id信息写到了regions里面
        stl_surface_out.attribute_int = stl_surface_out.regions;
        //stl_surface_out.regions.clear();
        if( output_debug )
        {
            outputSurfaceMesh(stl_surface_out, stl_output_vtk);
        }
    }

    //                                                             //
    //              PART 3: Use AutoGrid to repair stl            //
    //                 stl->repaired surface mesh                //
    spdlog::info("----Use AutoGrid to repair stl    stl->repaired surface mesh----");
    TiGER::SurfaceMesh autogrid_surf_in, autogrid_surf_out;
    TiGER::AutogridParameters autogrid_args;

    autogrid_args.seteps_rel(autogrid_epsilon);
    autogrid_args.setmax_num_passes(autogrid_max_passes);
    autogrid_args.setuse_local(autogrid_use_local);
    autogrid_args.setuse_round(autogrid_use_round);

    autogrid_surf_in = stl_surface_out;

    spdlog::info("before TMeshSurf_Repair");
    TiGER::discrete_geometry_repair::TMeshSurf_Repair(autogrid_surf_in, autogrid_args, autogrid_surf_out);
    if( output_debug )
    {
        outputSurfaceMesh(autogrid_surf_out, autogrid_output_vtk);
    }
    //                                                              //
    //                  PART 4: Generate size field                //
    //                repaired surface mesh->size field           //

    TiGER::SizingManager sf;
    if( !no_remesh )
    {
        spdlog::info("----Generate size field    repaired surface mesh->size field----");
        spdlog::info("size field...\n");
        TiGER::BackGroundParameters bkg_args;
        bkg_args.sethmin(size_hmin);
        bkg_args.sethmax(size_hmax);
        bkg_args.settheta(size_theta);
        bkg_args.setbeta(size_beta);
        bkg_args.setproximity(size_proximity);

        spdlog::info(sf.sf.size());
        spdlog::info("before setSizeField\n");

        std::cout << autogrid_surf_out.regions.size() << std::endl;
        std::cout << autogrid_surf_out.attribute_int.size() << std::endl;
        autogrid_surf_out.regions = autogrid_surf_out.attribute_int;

        TiGER::size_field::setSizeField(autogrid_surf_out, bkg_args, sf);
        spdlog::info("after setSizeField\n");
        spdlog::info(sf.sf.size());
    }

    //                                                               //
    //                 PART 5:Remesh the surface mesh               //
    //         (surafce mesh+size field)  ->  smoothed tet mesh    //
    TiGER::SurfaceMesh remesh_surf_in = autogrid_surf_out;
    TiGER::SurfaceMesh remesh_surf_out;
    if( !no_remesh )
    {
        spdlog::info("----Remesh the surface mesh (surafce mesh+size field)->smoothed tet ");
        TiGER::RemeshParameters args_remesh;
        args_remesh.setiteration_number(1);
        args_remesh.setb_sizefunction(true);
        args_remesh.seteps_rel(remesh_epsilon);

        std::vector<std::vector<int>> constrained_edge;
        std::vector<std::vector<int>> conformaing_edge;
        std::vector<std::vector<int>> constrained_edge_out;
        std::function<double(double, double, double)> scaled_size_function = [&](double x, double y, double z) {
            return size_scale_factor * TiGER::size_field::SizingFunction_getSizeAtPoint(x, y, z, sf);
        };
        TiGER::discretize_surface_remesh::TMeshSurf_Opt(remesh_surf_in, constrained_edge, conformaing_edge, args_remesh,
                                                        scaled_size_function, remesh_surf_out, constrained_edge_out);
    }
    if( output_debug )
    {
        outputSurfaceMesh(remesh_surf_out, remesh_output_vtk);
    }

    //                                                           //
    //             PART 6:  Generate tetrahedral mesh           //
    //             repaired surface mesh  ->  tet mesh         //
    TiGER::SurfaceMesh dt_surf_in;
    if( no_remesh )
    {
        dt_surf_in = autogrid_surf_out;
    }
    else
    {
        dt_surf_in = remesh_surf_out;
    }
    spdlog::info("----Generate tetrahedral mesh   repaired surface mesh->tet mesh----");
    TiGER::VolumeMesh vol;
    TiGER::tetrahedral_mesh::TetGrid_Constrained(dt_surf_in, args, vol, NULL);

    auto afterTime = std::chrono::steady_clock::now();
    double duration_second = std::chrono::duration<double>(afterTime - beforeTime).count();

    //                                                            //
    //            PART 7: Write the tet mesh to vtk file         //
    //                  tet mesh    ->   vtk                    //
    TiGER::Mesh m;
    TiGER::list_to_matrix<3>(vol.coords, m.Vertex);
    TiGER::list_to_matrix<4>(vol.tetras, m.Topo);
    TiGER::MESHIO::writeVTK(output_vtk, m);

    //          PART 8: Calculate the quality of the tet mesh       //
    //                tet mesh  ->  quality                        //
    spdlog::info("----Calculate the quality of the tet mesh  tet mesh->quality----");
    const TiGER::VolMeshQualityType type = TiGER::VolMeshQualityType::TET_VOLUME;
    std::vector<TiGER::QualityResult> quality;
    TiGER::Tool_VolQuality(vol, type, quality);
    for( auto i : quality )
    {
        spdlog::info(i.ave_value);
        spdlog::info(i.max_value);
        spdlog::info(i.min_value);
    }

    //                                                            //

    spdlog::info("Total time: {}s", duration_second);
    return 0;
}

void outputSurfaceMesh(const TiGER::SurfaceMesh& surf, const std::string& filename)
{
    TiGER::Mesh m;
    TiGER::list_to_matrix<3>(surf.coords, m.Vertex);
    TiGER::list_to_matrix<3>(surf.tris, m.Topo);
    TiGER::list_to_matrix<1>(surf.attribute_int, m.Masks);
    TiGER::MESHIO::writeVTK(filename, m, "surface_id");
}