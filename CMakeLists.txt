cmake_minimum_required(VERSION 3.6)
if(MSVC)
    add_compile_options(/W0)
    add_compile_options(/Z7)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Z7")
endif()
if (MSVC)
    # 对于MSVC编译器，可以禁用所有警告
    add_compile_options(/W0)
elseif (CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    # 对于GNU编译器，可以使用 -w 参数来关闭所有警告
    add_compile_options(-w)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # 对于Clang编译器，可以使用 -Wno-all 来关闭所有警告
    add_compile_options(-Wno-all)
endif()
project(TiGER)
set(CMAKE_DEBUG_POSTFIX,"_d")
set(CMAKE_RELEASE_POSTFIX,"_r")
include_directories(${CMAKE_SOURCE_DIR}/tools)
set(EXTERN_DIR ${CMAKE_SOURCE_DIR}/extern)
include_directories(${EXTERN_DIR}/eigen)
include_directories(${EXTERN_DIR}/cli11)
include_directories(${EXTERN_DIR}/libigl/include)
include_directories(${EXTERN_DIR}/catch2)
# include_directories(${EXTERN_DIR}/geom)
set(server_IP "10.12.220.129")
include(${CMAKE_SOURCE_DIR}/cmake/add_module.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/add_internal_module.cmake)

#include(${CMAKE_SOURCE_DIR}/cmake/update_submodule.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/conan_config.cmake)


if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR})
include(${CMAKE_SOURCE_DIR}/cmake/add_extern.cmake)
add_extern_git(catch2 http://${server_IP}/opensource/catch2.git 78e71924)
add_extern_git(cli11 http://${server_IP}/opensource/cli11.git e608a13a)
add_extern_git(eigen http://${server_IP}/opensource/eigen.git 894d9c0b)
add_extern_git(libigl http://${server_IP}/opensource/libigl.git 79dc4f68)
add_extern_git(spdlog http://${server_IP}/opensource/spdlog.git 729d7f6d)
add_extern_git(simpleini http://${server_IP}/opensource/simpleini.git f7862c3dd7)
add_extern_git(googletest http://${server_IP}/opensource/googletest.git c4b7c6c3)
add_subdirectory(${CMAKE_SOURCE_DIR}/extern/googletest)
add_extern_git(geom http://${server_IP}/opensource/geom.git 9e690b6e89d6090f371)
add_subdirectory(${CMAKE_SOURCE_DIR}/extern/geom)

add_module(ALM_Hybrid blmesh http://${server_IP}/alm_hybrid/light_blmesh.git   f728cb1242217c246112488429464fa22f64d2d1)
add_module(DT_Tetra wyfdt http://${server_IP}/dt_tetra/wyfdt.git  bdb69e88287f57eacc2cab43cc0c0fe8dfab8cf4)
add_module(DT_Tri wyfdt2d http://${server_IP}/opensource/wyfdt2d.git db8beee2a2)

add_module(TRANS2_Geo trans2_geo http://${server_IP}/LiwGeom/trans2_geo.git  ffe5eba0e08e20acc9b8693e0bf3e75338c71430)

add_module(REMesh_Triangle remesh http://${server_IP}/remesh_triangle/ltremesh.git  1163f6e9ada8a3f7ffff1441595)

add_module(Extrude_Struct extrude3d http://${server_IP}/tiger/extrude3d.git 2545f4b)

add_module(Quality_Data meshmetric http://${server_IP}/tiger/meshmetric.git 1914c360)

add_module(AFT_Tri meshgen http://${server_IP}/aft_tri/meshgen.git 6bd0c19f)

add_module(SizingFunc sizingfunction http://${server_IP}/sizefunction/emeshgensizingfunction.git 68363086)
add_internal_module(Mesh_Repair http://${server_IP}/tiger/mesh_repair.git  1679152f)
add_module(DT_Remesh DT_remesh http://${server_IP}/tiger/dt_remesh.git da348b993a7805351ef537e6cc67faccd7239f20)

add_module(Autogrid_Triangle autogrid http://${server_IP}/tmeshsurf/autogrid.git 876601ff97f)


# extern
#include(${CMAKE_SOURCE_DIR}/cmake/TRANS2_Geo.cmake)
# include(${CMAKE_SOURCE_DIR}/cmake/dt_remesh.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/Generate_Line.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/cgns_write.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/tigerdemo_legacy2.cmake)
# include(${CMAKE_SOURCE_DIR}/cmake/test_api.cmake)
# include(${CMAKE_SOURCE_DIR}/cmake/test_aft_tri.cmake)
# include(${CMAKE_SOURCE_DIR}/cmake/test_dt_remesh.cmake)
# include(${CMAKE_SOURCE_DIR}/cmake/REMesh_Triangle.cmake)
# include(${CMAKE_SOURCE_DIR}/cmake/SizingFunc.cmake)
#include(${CMAKE_SOURCE_DIR}/cmake/turbomachinery.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/tigerdemo_AutoGrid.cmake)



install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/
        DESTINATION include
        FILES_MATCHING PATTERN "*.h")
