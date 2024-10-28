message("building test example for tiger 1.5")
project(AFT_Tri)
set(CMAKE_DEBUG_POSTFIX,"_d")
set(CMAKE_RELEASE_POSTFIX,"_r")

include(${CMAKE_SOURCE_DIR}/cmake/OpenCasCadeConfig.cmake)

set(EXTERNAL_DIR "${CMAKE_SOURCE_DIR}/extern")

set(Smesh_SOURCE ${CMAKE_SOURCE_DIR}/src/meshgen)

if(OPENMP_FOUND)
    message("OpenMP found")

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
include_directories(${CMAKE_SOURCE_DIR}/include)
set(DIR ${CMAKE_SOURCE_DIR}/src/module/AFT_Tri)
set(SOURCE
    ${DIR}/AFT_triangle_mesh.cpp
)

add_executable(${PROJECT_NAME}_test ${CMAKE_SOURCE_DIR}/test/test_AFT_Tri.cpp)

if(MSVC)
    target_link_directories(${PROJECT_NAME}_test PUBLIC ${OCC_DIR}/lib)
    target_link_libraries(${PROJECT_NAME}_test PUBLIC ${OCC_LIBS})
else()
    #if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm|ARM")
    #target_link_directories(${PROJECT_NAME} PUBLIC ${OCC_DIR}/lib/ARM)
    #else()
    #target_link_directories(${PROJECT_NAME} PUBLIC ${OCC_DIR}/lib/linux)
    #message("OCC LIB PATH=${OCC_DIR}/lib/linux")
    #endif()
    #target_link_libraries(
    ##    ${PROJECT_NAME}
    #    -Wl,--start-group
    #    ${OCC_LIBS}
    #    -Wl,--end-group)
    #find_package(OpenCasCade REQUIRED)
    #target_include_directories(${PROJECT_NAME}_test PUBLIC ${opencascade_INCLUDE_DIRS})
    target_include_directories(${PROJECT_NAME}_test PUBLIC ${opencascade_INCLUDE_DIRS})

    target_link_libraries(${PROJECT_NAME}_test PUBLIC 
        -Wl,--start-group 
        ${OCC_LIBS}
        -Wl,--end-group)
endif()


message("source dir=${DIR}/AFT_triangle_mesh.cpp")

add_library(${PROJECT_NAME} SHARED ${SOURCE})
target_link_directories(${PROJECT_NAME} PUBLIC ${Smesh_SOURCE}/build/lib)
target_include_directories(${PROJECT_NAME} PUBLIC ${Smesh_SOURCE}/include)

target_include_directories(${PROJECT_NAME} PUBLIC ${Smesh_SOURCE}/include/SmeshGen/MeshData)
target_include_directories(${PROJECT_NAME} PUBLIC ${Smesh_SOURCE}/include/SmeshGen/MeshGen) 

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/eigen)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/geom)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/libigl/include)

target_link_libraries(${PROJECT_NAME} PUBLIC REMesh_Triangle)
target_link_libraries(${PROJECT_NAME} PUBLIC SmeshGen_lib)
target_link_libraries(${PROJECT_NAME} PUBLIC tiger_geom)


target_include_directories(${PROJECT_NAME}_test PUBLIC ${Smesh_SOURCE}/include)
target_link_libraries(${PROJECT_NAME}_test PUBLIC ${PROJECT_NAME}) 
target_link_libraries(${PROJECT_NAME}_test PUBLIC ${OpenMP_CXX_FLAGS})
target_link_libraries(${PROJECT_NAME}_test PUBLIC TRANS2_Geo)
target_link_libraries(${PROJECT_NAME}_test PUBLIC Quality_Data)
target_link_libraries(${PROJECT_NAME}_test PUBLIC DT_Remesh)
target_link_libraries(${PROJECT_NAME}_test PUBLIC SizingFunc)
target_include_directories(${PROJECT_NAME}_test PUBLIC ${CMAKE_SOURCE_DIR}/extern/catch2)
target_include_directories(${PROJECT_NAME}_test PUBLIC ${CMAKE_SOURCE_DIR}/extern/eigen)
target_include_directories(${PROJECT_NAME}_test PUBLIC ${CMAKE_SOURCE_DIR}/extern/libigl/include)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG "${PROJECT_SOURCE_DIR}/build/Debug")
