message("building test example for tiger 1.5")
project(ALM_Hybrid)
FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OpenMP found")

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
add_executable(${PROJECT_NAME}_test ${CMAKE_SOURCE_DIR}/test/test_ALM_Hybrid.cpp)
target_link_libraries(${PROJECT_NAME}_test PUBLIC ${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME}_test PUBLIC DT_Tetra)
target_include_directories(${PROJECT_NAME}_test PUBLIC ${CMAKE_SOURCE_DIR}/extern/catch2)
target_include_directories(${PROJECT_NAME}_test PUBLIC ${CMAKE_SOURCE_DIR}/extern/eigen)
target_link_libraries(${PROJECT_NAME}_test PUBLIC ${OpenMP_CXX_FLAGS})


set(BLMESH_SOURCE ${CMAKE_SOURCE_DIR}/src/blmesh)
include_directories(${BLMESH_SOURCE}/vmesh/include)
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/extern/spdlog/include)
set(DIR ${CMAKE_SOURCE_DIR}/src/module/ALM_Hybrid)
set(SOURCE
    ${DIR}/boundary_layer_mesh.cpp
)
add_library(${PROJECT_NAME} SHARED ${SOURCE})
target_link_libraries(${PROJECT_NAME} PUBLIC vmesh) 
#target_link_libraries(${PROJECT_NAME} PUBLIC TigerAdasize)
target_link_libraries(${PROJECT_NAME} PUBLIC libblmesh)
#target_link_libraries(${PROJECT_NAME} PUBLIC create_sol)
target_link_libraries(${PROJECT_NAME} PUBLIC tetgen)
target_link_libraries(${PROJECT_NAME} PUBLIC DT_Tetra)
target_link_libraries(${PROJECT_NAME} PUBLIC triangle)
target_link_libraries(${PROJECT_NAME} PUBLIC 2dremesh)
target_link_libraries(${PROJECT_NAME} PUBLIC blpre libdtiso3d
 libdtiso2d libgeom meshOrientLib libspr libMeshChamfer
  DistanceCaculator libpostprocess)

target_link_directories(${PROJECT_NAME} PUBLIC ${BLMESH_SOURCE}/tiger/lib)
target_link_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/lib)



