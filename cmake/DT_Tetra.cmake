message("building test example for tiger 1.5")
project(DT_Tetra)

set(EXTERNAL_DIR "${CMAKE_SOURCE_DIR}/extern")
set(DT_SOURCE ${CMAKE_SOURCE_DIR}/src/wyfdt)
find_package(OpenMP)
set(DIR ${CMAKE_SOURCE_DIR}/src/module/DT_Tetra)
set(SOURCE ${DIR}/dt_mesh.cpp)
include_directories(
  ${CMAKE_SOURCE_DIR}/include
  ${DT_SOURCE}/src
)

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OpenMP found")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_library(${PROJECT_NAME} SHARED ${SOURCE})
target_link_libraries(${PROJECT_NAME} PUBLIC dt) 
target_link_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/lib) 

add_executable(${PROJECT_NAME}_test ${CMAKE_SOURCE_DIR}/test/test_DT_Tetra.cpp)
target_link_libraries(${PROJECT_NAME}_test PUBLIC ${OpenMP_CXX_FLAGS})
target_link_libraries(${PROJECT_NAME}_test PUBLIC TRANS2_Geo)
target_link_libraries(${PROJECT_NAME}_test PUBLIC AFT_Tri)
target_link_libraries(${PROJECT_NAME}_test PUBLIC DT_Remesh)
target_link_libraries(${PROJECT_NAME}_test PUBLIC Quality_Data)
target_link_libraries(${PROJECT_NAME}_test PUBLIC SizingFunc)
target_link_libraries(${PROJECT_NAME}_test PUBLIC ${PROJECT_NAME})
target_include_directories(${PROJECT_NAME}_test PUBLIC ${CMAKE_SOURCE_DIR}/extern/catch2)
target_include_directories(${PROJECT_NAME}_test PUBLIC ${CMAKE_SOURCE_DIR}/test)
target_include_directories(${PROJECT_NAME}_test PUBLIC ${CMAKE_SOURCE_DIR}/extern/eigen)



