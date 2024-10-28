message("building test example for tiger 1.5")
project(REMesh_Triangle)

set(EXTERNAL_DIR "${CMAKE_SOURCE_DIR}/extern")

set(Remesh_SOURCE ${CMAKE_SOURCE_DIR}/src/remesh)
include_directories(${CMAKE_SOURCE_DIR}/include)
set(DIR ${CMAKE_SOURCE_DIR}/src/module/REMesh_Triangle)
set(SOURCE
    ${DIR}/discretize_surface_remesh.cpp
)
add_library(${PROJECT_NAME} SHARED ${SOURCE})
target_compile_definitions(${PROJECT_NAME} PRIVATE Tiger_EXPORTS)
target_include_directories(${PROJECT_NAME} PUBLIC ${Remesh_SOURCE}/src)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/cli11/include)
target_link_libraries(${PROJECT_NAME} PUBLIC LTRemesh) 
find_package(OpenMP)
target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)

add_executable(${PROJECT_NAME}_test ${CMAKE_SOURCE_DIR}/test/test_${PROJECT_NAME}.cpp)
target_link_libraries(${PROJECT_NAME}_test PUBLIC ${PROJECT_NAME})
target_include_directories(${PROJECT_NAME}_test PUBLIC ${CMAKE_SOURCE_DIR}/extern/catch2)
target_include_directories(${PROJECT_NAME}_test PUBLIC ${CMAKE_SOURCE_DIR}/extern/eigen)
