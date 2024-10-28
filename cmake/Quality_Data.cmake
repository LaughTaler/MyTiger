cmake_minimum_required(VERSION 3.0)

project(Quality_Data)


set(MM_SOURCE ${CMAKE_SOURCE_DIR}/src/meshmetric)

include_directories(
${MM_SOURCE}/MeshMetric
${MM_SOURCE}/verdict-1.0.1
${CMAKE_SOURCE_DIR}/include
)

set(DIR ${CMAKE_SOURCE_DIR}/src/module/Quality_Data)
set(SOURCE
    ${DIR}/mesh_quality.cpp
)
add_library(${PROJECT_NAME} SHARED ${SOURCE})
target_link_libraries(${PROJECT_NAME} PUBLIC verdict101)

#test
add_executable(${PROJECT_NAME}_test ${CMAKE_SOURCE_DIR}/test/test_${PROJECT_NAME}.cpp)
target_include_directories(${PROJECT_NAME}_test PUBLIC ${CMAKE_SOURCE_DIR}/extern/catch2)
target_include_directories(${PROJECT_NAME}_test PUBLIC ${CMAKE_SOURCE_DIR}/extern/libigl/include)
target_link_libraries(${PROJECT_NAME}_test PUBLIC ${PROJECT_NAME})

