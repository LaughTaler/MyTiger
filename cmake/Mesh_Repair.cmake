cmake_minimum_required(VERSION 3.0)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OpenMP found")

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

project(Mesh_Repair)

if(UNIX)
add_compile_options(-fPIC)
endif()

find_package(libtiger REQUIRED)
set(REPAIR ${CMAKE_SOURCE_DIR}/src/${PROJECT_NAME})
set(SOURCE
${REPAIR}/src/discrete_geometry_repair.cpp
)
add_library(${PROJECT_NAME} SHARED ${SOURCE})
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/eigen)
target_link_directories(${PROJECT_NAME} PUBLIC ${libtiger_INCLUDE_DIR}/release/lib)


target_compile_definitions(${PROJECT_NAME} PRIVATE Tiger_EXPORTS)
# target_link_libraries(${PROJECT_NAME} PUBLIC Quality_Data)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/catch2)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/eigen)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/geom)

target_link_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/tiger)

# add_subdirectory(${CMAKE_SOURCE_DIR}/extern/geom) 

target_link_libraries(${PROJECT_NAME} PUBLIC tiger_geom)  # 需要链接geom库
target_link_libraries(${PROJECT_NAME} PUBLIC Autogrid_Triangle)
# target_link_libraries(${PROJECT_NAME} PUBLIC Quality_Data)


# add_executable(${PROJECT_NAME}_test ${SOURCE} ${CMAKE_SOURCE_DIR}/test/test_mesh_repair.cpp)
add_executable(${PROJECT_NAME}_test ${CMAKE_SOURCE_DIR}/test/test_Mesh_Repair.cpp)
target_link_libraries(${PROJECT_NAME}_test PUBLIC ${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME}_test PUBLIC ${OpenMP_CXX_FLAGS})
# # target_link_libraries(${PROJECT_NAME}_test PUBLIC Quality_Data)
target_include_directories(${PROJECT_NAME}_test PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_include_directories(${PROJECT_NAME}_test PUBLIC ${CMAKE_SOURCE_DIR}/extern/catch2)
target_include_directories(${PROJECT_NAME}_test PUBLIC ${CMAKE_SOURCE_DIR}/extern/eigen)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/geom)