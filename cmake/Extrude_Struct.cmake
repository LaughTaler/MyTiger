message("building test example for tiger 1.5")
# 设置项目名称
project(Extrude_Struct)

include_directories(${CMAKE_SOURCE_DIR}/include)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# 添加可执行文件
#aux_source_directory(${PROJECT_SOURCE_DIR}/src SRC_LIST)
set(FILE_PATH ${PROJECT_SOURCE_DIR}/src/extrude3d/src)
set(SRC_LIST 
${FILE_PATH}/extrude.cpp)
#file(GLOB_RECURSE extrude3d_lib "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

add_library(${PROJECT_NAME} SHARED  ${SRC_LIST})

# 如果有其他头文件目录，可以使用include_directories命令添加
include_directories(${CMAKE_SOURCE_DIR}/extern/eigen)
include_directories(${CMAKE_SOURCE_DIR}/extern/cli11/include)
include_directories(${CMAKE_SOURCE_DIR}/extern/libigl/include)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/geom)

# 如果有其他依赖库，可以使用target_link_libraries命令添加
target_link_libraries(${PROJECT_NAME} PUBLIC tiger_geom)
target_link_libraries(${PROJECT_NAME} PUBLIC Quality_Data)

add_executable(${PROJECT_NAME}_test ${CMAKE_SOURCE_DIR}/test/test_${PROJECT_NAME}.cpp)
target_include_directories(${PROJECT_NAME}_test PUBLIC ${CMAKE_SOURCE_DIR}/extern/catch2)
target_link_libraries(${PROJECT_NAME}_test PUBLIC ${PROJECT_NAME})
target_include_directories(${PROJECT_NAME}_test PUBLIC ${CMAKE_SOURCE_DIR}/extern/geom)

if(UNIX)
    target_link_libraries(${PROJECT_NAME}_test PUBLIC -lpthread)
endif()


