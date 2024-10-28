message("building test example for tiger 1.5")
project(DT_Tri)

set(EXTERNAL_DIR ${CMAKE_SOURCE_DIR}/extern)
set(DT_SOURCE ${CMAKE_SOURCE_DIR}/src/wyfdt2d)
set(DIR ${CMAKE_SOURCE_DIR}/src/module/DT_Tri)
set(SOURCE ${DIR}/dt_2D.cpp)
include_directories(
  ${CMAKE_SOURCE_DIR}/include
  ${DT_SOURCE}/src
  ${DT_SOURCE}/geom
)

add_library(${PROJECT_NAME} SHARED ${SOURCE})
target_link_libraries(${PROJECT_NAME} PUBLIC dt2D) 

add_executable(${PROJECT_NAME}_test ${CMAKE_SOURCE_DIR}/test/test_DT_Tri.cpp)
target_link_libraries(${PROJECT_NAME}_test PUBLIC ${PROJECT_NAME})
target_include_directories(${PROJECT_NAME}_test PUBLIC ${CMAKE_SOURCE_DIR}/extern/catch2)
target_include_directories(${PROJECT_NAME}_test PUBLIC ${CMAKE_SOURCE_DIR}/extern/eigen)