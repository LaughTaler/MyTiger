message("building test example for tiger 1.5")
project(Autogrid_Triangle)

set(EXTERNAL_DIR "${CMAKE_SOURCE_DIR}/extern")

set(Autogrid_SOURCE ${CMAKE_SOURCE_DIR}/src/autogrid)
include_directories(${CMAKE_SOURCE_DIR}/include)
set(DIR ${CMAKE_SOURCE_DIR}/src/module/Autogrid_Triangle)
set(SOURCE
    ${DIR}/autogrid.cpp
)

find_package(opencascade REQUIRED)
find_package(libtiger REQUIRED)
message(${opencascade_INCLUDE_DIRS})
message(${libtiger_INCLUDE_DIRS})
add_library(${PROJECT_NAME} SHARED ${SOURCE})
target_compile_definitions(${PROJECT_NAME} PRIVATE Tiger_EXPORTS)
target_include_directories(${PROJECT_NAME} PUBLIC ${Autogrid_SOURCE}/src)
target_include_directories(${PROJECT_NAME} PUBLIC ${Autogrid_SOURCE}/src/tetwild)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/tools)
# target_link_directories(${PROJECT_NAME} PUBLIC ${Autogrid_SOURCE}/lib)
target_link_libraries(${PROJECT_NAME} PUBLIC libAutoGrid) 
# find_package(OpenMP)
# target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)

add_executable(${PROJECT_NAME}_test ${CMAKE_SOURCE_DIR}/test/test_${PROJECT_NAME}.cpp)
target_link_directories(${PROJECT_NAME}_test PUBLIC ${Autogrid_SOURCE}/lib)
target_link_libraries(${PROJECT_NAME}_test PUBLIC ${PROJECT_NAME} libAutoGrid)
target_include_directories(${PROJECT_NAME}_test PUBLIC ${CMAKE_SOURCE_DIR}/extern/catch2)
target_include_directories(${PROJECT_NAME}_test PUBLIC ${CMAKE_SOURCE_DIR}/extern/eigen)
