message("building test example for tiger 1.5")
#include("${CMAKE_SOURCE_DIR}/cmake/DownloadProject.cmake")
project(DT_Remesh)

set(EXTERNAL_DIR "${CMAKE_SOURCE_DIR}/extern")
include(${CMAKE_SOURCE_DIR}/cmake/OpenCasCadeConfig.cmake)

add_definitions(-DBUILDING_DT_REMESH)

if(OPENMP_FOUND)
    message("OpenMP found")
    
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
find_package(opencascade REQUIRED)

include_directories(${CMAKE_SOURCE_DIR}/include)
set(DIR ${CMAKE_SOURCE_DIR}/src/DT_remesh)
set(SOURCE
    ${DIR}/DT_remesh.cpp
)

add_library(${PROJECT_NAME} SHARED ${SOURCE})
add_executable(${PROJECT_NAME}_test ${CMAKE_SOURCE_DIR}/test/test_dt_remesh.cpp )

if(MSVC)
    target_link_directories(${PROJECT_NAME} PUBLIC ${opencascade_LIBRARIES_DEBUG})
    target_link_libraries(${PROJECT_NAME} PUBLIC ${opencascade_LIBRARIES_DEBUG})
    target_link_directories(${PROJECT_NAME} PUBLIC ${opencascade_INCLUDE_DIRS}/../lib)
        target_link_libraries(${PROJECT_NAME} PUBLIC 
        -Wl,--start-group 
        ${OCC_LIBS}
        -Wl,--end-group)

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
    #target_include_directories(${PROJECT_NAME}_test PUBLIC ${opencascade_INCLUDE_DIRS})
    target_include_directories(${PROJECT_NAME} PUBLIC ${opencascade_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} PUBLIC pthread)
    target_link_libraries(${PROJECT_NAME} PUBLIC dl)
    target_link_directories(${PROJECT_NAME} PUBLIC ${opencascade_INCLUDE_DIRS}/../lib)
    target_link_directories(${PROJECT_NAME} PUBLIC ${opencascade_INCLUDE_DIRS}/../bin)
    target_link_libraries(${PROJECT_NAME} PUBLIC 
        -Wl,--start-group 
        ${OCC_LIBS}
        -Wl,--end-group)
endif()
find_package(libtiger REQUIRED)
target_link_directories(${PROJECT_NAME} PUBLIC ${libtiger_INCLUDE_DIRS}/release/lib)

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/geom)
target_link_libraries(${PROJECT_NAME} PUBLIC tiger_geom) 

find_package(libtiger REQUIRED)
message(${libtiger_INCLUDE_DIRS})
target_link_directories(${PROJECT_NAME} PUBLIC ${libtiger_INCLUDE_DIRS}/release/lib)

target_link_directories(${PROJECT_NAME} PUBLIC ${Smesh_SOURCE}/build/lib)
target_include_directories(${PROJECT_NAME} PUBLIC ${Smesh_SOURCE}/include)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/eigen)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/libigl/include)
target_link_libraries(${PROJECT_NAME} PUBLIC DT_Tri)
target_link_libraries(${PROJECT_NAME} PUBLIC REMesh_Triangle)
target_link_libraries(${PROJECT_NAME} PUBLIC TRANS2_Geo)
target_link_libraries(${PROJECT_NAME} PUBLIC AFT_Tri)


target_link_libraries(${PROJECT_NAME}_test PUBLIC ${PROJECT_NAME})
target_include_directories(${PROJECT_NAME}_test PUBLIC ${Smesh_SOURCE}/include)
target_link_libraries(${PROJECT_NAME}_test PUBLIC ${OpenMP_CXX_FLAGS})
target_link_libraries(${PROJECT_NAME}_test PUBLIC Quality_Data)

target_include_directories(${PROJECT_NAME}_test PUBLIC ${CMAKE_SOURCE_DIR}/extern/catch2)
target_include_directories(${PROJECT_NAME}_test PUBLIC ${CMAKE_SOURCE_DIR}/extern/eigen)
target_include_directories(${PROJECT_NAME}_test PUBLIC ${CMAKE_SOURCE_DIR}/extern/libigl/include)

message(STATUS "subProject Name: ${PROJECT_NAME}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG "${PROJECT_SOURCE_DIR}/build/Debug")
