message("building tigerdemo for tiger 1.5")
project(tigerdemo_legacy)

set(EXTERNAL_DIR "${CMAKE_SOURCE_DIR}/extern")
include_directories(${EXTERN_DIR}/geom)

set(OCC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/extern/OpenCasCade)
# set(Tran2fli_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/trans2_geo)
set(SOURCE 
${CMAKE_SOURCE_DIR}/test/test_legacy.cpp
)
# set(Tran2fli_SOURCE ${CMAKE_SOURCE_DIR}/src/trans2_geo)

if(OPENMP_FOUND)
    message("OpenMP found")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(opencascade REQUIRED)
find_package(libtiger REQUIRED)
add_executable(${PROJECT_NAME} ${SOURCE})
include(${CMAKE_SOURCE_DIR}/cmake/OpenCasCadeConfig.cmake)

target_include_directories(
    ${PROJECT_NAME} PUBLIC 
    ${CMAKE_SOURCE_DIR}/include
    ${OCC_DIR}/include
    ${Tran2fli_SOURCE}
)

# link_directories(${CMAKE_SOURCE_DIR}/extern/tiger)  # 假设静态库位于构建目录下


if(MSVC)
    target_link_directories(${PROJECT_NAME} PUBLIC ${opencascade_LIBRARIES_DEBUG})
    target_link_libraries(${PROJECT_NAME} PUBLIC ${opencascade_LIBRARIES_DEBUG})
    target_link_directories(${PROJECT_NAME} PUBLIC ${opencascade_INCLUDE_DIRS}/../lib)
    target_link_directories(${PROJECT_NAME} PUBLIC ${libtiger_INCLUDE_DIR}/release/lib)
    message("${PROJECT_NAME} PUBLIC ${libtiger_INCLUDE_DIR}/release/lib")
    target_link_directories(${PROJECT_NAME} PUBLIC ${OCC_DIR}/lib)
    message("${PROJECT_NAME} PUBLIC ${OCC_DIR}/lib")
    target_link_libraries(${PROJECT_NAME} PUBLIC ${OCC_LIBS})

    target_link_libraries(${PROJECT_NAME} PUBLIC TRANS2_Geo)
    target_link_libraries(${PROJECT_NAME} PUBLIC AFT_Tri)
    # target_link_libraries(${PROJECT_NAME} PUBLIC ALM_Hybrid)
    target_link_libraries(${PROJECT_NAME} PUBLIC Quality_Data)
    target_link_libraries(${PROJECT_NAME} PUBLIC SizingFunc)
    target_link_libraries(${PROJECT_NAME} PUBLIC REMesh_Triangle)
    target_link_libraries(${PROJECT_NAME} PUBLIC DT_Tetra)
    target_link_libraries(${PROJECT_NAME} PUBLIC DT_Tri)
    target_link_libraries(${PROJECT_NAME} PUBLIC Mesh_Repair)
    target_link_libraries(${PROJECT_NAME} PUBLIC tiger_geom)  # 需要链接geom库

else()
    find_package(OpenCasCade REQUIRED)
    # #target_include_directories(${PROJECT_NAME}_test PUBLIC ${opencascade_INCLUDE_DIRS})
    target_include_directories(${PROJECT_NAME} PUBLIC ${opencascade_INCLUDE_DIRS})

    target_link_libraries(${PROJECT_NAME} PUBLIC pthread)
    target_link_libraries(${PROJECT_NAME} PUBLIC dl)
    target_link_directories(${PROJECT_NAME} PUBLIC ${opencascade_INCLUDE_DIRS}/../lib)
    message("${PROJECT_NAME} PUBLIC ${opencascade_INCLUDE_DIRS}/../lib")
    target_link_directories(${PROJECT_NAME} PUBLIC ${opencascade_INCLUDE_DIRS}/../bin)
    target_link_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/extern/OpenCasCade/lib)
    
    target_link_libraries(${PROJECT_NAME} PUBLIC TRANS2_Geo)

    target_link_libraries(${PROJECT_NAME} PUBLIC 
        -Wl,--start-group 
        ${OCC_LIBS}
        -Wl,--end-group
        )
   
    
    # target_link_libraries(${PROJECT_NAME} PUBLIC zgrid_geometry0821)
    target_link_directories(${PROJECT_NAME} PUBLIC ${libtiger_INCLUDE_DIR}/release/lib)
    message("${PROJECT_NAME} PUBLIC ${libtiger_INCLUDE_DIR}/release/lib")
    target_link_directories(${PROJECT_NAME} PUBLIC ${OCC_DIR}/lib)
    message("${PROJECT_NAME} PUBLIC ${OCC_DIR}/lib")
    target_link_libraries(${PROJECT_NAME} PUBLIC ${OCC_LIBS})



    target_link_libraries(${PROJECT_NAME} PUBLIC TRANS2_Geo)
    target_link_libraries(${PROJECT_NAME} PUBLIC AFT_Tri)
    # target_link_libraries(${PROJECT_NAME} PUBLIC ALM_Hybrid)
    target_link_libraries(${PROJECT_NAME} PUBLIC Quality_Data)
    target_link_libraries(${PROJECT_NAME} PUBLIC SizingFunc)
    target_link_libraries(${PROJECT_NAME} PUBLIC REMesh_Triangle)
    target_link_libraries(${PROJECT_NAME} PUBLIC DT_Tetra)
    target_link_libraries(${PROJECT_NAME} PUBLIC DT_Tri)
    target_link_libraries(${PROJECT_NAME} PUBLIC Mesh_Repair)
    target_link_libraries(${PROJECT_NAME} PUBLIC tiger_geom)  # 需要链接geom库

    
    
    
    
    # target_link_libraries(${PROJECT_NAME} PUBLIC opencascade::opencascade)

    target_include_directories(${PROJECT_NAME} PUBLIC ${opencascade_INCLUDE_DIRS})
    target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern)
    target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/simpleini)
    target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/spdlog/include)
    # target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/cli11/include)
    target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/libigl/include)


    # link_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/tiger)
    # target_link_libraries(${PROJECT_NAME} PUBLIC adaSurfSizing_new)
    # target_link_libraries(${PROJECT_NAME} PUBLIC verdict101)
    # target_link_libraries(${PROJECT_NAME} PUBLIC Autogrid_Triangle)
    # target_link_libraries(${PROJECT_NAME} PUBLIC Mesh_Repair_r)
    target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/OpenCasCade/include) 
    


    target_link_libraries(${PROJECT_NAME} PUBLIC ${OpenMP_CXX_FLAGS})
    target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/simpleini)
    target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/spdlog/include)

    # target_link_libraries(${PROJECT_NAME} PUBLIC /home/gridteam/.conan2/p/opence016dd5ecaa39/p/include/../lib/libTKShHealing.so)
endif()

install(TARGETS ${PROJECT_NAME}
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib
RUNTIME DESTINATION bin
)

