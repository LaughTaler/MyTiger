message("building turbomachinery example for tiger 1.5")
project(turbomachinery_test)

set(EXTERN_DIR ${CMAKE_SOURCE_DIR}/extern)
set(SOURCE ${CMAKE_SOURCE_DIR}/test/test_turbomachinery.cpp)

find_package(opencascade REQUIRED)
find_package(libtiger REQUIRED)
add_executable(${PROJECT_NAME} ${SOURCE})



if(WIN32)

# 获取所有 DLL 文件
file(GLOB DLL_FILES "${libtiger_INCLUDE_DIR}/release/bin/*.dll" "${opencascade_INCLUDE_DIRS}/../bin/*.dll")

# 复制所有 DLL 文件到目标目录
foreach(DLL_FILE ${DLL_FILES})
    file(COPY "${DLL_FILE}" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
endforeach()
    
endif()


if(WIN32)
    find_package(cgns REQUIRED)
    target_include_directories(${PROJECT_NAME} PUBLIC ${cgns_INCLUDE_DIRS})
    target_link_directories(${PROJECT_NAME} PUBLIC ${cgns_INCLUDE_DIRS}/../lib)
    # 获取Windows下的库文件列表
    file(GLOB LIBRARY_FILES "${cgns_INCLUDE_DIRS}/../lib/*.lib")
    foreach(lib_file ${LIBRARY_FILES})
        get_filename_component(lib_name ${lib_file} NAME_WE)
        target_link_libraries(${PROJECT_NAME} PUBLIC ${lib_name})
    endforeach()
elseif(UNIX)
    find_package(cgns REQUIRED)
    target_include_directories(${PROJECT_NAME} PUBLIC ${cgns_INCLUDE_DIRS})
	target_link_directories(${PROJECT_NAME} PUBLIC ${cgns_INCLUDE_DIRS}/../lib)
	target_link_libraries(${PROJECT_NAME} PUBLIC libcgns.a)
    target_link_libraries(${PROJECT_NAME} PUBLIC libhdf5.a)
    install(TARGETS test_CGNS_Write DESTINATION bin)
endif()



target_include_directories(${PROJECT_NAME} PUBLIC 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/tools
    ${EXTERN_DIR}/eigen
    ${EXTERN_DIR}/cli11
    ${EXTERN_DIR}/libigl/include
    ${EXTERN_DIR}/spdlog/include
    ${EXTERN_DIR}/cgns/include
    ${EXTERN_DIR}/geom
    ${EXTERN_DIR}/MeshOrient
) 





#TODO 标准化extern
add_subdirectory(extern/MeshOrient)

target_link_libraries(${PROJECT_NAME}
PUBLIC
        $<$<CONFIG:Debug>:${opencascade_LIBRARIES_DEBUG}>     
        $<$<CONFIG:Release>:${opencascade_LIBRARIES_RELEASE}> 
)



target_link_directories(${PROJECT_NAME} PUBLIC ${OCC_DIR}/lib )
target_link_directories(${PROJECT_NAME} PUBLIC ${opencascade_INCLUDE_DIRS}/../lib)
target_link_directories(${PROJECT_NAME} PUBLIC ${libtiger_INCLUDE_DIR}/release/lib)

target_link_libraries(${PROJECT_NAME} PUBLIC ${opencascade_LIBRARIES_DEBUG})
target_link_libraries(${PROJECT_NAME} PUBLIC ${OCC_LIBS})
target_link_libraries(${PROJECT_NAME} PUBLIC DT_Tetra)
target_link_libraries(${PROJECT_NAME} PUBLIC DT_Tri)
target_link_libraries(${PROJECT_NAME} PUBLIC TRANS2_Geo)
target_link_libraries(${PROJECT_NAME} PUBLIC AFT_Tri)
target_link_libraries(${PROJECT_NAME} PUBLIC Quality_Data)
target_link_libraries(${PROJECT_NAME} PUBLIC SizingFunc)
target_link_libraries(${PROJECT_NAME} PUBLIC REMesh_Triangle)
target_link_libraries(${PROJECT_NAME} PUBLIC Mesh_Repair)
target_link_libraries(${PROJECT_NAME} PUBLIC ALM_Hybrid)
target_link_libraries(${PROJECT_NAME} PUBLIC tiger_geom)
target_link_libraries(${PROJECT_NAME} PUBLIC MeshOrient)



