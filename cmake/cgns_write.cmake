message("building cgns test example for tiger 1.5")
project(test_CGNS_Write)

include_directories(${CMAKE_SOURCE_DIR}/tools)
include_directories(${CMAKE_SOURCE_DIR}/include)

set(EXTERNAL_DIR "${CMAKE_SOURCE_DIR}/extern")
include_directories(${EXTERN_DIR}/eigen)

set(SOURCE ${CMAKE_SOURCE_DIR}/test/test_CGNS_Write.cpp)
add_executable(${PROJECT_NAME} ${SOURCE})

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

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/eigen)

install(TARGETS ${PROJECT_NAME}
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib
RUNTIME DESTINATION bin
)


