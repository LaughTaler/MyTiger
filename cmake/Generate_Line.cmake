message("building line mesh test example for tiger 1.5")
project(test_Generate_Line)

include_directories(${CMAKE_SOURCE_DIR}/tools)
include_directories(${CMAKE_SOURCE_DIR}/include)

set(EXTERNAL_DIR "${CMAKE_SOURCE_DIR}/extern")
include_directories(${EXTERN_DIR}/eigen)

set(SOURCE ${CMAKE_SOURCE_DIR}/test/test_Generate_Line.cpp)
add_executable(${PROJECT_NAME} ${SOURCE})

find_package(libtiger REQUIRED)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_link_directories(${PROJECT_NAME} PUBLIC ${libtiger_INCLUDE_DIRS}/debug/lib)
    target_link_libraries(${PROJECT_NAME} PUBLIC AFT_Tri)
    set_target_properties(${PROJECT_NAME} PROPERTIES 
        VS_DEBUGGER_ENVIRONMENT "PATH=${libtiger_INCLUDE_DIRS}/debug/bin;$ENV{PATH}")
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_link_directories(${PROJECT_NAME} PUBLIC ${libtiger_INCLUDE_DIRS}/release/lib)
    target_link_libraries(${PROJECT_NAME} PUBLIC AFT_Tri)
    set_target_properties(${PROJECT_NAME} PROPERTIES 
        VS_DEBUGGER_ENVIRONMENT "PATH=${libtiger_INCLUDE_DIRS}/release/bin;$ENV{PATH}")
endif()

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/eigen)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/catch2)