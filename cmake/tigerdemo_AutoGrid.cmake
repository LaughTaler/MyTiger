message("building tigerdemo for tiger 1.5")
project(tigerdemo_AutoGrid)

set(EXTERNAL_DIR "${CMAKE_SOURCE_DIR}/extern")

set(OCC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/extern/OpenCasCade)
set(Tran2fli_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/trans2_geo)

set(SOURCE 
${CMAKE_SOURCE_DIR}/test/test_AutoGrid.cpp
)

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OpenMP found")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include(${CMAKE_SOURCE_DIR}/cmake/OpenCasCadeConfig.cmake)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(${PROJECT_NAME} ${SOURCE})

find_package(OpenCasCade REQUIRED)
find_package(libtiger REQUIRED)

if(NOT libtiger_FOUND)
    message(FATAL_ERROR "libtiger_FOUND not found. Please install libtiger_FOUND and set libtiger_FOUND_DIR.")
endif()
message("libtiger_INCLUDE_DIRS:=${libtiger_INCLUDE_DIRS}") 
# message("libtiger_FOUND_LIB_DIRS_DEBUG:=${libtiger_FOUND_LIB_DIRS_DEBUG}") 

message(DIRS:${libtiger_INCLUDE_DIRS})

target_include_directories(${PROJECT_NAME} PUBLIC ${opencascade_INCLUDE_DIRS})
target_link_directories(${PROJECT_NAME} PUBLIC ${libtiger_INCLUDE_DIRS}/release/lib)
if(UNIX)
target_link_libraries(${PROJECT_NAME} PUBLIC pthread)
target_link_libraries(${PROJECT_NAME} PUBLIC dl)
endif()

target_link_directories(${PROJECT_NAME} PUBLIC ${opencascade_INCLUDE_DIRS}/../lib)

message("opencascade_INCLUDE_DIRS:=${opencascade_INCLUDE_DIRS}")
target_link_libraries(${PROJECT_NAME} PUBLIC 
    -Wl,--start-group 
    ${OCC_LIBS}
    -Wl,--end-group
)
target_compile_definitions(${PROJECT_NAME} PRIVATE Tiger_EXPORTS)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/simpleini)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/spdlog/include)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/cli11)
# target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/cxxopts/include)
# target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/cli11)

target_link_libraries(${PROJECT_NAME} PUBLIC DT_Tetra)
target_link_libraries(${PROJECT_NAME} PUBLIC DT_Tri)
target_link_libraries(${PROJECT_NAME} PUBLIC TRANS2_Geo)
target_link_libraries(${PROJECT_NAME} PUBLIC Quality_Data)
target_link_libraries(${PROJECT_NAME} PUBLIC SizingFunc)
target_link_libraries(${PROJECT_NAME} PUBLIC ${OpenMP_CXX_FLAGS})
# target_link_libraries(${PROJECT_NAME} PUBLIC adaSurfSizing_new)
target_link_libraries(${PROJECT_NAME} PUBLIC Autogrid_Triangle)
target_link_libraries(${PROJECT_NAME} PUBLIC REMesh_Triangle)

install(TARGETS ${PROJECT_NAME}
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib
RUNTIME DESTINATION bin
)