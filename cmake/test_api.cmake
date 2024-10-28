message("building tigerdemo for tiger 1.5")
project(test_api)

set(EXTERNAL_DIR "${CMAKE_SOURCE_DIR}/extern")


# set(OCC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/extern/OpenCasCade)
set(Tran2fli_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/trans2_geo)
set(SOURCE 
${CMAKE_SOURCE_DIR}/test/test_api.cpp
${CMAKE_SOURCE_DIR}/test/test_api.h
)
set(Tran2fli_SOURCE ${CMAKE_SOURCE_DIR}/src/trans2_geo)

if(OPENMP_FOUND)
    message("OpenMP found")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


find_package(OpenCasCade REQUIRED)
# 这行加不加存疑
include(${CMAKE_SOURCE_DIR}/cmake/OpenCasCadeConfig.cmake)


include_directories(
  ${CMAKE_SOURCE_DIR}/include
 # ${OCC_DIR}/include
 ${opencascade_INCLUDE_DIRS}
  ${Tran2fli_SOURCE}
)

message("${opencascade_INCLUDE_DIRS}")
message("${opencascade_LIBRARIES_DEBUG}")


# 这个可能要改成cpp
add_executable(${PROJECT_NAME} ${SOURCE})


if(MSVC)

    target_link_directories(${PROJECT_NAME} PUBLIC ${opencascade_LIBRARIES_DEBUG})
    target_link_libraries(${PROJECT_NAME} PUBLIC ${opencascade_LIBRARIES_DEBUG})
    target_link_directories(${PROJECT_NAME} PUBLIC ${opencascade_INCLUDE_DIRS}/../lib)
    target_link_libraries(${PROJECT_NAME} PUBLIC 
        -Wl,--start-group 
        ${OCC_LIBS}
        -Wl,--end-group)
        
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
    # find_package(OpenCasCade REQUIRED)
    # #target_include_directories(${PROJECT_NAME}_test PUBLIC ${opencascade_INCLUDE_DIRS})
    # target_include_directories(${PROJECT_NAME} PUBLIC ${opencascade_INCLUDE_DIRS})

    # target_link_libraries(${PROJECT_NAME} PUBLIC pthread)
    # target_link_libraries(${PROJECT_NAME} PUBLIC dl)
    # target_link_directories(${PROJECT_NAME} PUBLIC ${opencascade_INCLUDE_DIRS}/../lib)
    # target_link_libraries(${PROJECT_NAME} PUBLIC 
    #     -Wl,--start-group 
    #     ${OCC_LIBS}
    #     -Wl,--end-group)
    target_include_directories(${PROJECT_NAME} PUBLIC ${opencascade_INCLUDE_DIRS})
    target_link_directories(${PROJECT_NAME} PUBLIC ${opencascade_INCLUDE_DIRS}/../lib)
    target_link_libraries(${PROJECT_NAME} PUBLIC 
        -Wl,--start-group 
        ${OCC_LIBS}
        -Wl,--end-group)
endif()


    #target_link_libraries(
    #    ${PROJECT_NAME}
    #    PUBLIC
    #    -Wl,--start-group
    #    ${OCC_LIBS}
    #    -Wl,--end-group)


# target_link_libraries(${PROJECT_NAME} PUBLIC ${PROJECT_NAME})


# target_link_libraries(${PROJECT_NAME} PUBLIC DT_Tetra)
target_link_libraries(${PROJECT_NAME} PUBLIC TRANS2_Geo)
target_link_libraries(${PROJECT_NAME} PUBLIC AFT_Tri)
target_link_libraries(${PROJECT_NAME} PUBLIC Quality_Data)
target_link_libraries(${PROJECT_NAME} PUBLIC SizingFunc)
target_link_libraries(${PROJECT_NAME} PUBLIC ${OpenMP_CXX_FLAGS})
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/catch2)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/simpleini)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/spdlog/include)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/cli11)