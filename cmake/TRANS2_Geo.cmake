message("building test example for tiger 1.5")
project(TRANS2_Geo)

set(EXTERNAL_DIR "${CMAKE_SOURCE_DIR}/extern")

set(Tran2fli_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/trans2_geo)

file(GLOB SOURCE "${Tran2fli_DIR}/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/src/module/TRANS2_Geo/*.cpp" 
# "${Tran2fli_DIR}/wyfdt2d/src/*.cpp"
#  "${Tran2fli_DIR}/wyfdt2d/src/*.h"
#   "${Tran2fli_DIR}/wyfdt2d/extern/geom/*.h"
#   "${Tran2fli_DIR}/wyfdt2d/extern/geom/*.cpp"
#   "${Tran2fli_DIR}/wyfdt2d/extern/geom/*.c"
 )

find_package(opencascade REQUIRED)

set(Tran2fli_SOURCE ${CMAKE_SOURCE_DIR}/src/trans2_geo)


include(${CMAKE_SOURCE_DIR}/cmake/OpenCasCadeConfig.cmake)
# include_directories (
#   ${Tran2fli_DIR}/wyfdt2d/src
#   ${Tran2fli_DIR}/wyfdt2d/extern/geom
#   ${Tran2fli_DIR}/wyfdt2d/extern/cli11
#   ${Tran2fli_DIR}/wyfdt2d/extern/spdlog/include
# )

include_directories(
  ${opencascade_INCLUDE_DIRS}
  ${Tran2fli_SOURCE}
)
message("${opencascade_INCLUDE_DIRS}")
add_library(${PROJECT_NAME} SHARED ${Tran2fli_SOURCE} ${SOURCE})
add_executable(${PROJECT_NAME}_test ${CMAKE_SOURCE_DIR}/test/test_TRANS2_Geo.cpp )

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


      #  get_target_property(result ${opencascade_LIBRARIES}  INTERFACE_LINK_LIBRARIES)

    #target_link_libraries(
    #    ${PROJECT_NAME}
    #    PUBLIC
    #    -Wl,--start-group
    #    ${OCC_LIBS}
    #    -Wl,--end-group)
target_link_libraries(${PROJECT_NAME} PUBLIC DT_Tri)
target_link_libraries(${PROJECT_NAME}_test PUBLIC ${PROJECT_NAME})
target_include_directories(${PROJECT_NAME}_test PUBLIC ${CMAKE_SOURCE_DIR}/extern/catch2)
