message("building test example for tiger 1.5")
project(SizingFunc)

set(EXTERNAL_DIR "${CMAKE_SOURCE_DIR}/extern")

set(Sizefunction_SOURCE ${CMAKE_SOURCE_DIR}/src/SizingFunc)
include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/include)
set(DIR ${CMAKE_SOURCE_DIR}/src/module/SizingFunc)
set(SOURCE
    ${DIR}/sizing_field.cpp
)
add_library(${PROJECT_NAME} SHARED ${SOURCE})
target_link_libraries(${PROJECT_NAME} PUBLIC adaSurfSizing_new)
target_link_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/build/lib/Release) 
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/eigen)


add_executable(${PROJECT_NAME}_test ${CMAKE_SOURCE_DIR}/test/test_${PROJECT_NAME}.cpp)
target_link_libraries(${PROJECT_NAME}_test PUBLIC ${PROJECT_NAME})
target_include_directories(${PROJECT_NAME}_test PUBLIC ${CMAKE_SOURCE_DIR}/extern/catch2)
target_include_directories(${PROJECT_NAME}_test PUBLIC ${CMAKE_SOURCE_DIR}/extern/eigen)
