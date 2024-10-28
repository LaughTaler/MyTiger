execute_process(COMMAND ${CMAKE_SOURCE_DIR}/extern/geogram/run.sh)
execute_process(COMMAND cp ${CMAKE_SOURCE_DIR}/extern/geogram/build/Linux64-gcc-Debug/lib/libgeogram.a ${CMAKE_SOURCE_DIR}/lib/libgeogram_d.a)
execute_process(COMMAND cp ${CMAKE_SOURCE_DIR}/extern/geogram/build/Linux64-gcc-Release/lib/libgeogram.a ${CMAKE_SOURCE_DIR}/lib/libgeogram_r.a)



