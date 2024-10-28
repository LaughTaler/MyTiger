# 预处理柯南
set(CON "conan")
execute_process(
    COMMAND ${CON} --version
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/cmake
    RESULT_VARIABLE CONAN_CHECK_RESULT
)
if (NOT CONAN_CHECK_RESULT EQUAL 0)
    set(CON "/home/gitlab-runner/anaconda3/envs/conan/bin/conan")
    execute_process(
        COMMAND ${CON} --version
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/cmake
        RESULT_VARIABLE CONAN_CHECK_RESULT
    )
endif()




if (NOT CONAN_CHECK_RESULT EQUAL 0)
    message(FATAL_ERROR "Conan is not installed or not found in PATH")
endif()



execute_process(
    COMMAND ${CON} remote add lab http://10.12.220.129:18081/artifactory/api/conan/conan-local --force
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/cmake
)


execute_process(
    COMMAND ${CON} remote login lab admin -p Test1234888
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/cmake
)

execute_process(
    COMMAND ${CON} profile detect
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/cmake
)


if (CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(platform "win64")
elseif (CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(platform "linux")
else ()
    set(platform "linux")
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(build_type "debug")
elseif (CMAKE_BUILD_TYPE STREQUAL "Release")
    set(build_type "release")
else ()
    set(build_type "debug")
endif()


execute_process(
    COMMAND ${CON} install . --profile=conan_profile_${platform}_${build_type} --remote=lab --output-folder=${CMAKE_BINARY_DIR} --build=missing
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/cmake
)



