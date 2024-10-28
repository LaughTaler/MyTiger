########## MACROS ###########################################################################
#############################################################################################

# Requires CMake > 3.15
if(${CMAKE_VERSION} VERSION_LESS "3.15")
    message(FATAL_ERROR "The 'CMakeDeps' generator only works with CMake >= 3.15")
endif()

if(opencascade_FIND_QUIETLY)
    set(opencascade_MESSAGE_MODE VERBOSE)
else()
    set(opencascade_MESSAGE_MODE STATUS)
endif()

include(${CMAKE_CURRENT_LIST_DIR}/cmakedeps_macros.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/opencascadeTargets.cmake)
include(CMakeFindDependencyMacro)

check_build_type_defined()

foreach(_DEPENDENCY ${opencascade_FIND_DEPENDENCY_NAMES} )
    # Check that we have not already called a find_package with the transitive dependency
    if(NOT ${_DEPENDENCY}_FOUND)
        find_dependency(${_DEPENDENCY} REQUIRED ${${_DEPENDENCY}_FIND_MODE})
    endif()
endforeach()

set(opencascade_VERSION_STRING "7.7.0")
set(opencascade_INCLUDE_DIRS ${opencascade_INCLUDE_DIRS_DEBUG} )
set(opencascade_INCLUDE_DIR ${opencascade_INCLUDE_DIRS_DEBUG} )
set(opencascade_LIBRARIES ${opencascade_LIBRARIES_DEBUG} )
set(opencascade_DEFINITIONS ${opencascade_DEFINITIONS_DEBUG} )


# Only the last installed configuration BUILD_MODULES are included to avoid the collision
foreach(_BUILD_MODULE ${opencascade_BUILD_MODULES_PATHS_DEBUG} )
    message(${opencascade_MESSAGE_MODE} "Conan: Including build module from '${_BUILD_MODULE}'")
    include(${_BUILD_MODULE})
endforeach()


