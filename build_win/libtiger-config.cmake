########## MACROS ###########################################################################
#############################################################################################

# Requires CMake > 3.15
if(${CMAKE_VERSION} VERSION_LESS "3.15")
    message(FATAL_ERROR "The 'CMakeDeps' generator only works with CMake >= 3.15")
endif()

if(libtiger_FIND_QUIETLY)
    set(libtiger_MESSAGE_MODE VERBOSE)
else()
    set(libtiger_MESSAGE_MODE STATUS)
endif()

include(${CMAKE_CURRENT_LIST_DIR}/cmakedeps_macros.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/libtigerTargets.cmake)
include(CMakeFindDependencyMacro)

check_build_type_defined()

foreach(_DEPENDENCY ${libtiger_FIND_DEPENDENCY_NAMES} )
    # Check that we have not already called a find_package with the transitive dependency
    if(NOT ${_DEPENDENCY}_FOUND)
        find_dependency(${_DEPENDENCY} REQUIRED ${${_DEPENDENCY}_FIND_MODE})
    endif()
endforeach()

set(libtiger_VERSION_STRING "10.23")
set(libtiger_INCLUDE_DIRS ${libtiger_INCLUDE_DIRS_DEBUG} )
set(libtiger_INCLUDE_DIR ${libtiger_INCLUDE_DIRS_DEBUG} )
set(libtiger_LIBRARIES ${libtiger_LIBRARIES_DEBUG} )
set(libtiger_DEFINITIONS ${libtiger_DEFINITIONS_DEBUG} )


# Only the last installed configuration BUILD_MODULES are included to avoid the collision
foreach(_BUILD_MODULE ${libtiger_BUILD_MODULES_PATHS_DEBUG} )
    message(${libtiger_MESSAGE_MODE} "Conan: Including build module from '${_BUILD_MODULE}'")
    include(${_BUILD_MODULE})
endforeach()

