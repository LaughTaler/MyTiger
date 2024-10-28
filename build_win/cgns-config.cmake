########## MACROS ###########################################################################
#############################################################################################

# Requires CMake > 3.15
if(${CMAKE_VERSION} VERSION_LESS "3.15")
    message(FATAL_ERROR "The 'CMakeDeps' generator only works with CMake >= 3.15")
endif()

if(cgns_FIND_QUIETLY)
    set(cgns_MESSAGE_MODE VERBOSE)
else()
    set(cgns_MESSAGE_MODE STATUS)
endif()

include(${CMAKE_CURRENT_LIST_DIR}/cmakedeps_macros.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/cgnsTargets.cmake)
include(CMakeFindDependencyMacro)

check_build_type_defined()

foreach(_DEPENDENCY ${cgns_FIND_DEPENDENCY_NAMES} )
    # Check that we have not already called a find_package with the transitive dependency
    if(NOT ${_DEPENDENCY}_FOUND)
        find_dependency(${_DEPENDENCY} REQUIRED ${${_DEPENDENCY}_FIND_MODE})
    endif()
endforeach()

set(cgns_VERSION_STRING "3.3")
set(cgns_INCLUDE_DIRS ${cgns_INCLUDE_DIRS_DEBUG} )
set(cgns_INCLUDE_DIR ${cgns_INCLUDE_DIRS_DEBUG} )
set(cgns_LIBRARIES ${cgns_LIBRARIES_DEBUG} )
set(cgns_DEFINITIONS ${cgns_DEFINITIONS_DEBUG} )


# Only the last installed configuration BUILD_MODULES are included to avoid the collision
foreach(_BUILD_MODULE ${cgns_BUILD_MODULES_PATHS_DEBUG} )
    message(${cgns_MESSAGE_MODE} "Conan: Including build module from '${_BUILD_MODULE}'")
    include(${_BUILD_MODULE})
endforeach()


