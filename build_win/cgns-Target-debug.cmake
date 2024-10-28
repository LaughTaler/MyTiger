# Avoid multiple calls to find_package to append duplicated properties to the targets
include_guard()########### VARIABLES #######################################################################
#############################################################################################
set(cgns_FRAMEWORKS_FOUND_DEBUG "") # Will be filled later
conan_find_apple_frameworks(cgns_FRAMEWORKS_FOUND_DEBUG "${cgns_FRAMEWORKS_DEBUG}" "${cgns_FRAMEWORK_DIRS_DEBUG}")

set(cgns_LIBRARIES_TARGETS "") # Will be filled later


######## Create an interface target to contain all the dependencies (frameworks, system and conan deps)
if(NOT TARGET cgns_DEPS_TARGET)
    add_library(cgns_DEPS_TARGET INTERFACE IMPORTED)
endif()

set_property(TARGET cgns_DEPS_TARGET
             APPEND PROPERTY INTERFACE_LINK_LIBRARIES
             $<$<CONFIG:Debug>:${cgns_FRAMEWORKS_FOUND_DEBUG}>
             $<$<CONFIG:Debug>:${cgns_SYSTEM_LIBS_DEBUG}>
             $<$<CONFIG:Debug>:>)

####### Find the libraries declared in cpp_info.libs, create an IMPORTED target for each one and link the
####### cgns_DEPS_TARGET to all of them
conan_package_library_targets("${cgns_LIBS_DEBUG}"    # libraries
                              "${cgns_LIB_DIRS_DEBUG}" # package_libdir
                              "${cgns_BIN_DIRS_DEBUG}" # package_bindir
                              "${cgns_LIBRARY_TYPE_DEBUG}"
                              "${cgns_IS_HOST_WINDOWS_DEBUG}"
                              cgns_DEPS_TARGET
                              cgns_LIBRARIES_TARGETS  # out_libraries_targets
                              "_DEBUG"
                              "cgns"    # package_name
                              "${cgns_NO_SONAME_MODE_DEBUG}")  # soname

# FIXME: What is the result of this for multi-config? All configs adding themselves to path?
set(CMAKE_MODULE_PATH ${cgns_BUILD_DIRS_DEBUG} ${CMAKE_MODULE_PATH})

########## GLOBAL TARGET PROPERTIES Debug ########################################
    set_property(TARGET cgns::cgns
                 APPEND PROPERTY INTERFACE_LINK_LIBRARIES
                 $<$<CONFIG:Debug>:${cgns_OBJECTS_DEBUG}>
                 $<$<CONFIG:Debug>:${cgns_LIBRARIES_TARGETS}>
                 )

    if("${cgns_LIBS_DEBUG}" STREQUAL "")
        # If the package is not declaring any "cpp_info.libs" the package deps, system libs,
        # frameworks etc are not linked to the imported targets and we need to do it to the
        # global target
        set_property(TARGET cgns::cgns
                     APPEND PROPERTY INTERFACE_LINK_LIBRARIES
                     cgns_DEPS_TARGET)
    endif()

    set_property(TARGET cgns::cgns
                 APPEND PROPERTY INTERFACE_LINK_OPTIONS
                 $<$<CONFIG:Debug>:${cgns_LINKER_FLAGS_DEBUG}>)
    set_property(TARGET cgns::cgns
                 APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                 $<$<CONFIG:Debug>:${cgns_INCLUDE_DIRS_DEBUG}>)
    # Necessary to find LINK shared libraries in Linux
    set_property(TARGET cgns::cgns
                 APPEND PROPERTY INTERFACE_LINK_DIRECTORIES
                 $<$<CONFIG:Debug>:${cgns_LIB_DIRS_DEBUG}>)
    set_property(TARGET cgns::cgns
                 APPEND PROPERTY INTERFACE_COMPILE_DEFINITIONS
                 $<$<CONFIG:Debug>:${cgns_COMPILE_DEFINITIONS_DEBUG}>)
    set_property(TARGET cgns::cgns
                 APPEND PROPERTY INTERFACE_COMPILE_OPTIONS
                 $<$<CONFIG:Debug>:${cgns_COMPILE_OPTIONS_DEBUG}>)

########## For the modules (FindXXX)
set(cgns_LIBRARIES_DEBUG cgns::cgns)
