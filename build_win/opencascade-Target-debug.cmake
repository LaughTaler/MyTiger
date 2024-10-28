# Avoid multiple calls to find_package to append duplicated properties to the targets
include_guard()########### VARIABLES #######################################################################
#############################################################################################
set(opencascade_FRAMEWORKS_FOUND_DEBUG "") # Will be filled later
conan_find_apple_frameworks(opencascade_FRAMEWORKS_FOUND_DEBUG "${opencascade_FRAMEWORKS_DEBUG}" "${opencascade_FRAMEWORK_DIRS_DEBUG}")

set(opencascade_LIBRARIES_TARGETS "") # Will be filled later


######## Create an interface target to contain all the dependencies (frameworks, system and conan deps)
if(NOT TARGET opencascade_DEPS_TARGET)
    add_library(opencascade_DEPS_TARGET INTERFACE IMPORTED)
endif()

set_property(TARGET opencascade_DEPS_TARGET
             APPEND PROPERTY INTERFACE_LINK_LIBRARIES
             $<$<CONFIG:Debug>:${opencascade_FRAMEWORKS_FOUND_DEBUG}>
             $<$<CONFIG:Debug>:${opencascade_SYSTEM_LIBS_DEBUG}>
             $<$<CONFIG:Debug>:>)

####### Find the libraries declared in cpp_info.libs, create an IMPORTED target for each one and link the
####### opencascade_DEPS_TARGET to all of them
conan_package_library_targets("${opencascade_LIBS_DEBUG}"    # libraries
                              "${opencascade_LIB_DIRS_DEBUG}" # package_libdir
                              "${opencascade_BIN_DIRS_DEBUG}" # package_bindir
                              "${opencascade_LIBRARY_TYPE_DEBUG}"
                              "${opencascade_IS_HOST_WINDOWS_DEBUG}"
                              opencascade_DEPS_TARGET
                              opencascade_LIBRARIES_TARGETS  # out_libraries_targets
                              "_DEBUG"
                              "opencascade"    # package_name
                              "${opencascade_NO_SONAME_MODE_DEBUG}")  # soname

# FIXME: What is the result of this for multi-config? All configs adding themselves to path?
set(CMAKE_MODULE_PATH ${opencascade_BUILD_DIRS_DEBUG} ${CMAKE_MODULE_PATH})

########## GLOBAL TARGET PROPERTIES Debug ########################################
    set_property(TARGET opencascade::opencascade
                 APPEND PROPERTY INTERFACE_LINK_LIBRARIES
                 $<$<CONFIG:Debug>:${opencascade_OBJECTS_DEBUG}>
                 $<$<CONFIG:Debug>:${opencascade_LIBRARIES_TARGETS}>
                 )

    if("${opencascade_LIBS_DEBUG}" STREQUAL "")
        # If the package is not declaring any "cpp_info.libs" the package deps, system libs,
        # frameworks etc are not linked to the imported targets and we need to do it to the
        # global target
        set_property(TARGET opencascade::opencascade
                     APPEND PROPERTY INTERFACE_LINK_LIBRARIES
                     opencascade_DEPS_TARGET)
    endif()

    set_property(TARGET opencascade::opencascade
                 APPEND PROPERTY INTERFACE_LINK_OPTIONS
                 $<$<CONFIG:Debug>:${opencascade_LINKER_FLAGS_DEBUG}>)
    set_property(TARGET opencascade::opencascade
                 APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                 $<$<CONFIG:Debug>:${opencascade_INCLUDE_DIRS_DEBUG}>)
    # Necessary to find LINK shared libraries in Linux
    set_property(TARGET opencascade::opencascade
                 APPEND PROPERTY INTERFACE_LINK_DIRECTORIES
                 $<$<CONFIG:Debug>:${opencascade_LIB_DIRS_DEBUG}>)
    set_property(TARGET opencascade::opencascade
                 APPEND PROPERTY INTERFACE_COMPILE_DEFINITIONS
                 $<$<CONFIG:Debug>:${opencascade_COMPILE_DEFINITIONS_DEBUG}>)
    set_property(TARGET opencascade::opencascade
                 APPEND PROPERTY INTERFACE_COMPILE_OPTIONS
                 $<$<CONFIG:Debug>:${opencascade_COMPILE_OPTIONS_DEBUG}>)

########## For the modules (FindXXX)
set(opencascade_LIBRARIES_DEBUG opencascade::opencascade)
