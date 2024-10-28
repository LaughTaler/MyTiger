# Avoid multiple calls to find_package to append duplicated properties to the targets
include_guard()########### VARIABLES #######################################################################
#############################################################################################
set(libtiger_FRAMEWORKS_FOUND_DEBUG "") # Will be filled later
conan_find_apple_frameworks(libtiger_FRAMEWORKS_FOUND_DEBUG "${libtiger_FRAMEWORKS_DEBUG}" "${libtiger_FRAMEWORK_DIRS_DEBUG}")

set(libtiger_LIBRARIES_TARGETS "") # Will be filled later


######## Create an interface target to contain all the dependencies (frameworks, system and conan deps)
if(NOT TARGET libtiger_DEPS_TARGET)
    add_library(libtiger_DEPS_TARGET INTERFACE IMPORTED)
endif()

set_property(TARGET libtiger_DEPS_TARGET
             APPEND PROPERTY INTERFACE_LINK_LIBRARIES
             $<$<CONFIG:Debug>:${libtiger_FRAMEWORKS_FOUND_DEBUG}>
             $<$<CONFIG:Debug>:${libtiger_SYSTEM_LIBS_DEBUG}>
             $<$<CONFIG:Debug>:>)

####### Find the libraries declared in cpp_info.libs, create an IMPORTED target for each one and link the
####### libtiger_DEPS_TARGET to all of them
conan_package_library_targets("${libtiger_LIBS_DEBUG}"    # libraries
                              "${libtiger_LIB_DIRS_DEBUG}" # package_libdir
                              "${libtiger_BIN_DIRS_DEBUG}" # package_bindir
                              "${libtiger_LIBRARY_TYPE_DEBUG}"
                              "${libtiger_IS_HOST_WINDOWS_DEBUG}"
                              libtiger_DEPS_TARGET
                              libtiger_LIBRARIES_TARGETS  # out_libraries_targets
                              "_DEBUG"
                              "libtiger"    # package_name
                              "${libtiger_NO_SONAME_MODE_DEBUG}")  # soname

# FIXME: What is the result of this for multi-config? All configs adding themselves to path?
set(CMAKE_MODULE_PATH ${libtiger_BUILD_DIRS_DEBUG} ${CMAKE_MODULE_PATH})

########## GLOBAL TARGET PROPERTIES Debug ########################################
    set_property(TARGET libtiger::libtiger
                 APPEND PROPERTY INTERFACE_LINK_LIBRARIES
                 $<$<CONFIG:Debug>:${libtiger_OBJECTS_DEBUG}>
                 $<$<CONFIG:Debug>:${libtiger_LIBRARIES_TARGETS}>
                 )

    if("${libtiger_LIBS_DEBUG}" STREQUAL "")
        # If the package is not declaring any "cpp_info.libs" the package deps, system libs,
        # frameworks etc are not linked to the imported targets and we need to do it to the
        # global target
        set_property(TARGET libtiger::libtiger
                     APPEND PROPERTY INTERFACE_LINK_LIBRARIES
                     libtiger_DEPS_TARGET)
    endif()

    set_property(TARGET libtiger::libtiger
                 APPEND PROPERTY INTERFACE_LINK_OPTIONS
                 $<$<CONFIG:Debug>:${libtiger_LINKER_FLAGS_DEBUG}>)
    set_property(TARGET libtiger::libtiger
                 APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                 $<$<CONFIG:Debug>:${libtiger_INCLUDE_DIRS_DEBUG}>)
    # Necessary to find LINK shared libraries in Linux
    set_property(TARGET libtiger::libtiger
                 APPEND PROPERTY INTERFACE_LINK_DIRECTORIES
                 $<$<CONFIG:Debug>:${libtiger_LIB_DIRS_DEBUG}>)
    set_property(TARGET libtiger::libtiger
                 APPEND PROPERTY INTERFACE_COMPILE_DEFINITIONS
                 $<$<CONFIG:Debug>:${libtiger_COMPILE_DEFINITIONS_DEBUG}>)
    set_property(TARGET libtiger::libtiger
                 APPEND PROPERTY INTERFACE_COMPILE_OPTIONS
                 $<$<CONFIG:Debug>:${libtiger_COMPILE_OPTIONS_DEBUG}>)

########## For the modules (FindXXX)
set(libtiger_LIBRARIES_DEBUG libtiger::libtiger)
