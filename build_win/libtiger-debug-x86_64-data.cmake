########### AGGREGATED COMPONENTS AND DEPENDENCIES FOR THE MULTI CONFIG #####################
#############################################################################################

set(libtiger_COMPONENT_NAMES "")
if(DEFINED libtiger_FIND_DEPENDENCY_NAMES)
  list(APPEND libtiger_FIND_DEPENDENCY_NAMES )
  list(REMOVE_DUPLICATES libtiger_FIND_DEPENDENCY_NAMES)
else()
  set(libtiger_FIND_DEPENDENCY_NAMES )
endif()

########### VARIABLES #######################################################################
#############################################################################################
set(libtiger_PACKAGE_FOLDER_DEBUG "C:/Users/25439/.conan2/p/libti5f9e0ffcf824a/p")
set(libtiger_BUILD_MODULES_PATHS_DEBUG )


set(libtiger_INCLUDE_DIRS_DEBUG "${libtiger_PACKAGE_FOLDER_DEBUG}/include")
set(libtiger_RES_DIRS_DEBUG )
set(libtiger_DEFINITIONS_DEBUG )
set(libtiger_SHARED_LINK_FLAGS_DEBUG )
set(libtiger_EXE_LINK_FLAGS_DEBUG )
set(libtiger_OBJECTS_DEBUG )
set(libtiger_COMPILE_DEFINITIONS_DEBUG )
set(libtiger_COMPILE_OPTIONS_C_DEBUG )
set(libtiger_COMPILE_OPTIONS_CXX_DEBUG )
set(libtiger_LIB_DIRS_DEBUG "${libtiger_PACKAGE_FOLDER_DEBUG}/lib")
set(libtiger_BIN_DIRS_DEBUG "${libtiger_PACKAGE_FOLDER_DEBUG}/bin")
set(libtiger_LIBRARY_TYPE_DEBUG SHARED)
set(libtiger_IS_HOST_WINDOWS_DEBUG 1)
set(libtiger_LIBS_DEBUG )
set(libtiger_SYSTEM_LIBS_DEBUG )
set(libtiger_FRAMEWORK_DIRS_DEBUG )
set(libtiger_FRAMEWORKS_DEBUG )
set(libtiger_BUILD_DIRS_DEBUG )
set(libtiger_NO_SONAME_MODE_DEBUG FALSE)


# COMPOUND VARIABLES
set(libtiger_COMPILE_OPTIONS_DEBUG
    "$<$<COMPILE_LANGUAGE:CXX>:${libtiger_COMPILE_OPTIONS_CXX_DEBUG}>"
    "$<$<COMPILE_LANGUAGE:C>:${libtiger_COMPILE_OPTIONS_C_DEBUG}>")
set(libtiger_LINKER_FLAGS_DEBUG
    "$<$<STREQUAL:$<TARGET_PROPERTY:TYPE>,SHARED_LIBRARY>:${libtiger_SHARED_LINK_FLAGS_DEBUG}>"
    "$<$<STREQUAL:$<TARGET_PROPERTY:TYPE>,MODULE_LIBRARY>:${libtiger_SHARED_LINK_FLAGS_DEBUG}>"
    "$<$<STREQUAL:$<TARGET_PROPERTY:TYPE>,EXECUTABLE>:${libtiger_EXE_LINK_FLAGS_DEBUG}>")


set(libtiger_COMPONENTS_DEBUG )