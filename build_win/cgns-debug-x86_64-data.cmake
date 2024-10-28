########### AGGREGATED COMPONENTS AND DEPENDENCIES FOR THE MULTI CONFIG #####################
#############################################################################################

set(cgns_COMPONENT_NAMES "")
if(DEFINED cgns_FIND_DEPENDENCY_NAMES)
  list(APPEND cgns_FIND_DEPENDENCY_NAMES )
  list(REMOVE_DUPLICATES cgns_FIND_DEPENDENCY_NAMES)
else()
  set(cgns_FIND_DEPENDENCY_NAMES )
endif()

########### VARIABLES #######################################################################
#############################################################################################
set(cgns_PACKAGE_FOLDER_DEBUG "C:/Users/25439/.conan2/p/cgns0dc58f4d203dd/p")
set(cgns_BUILD_MODULES_PATHS_DEBUG )


set(cgns_INCLUDE_DIRS_DEBUG "${cgns_PACKAGE_FOLDER_DEBUG}/include")
set(cgns_RES_DIRS_DEBUG )
set(cgns_DEFINITIONS_DEBUG )
set(cgns_SHARED_LINK_FLAGS_DEBUG )
set(cgns_EXE_LINK_FLAGS_DEBUG )
set(cgns_OBJECTS_DEBUG )
set(cgns_COMPILE_DEFINITIONS_DEBUG )
set(cgns_COMPILE_OPTIONS_C_DEBUG )
set(cgns_COMPILE_OPTIONS_CXX_DEBUG )
set(cgns_LIB_DIRS_DEBUG "${cgns_PACKAGE_FOLDER_DEBUG}/lib")
set(cgns_BIN_DIRS_DEBUG "${cgns_PACKAGE_FOLDER_DEBUG}/bin")
set(cgns_LIBRARY_TYPE_DEBUG SHARED)
set(cgns_IS_HOST_WINDOWS_DEBUG 1)
set(cgns_LIBS_DEBUG cgns libhdf5 libszip libzlib)
set(cgns_SYSTEM_LIBS_DEBUG )
set(cgns_FRAMEWORK_DIRS_DEBUG )
set(cgns_FRAMEWORKS_DEBUG )
set(cgns_BUILD_DIRS_DEBUG )
set(cgns_NO_SONAME_MODE_DEBUG FALSE)


# COMPOUND VARIABLES
set(cgns_COMPILE_OPTIONS_DEBUG
    "$<$<COMPILE_LANGUAGE:CXX>:${cgns_COMPILE_OPTIONS_CXX_DEBUG}>"
    "$<$<COMPILE_LANGUAGE:C>:${cgns_COMPILE_OPTIONS_C_DEBUG}>")
set(cgns_LINKER_FLAGS_DEBUG
    "$<$<STREQUAL:$<TARGET_PROPERTY:TYPE>,SHARED_LIBRARY>:${cgns_SHARED_LINK_FLAGS_DEBUG}>"
    "$<$<STREQUAL:$<TARGET_PROPERTY:TYPE>,MODULE_LIBRARY>:${cgns_SHARED_LINK_FLAGS_DEBUG}>"
    "$<$<STREQUAL:$<TARGET_PROPERTY:TYPE>,EXECUTABLE>:${cgns_EXE_LINK_FLAGS_DEBUG}>")


set(cgns_COMPONENTS_DEBUG )