########### AGGREGATED COMPONENTS AND DEPENDENCIES FOR THE MULTI CONFIG #####################
#############################################################################################

set(opencascade_COMPONENT_NAMES "")
if(DEFINED opencascade_FIND_DEPENDENCY_NAMES)
  list(APPEND opencascade_FIND_DEPENDENCY_NAMES )
  list(REMOVE_DUPLICATES opencascade_FIND_DEPENDENCY_NAMES)
else()
  set(opencascade_FIND_DEPENDENCY_NAMES )
endif()

########### VARIABLES #######################################################################
#############################################################################################
set(opencascade_PACKAGE_FOLDER_DEBUG "C:/Users/25439/.conan2/p/openc5557227e8d14d/p")
set(opencascade_BUILD_MODULES_PATHS_DEBUG )


set(opencascade_INCLUDE_DIRS_DEBUG "${opencascade_PACKAGE_FOLDER_DEBUG}/include")
set(opencascade_RES_DIRS_DEBUG )
set(opencascade_DEFINITIONS_DEBUG )
set(opencascade_SHARED_LINK_FLAGS_DEBUG )
set(opencascade_EXE_LINK_FLAGS_DEBUG )
set(opencascade_OBJECTS_DEBUG )
set(opencascade_COMPILE_DEFINITIONS_DEBUG )
set(opencascade_COMPILE_OPTIONS_C_DEBUG )
set(opencascade_COMPILE_OPTIONS_CXX_DEBUG )
set(opencascade_LIB_DIRS_DEBUG "${opencascade_PACKAGE_FOLDER_DEBUG}/lib")
set(opencascade_BIN_DIRS_DEBUG "${opencascade_PACKAGE_FOLDER_DEBUG}/bin")
set(opencascade_LIBRARY_TYPE_DEBUG SHARED)
set(opencascade_IS_HOST_WINDOWS_DEBUG 1)
set(opencascade_LIBS_DEBUG TKBO TKBRep TKBin TKBinL TKBinTObj TKBinXCAF TKBool TKCAF TKCDF TKDCAF TKDraw TKExpress TKFeat TKFillet TKG2d TKG3d TKGeomAlgo TKGeomBase TKHLR TKIGES TKIVtk TKIVtkDraw TKLCAF TKMath TKMesh TKMeshVS TKOffset TKOpenGl TKOpenGlTest TKPrim TKQADraw TKRWMesh TKSTEP TKSTEP209 TKSTEPAttr TKSTEPBase TKSTL TKService TKShHealing TKStd TKStdL TKTObj TKTObjDRAW TKTopAlgo TKTopTest TKV3d TKVCAF TKVRML TKViewerTest TKXCAF TKXDE TKXDECascade TKXDEDRAW TKXDEIGES TKXDESTEP TKXMesh TKXSBase TKXSDRAW TKXml TKXmlL TKXmlTObj TKXmlXCAF TKernel)
set(opencascade_SYSTEM_LIBS_DEBUG )
set(opencascade_FRAMEWORK_DIRS_DEBUG )
set(opencascade_FRAMEWORKS_DEBUG )
set(opencascade_BUILD_DIRS_DEBUG )
set(opencascade_NO_SONAME_MODE_DEBUG FALSE)


# COMPOUND VARIABLES
set(opencascade_COMPILE_OPTIONS_DEBUG
    "$<$<COMPILE_LANGUAGE:CXX>:${opencascade_COMPILE_OPTIONS_CXX_DEBUG}>"
    "$<$<COMPILE_LANGUAGE:C>:${opencascade_COMPILE_OPTIONS_C_DEBUG}>")
set(opencascade_LINKER_FLAGS_DEBUG
    "$<$<STREQUAL:$<TARGET_PROPERTY:TYPE>,SHARED_LIBRARY>:${opencascade_SHARED_LINK_FLAGS_DEBUG}>"
    "$<$<STREQUAL:$<TARGET_PROPERTY:TYPE>,MODULE_LIBRARY>:${opencascade_SHARED_LINK_FLAGS_DEBUG}>"
    "$<$<STREQUAL:$<TARGET_PROPERTY:TYPE>,EXECUTABLE>:${opencascade_EXE_LINK_FLAGS_DEBUG}>")


set(opencascade_COMPONENTS_DEBUG )