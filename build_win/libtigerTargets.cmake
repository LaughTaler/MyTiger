# Load the debug and release variables
file(GLOB DATA_FILES "${CMAKE_CURRENT_LIST_DIR}/libtiger-*-data.cmake")

foreach(f ${DATA_FILES})
    include(${f})
endforeach()

# Create the targets for all the components
foreach(_COMPONENT ${libtiger_COMPONENT_NAMES} )
    if(NOT TARGET ${_COMPONENT})
        add_library(${_COMPONENT} INTERFACE IMPORTED)
        message(${libtiger_MESSAGE_MODE} "Conan: Component target declared '${_COMPONENT}'")
    endif()
endforeach()

if(NOT TARGET libtiger::libtiger)
    add_library(libtiger::libtiger INTERFACE IMPORTED)
    message(${libtiger_MESSAGE_MODE} "Conan: Target declared 'libtiger::libtiger'")
endif()
# Load the debug and release library finders
file(GLOB CONFIG_FILES "${CMAKE_CURRENT_LIST_DIR}/libtiger-Target-*.cmake")

foreach(f ${CONFIG_FILES})
    include(${f})
endforeach()