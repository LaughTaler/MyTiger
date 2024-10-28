# Load the debug and release variables
file(GLOB DATA_FILES "${CMAKE_CURRENT_LIST_DIR}/opencascade-*-data.cmake")

foreach(f ${DATA_FILES})
    include(${f})
endforeach()

# Create the targets for all the components
foreach(_COMPONENT ${opencascade_COMPONENT_NAMES} )
    if(NOT TARGET ${_COMPONENT})
        add_library(${_COMPONENT} INTERFACE IMPORTED)
        message(${opencascade_MESSAGE_MODE} "Conan: Component target declared '${_COMPONENT}'")
    endif()
endforeach()

if(NOT TARGET opencascade::opencascade)
    add_library(opencascade::opencascade INTERFACE IMPORTED)
    message(${opencascade_MESSAGE_MODE} "Conan: Target declared 'opencascade::opencascade'")
endif()
# Load the debug and release library finders
file(GLOB CONFIG_FILES "${CMAKE_CURRENT_LIST_DIR}/opencascade-Target-*.cmake")

foreach(f ${CONFIG_FILES})
    include(${f})
endforeach()