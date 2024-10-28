set(EXDIR_PATH "${CMAKE_SOURCE_DIR}/extern")

function(add_extern_submodule name_internal repository)
    message("try to add extern ${name_internal}")
    set(DIR_MODULE_DIR "${EXDIR_PATH}/${name_internal}")
    set(CMAKELISTS_PATH "${DIR_MODULE_DIR}/CMakeLists.txt")

    if(EXISTS ${CMAKELISTS_PATH})
        #add_subdirectory(${DIR_MODULE_DIR})
        message(STATUS "Directory ${CMAKELISTS_PATH} added.")
    else()
 
        # download_project(PROJ     ${name_internal}
        #     GIT_REPOSITORY      ${repository}
        #     GIT_TAG             master
        #     UPDATE_DISCONNECTED 1
        #     PREFIX             ""
        #    # DOWNLOAD_DIR         ${CMAKE_SOURCE_DIR}/extern/${name_internal}
        #     SOURCE_DIR         ${CMAKE_SOURCE_DIR}/extern/${name_internal}
        #     QUIET
        # )

        execute_process(
            COMMAND git clone ${repository} ${name_internal}
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/extern
        )

        message(STATUS "try to download ${name_internal} form ${repository}.")

        if(EXISTS ${CMAKELISTS_PATH})
            #include(${CMAKECONFIG_PATH})
        else()
            message(STATUS "downloading failed.")
            message(STATUS "Directory ${CMAKELISTS_PATH} does not exist, omitted.")
        endif()
        
    endif()
endfunction()


