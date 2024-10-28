function(add_internal_module name_internal repository commit_id)
    message("try to add moudle ${name_internal}")
    set(DIR_MODULE_DIR "${DIR_PATH}/${name_internal}")
    set(CMAKELISTS_PATH "${DIR_MODULE_DIR}/CMakeLists.txt")
    set(CMAKECONFIG_PATH "${CMAKE_SOURCE_DIR}/cmake/${name_internal}.cmake")

    if(EXISTS ${CMAKELISTS_PATH})


    execute_process(
        COMMAND git fetch origin master
        WORKING_DIRECTORY ${DIR_MODULE_DIR}
    )

    execute_process(
        COMMAND git merge origin/master
        OUTPUT_VARIABLE REMOTE_HEAD
        WORKING_DIRECTORY ${DIR_MODULE_DIR}
    )

    execute_process(
        COMMAND git checkout ${commit_id}
        WORKING_DIRECTORY ${DIR_MODULE_DIR}
    )
        # execute_process(
        #     COMMAND git fetch origin master
        #     WORKING_DIRECTORY ${DIR_MODULE_DIR}
        # )

        # execute_process(
        #     COMMAND git rev-parse HEAD
        #     OUTPUT_VARIABLE LOCAL_HEAD
        #     WORKING_DIRECTORY ${DIR_MODULE_DIR}
        # )

        # execute_process(
        #     COMMAND git rev-parse origin/master
        #     OUTPUT_VARIABLE REMOTE_HEAD
        #     WORKING_DIRECTORY ${DIR_MODULE_DIR}
        # )

        # if(NOT "${LOCAL_HEAD}" STREQUAL "${REMOTE_HEAD}")
        # execute_process(
        #     COMMAND git merge origin/master
        #     OUTPUT_VARIABLE REMOTE_HEAD
        #     WORKING_DIRECTORY ${DIR_MODULE_DIR}
        # )
        #     message(WARNING "Remote repository is not up to date!")
        # else()
        #     message(STATUS "Remote repository is up to date.")
        # endif()

        include(${CMAKECONFIG_PATH})
        add_subdirectory(${DIR_MODULE_DIR})
        message(STATUS "Directory ${CMAKELISTS_PATH} added.")
    else()
        message(STATUS "try to download ${name} from ${repository}")
        
        execute_process(
            COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_SOURCE_DIR}/src/${name_internal}
            RESULT_VARIABLE result
        )
        message("${CMAKE_SOURCE_DIR}/src/${name_internal}")
        execute_process(
            COMMAND git clone ${repository} ${name_internal}
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src
        )
        execute_process(
            COMMAND git checkout ${commit_id}
            WORKING_DIRECTORY ${DIR_MODULE_DIR}
        )

        message(STATUS "try to download ${name_internal} from ${repository}")

        if(EXISTS ${CMAKELISTS_PATH})
            include(${CMAKECONFIG_PATH})
            add_subdirectory(${DIR_MODULE_DIR})
            install(TARGETS ${name_internal}
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib
            RUNTIME DESTINATION bin
        )
           install(TARGETS ${name_internal}_test
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib
            RUNTIME DESTINATION bin
        )

        else()
            message(STATUS "download failed.")
        endif()
    endif()



endfunction()