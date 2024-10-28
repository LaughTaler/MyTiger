include("${CMAKE_SOURCE_DIR}/cmake/DownloadProject.cmake")
set(DIR_PATH "${CMAKE_SOURCE_DIR}/src")

function(add_module name name_internal repository commit_id)
    message("try to add moudle ${name_internal}")
    set(DIR_MODULE_DIR "${DIR_PATH}/${name_internal}")
    set(CMAKELISTS_PATH "${DIR_MODULE_DIR}/CMakeLists.txt")
    set(CMAKECONFIG_PATH "${CMAKE_SOURCE_DIR}/cmake/${name}.cmake")

    # file(MAKE_DIRECTORY ${DIR_PATH})
    # file(MAKE_DIRECTORY ${DIR_MODULE_DIR})

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
        if(MSVC)
                add_compile_options(/W0)
                add_compile_options(/Z7)
                set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Z7")
            endif()
        add_subdirectory(${DIR_MODULE_DIR})
        message(STATUS "Directory ${CMAKELISTS_PATH} added.")
    else()
        message(STATUS "try to download ${name} from ${repository}")
        
        execute_process(
            COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_SOURCE_DIR}/src/${name_internal}
            RESULT_VARIABLE result
        )
        message("${CMAKE_SOURCE_DIR}/src/${name_internal}")
        
        #clone_git_repository(${repository} ${commit_id} ${CMAKE_SOURCE_DIR}/src ${name_internal})
        #download_repository(${repository} ${commit_id} ${CMAKE_SOURCE_DIR}/src ${name_internal} )

        # download_project(PROJ ${name_internal}
        #     GIT_REPOSITORY ${repository}
        #     GIT_TAG master
        #     UPDATE_DISCONNECTED 1
        #     PREFIX ""
        #     #DOWNLOAD_DIR  ${CMAKE_SOURCE_DIR}/src/${name_internal}
        #     SOURCE_DIR ${CMAKE_SOURCE_DIR}/src/${name_internal}
        # )
        execute_process(
            COMMAND git clone ${repository} ${name_internal}
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src
        )
        execute_process(
            COMMAND git checkout ${commit_id}
            WORKING_DIRECTORY ${DIR_MODULE_DIR}
        )

        if(EXISTS ${CMAKELISTS_PATH})
            if(MSVC)
                add_compile_options(/W0)
                add_compile_options(/Z7)
                set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Z7")
            endif()
            add_subdirectory(${DIR_MODULE_DIR})
            include(${CMAKECONFIG_PATH})
    
        else()
            message(STATUS "download failed.")
           
        endif()
    endif()
    file(MAKE_DIRECTORY ${DIR_PATH})
    file(MAKE_DIRECTORY ${DIR_MODULE_DIR})
    if(EXISTS ${CMAKELISTS_PATH})
            install(TARGETS ${name}
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib
            RUNTIME DESTINATION bin
        )
           install(TARGETS ${name}_test
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib
            RUNTIME DESTINATION bin
        )
     endif()



endfunction()