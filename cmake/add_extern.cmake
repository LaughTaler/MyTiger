include("${CMAKE_SOURCE_DIR}/cmake/DownloadProject.cmake")
set(DIR_PATH_EXTERN "${CMAKE_SOURCE_DIR}/extern")

function(add_extern_git name_internal repository commit_id)
    message("try to add extern ${name_internal}")
    set(DIR_MODULE_DIR "${DIR_PATH_EXTERN}/${name_internal}")
    set(CMAKELISTS_PATH "${DIR_MODULE_DIR}/CMakeLists.txt")

    file(MAKE_DIRECTORY ${DIR_PATH_EXTERN})
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

        if(MSVC)
                add_compile_options(/W0)
                add_compile_options(/Z7)
                set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Z7")
            endif()
        #add_subdirectory(${DIR_MODULE_DIR})
        message(STATUS "Directory ${CMAKELISTS_PATH} added.")
    else()
        message(STATUS "try to download ${name} from ${repository}")
        
        execute_process(
            COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_SOURCE_DIR}/extern/${name_internal}
            RESULT_VARIABLE result
        )
        execute_process(
            COMMAND git clone ${repository} ${name_internal}
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/extern
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
        else()
            message(STATUS "download failed.")
           
        endif()
    endif()




endfunction()