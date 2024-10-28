execute_process(
        COMMAND git submodule update --init
        OUTPUT_VARIABLE LOCAL_HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        )
execute_process(
        COMMAND git submodule update 
        OUTPUT_VARIABLE LOCAL_HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        )