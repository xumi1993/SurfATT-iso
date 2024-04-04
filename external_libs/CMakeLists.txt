function(initialize_submodule DIRECTORY)
    find_package(Git QUIET)
    if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    # Update submodules as needed
        message(STATUS "Submodule update ${PROJECT_SOURCE_DIR}/external_libs/${DIRECTORY}")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endfunction(initialize_submodule)


initialize_submodule(yaml-cpp)

add_subdirectory(yaml-cpp EXCLUDE_FROM_ALL)