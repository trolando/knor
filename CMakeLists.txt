cmake_minimum_required(VERSION 3.0)
project(knor VERSION 1.0.0 LANGUAGES C CXX)

set(PROJECT_DESCRIPTION "Knor, a simple synthesis tool for HOA parity automata")
set(PROJECT_URL "https://github.com/trolando/knor")

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_FLAGS "-Wall -Wextra")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# Some info
message(STATUS "OS: ${CMAKE_SYSTEM_NAME}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID}")

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # require at least gcc 4.9
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9)
        message(FATAL_ERROR "GCC version must be at least 4.9!  " ${CMAKE_CXX_COMPILER_VERSION})
    endif()
endif()

include_directories(src)

find_package(FLEX)
find_package(BISON)
find_package(Boost REQUIRED filesystem iostreams regex system)

FLEX_TARGET(hoalexer src/hoa.l ${CMAKE_CURRENT_BINARY_DIR}/hoalexer.c DEFINES_FILE ${CMAKE_CURRENT_BINARY_DIR}/hoalexer.h)
BISON_TARGET(hoaparser src/hoa.y ${CMAKE_CURRENT_BINARY_DIR}/hoaparser.c DEFINES_FILE ${CMAKE_CURRENT_BINARY_DIR}/hoaparser.h)
ADD_FLEX_BISON_DEPENDENCY(hoalexer hoaparser)

find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    # Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(WARNING "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/sylvan/CMakeLists.txt" OR
   NOT EXISTS "${PROJECT_SOURCE_DIR}/oink/CMakeLists.txt")
    message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()

add_subdirectory("oink")
include_directories(oink/src)

add_subdirectory("sylvan")
include_directories(sylvan/src)

add_executable(knor src/knor.cpp src/symgame.cpp src/bisim.cpp src/simplehoa.c src/aiger.c ${BISON_hoaparser_OUTPUTS} ${FLEX_hoalexer_OUTPUTS})
target_link_libraries(knor PRIVATE oink sylvan)
target_link_libraries(knor PUBLIC "-static")
