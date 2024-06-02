cmake_minimum_required(VERSION 3.14)

project(knor
    VERSION 1.0.0
    DESCRIPTION "Knor, a synthesis tool for HOA parity automata using BDDs"
    HOMEPAGE_URL "https://github.com/trolando/knor"
    LANGUAGES C CXX
)

add_subdirectory(libs)


# Parser generator
find_package(FLEX)
find_package(BISON)
FLEX_TARGET(hoalexer src/hoa.l ${CMAKE_CURRENT_BINARY_DIR}/hoalexer.c DEFINES_FILE ${CMAKE_CURRENT_BINARY_DIR}/hoalexer.h)
BISON_TARGET(hoaparser src/hoa.y ${CMAKE_CURRENT_BINARY_DIR}/hoaparser.c DEFINES_FILE ${CMAKE_CURRENT_BINARY_DIR}/hoaparser.h)
ADD_FLEX_BISON_DEPENDENCY(hoalexer hoaparser)

find_package(Boost REQUIRED random)

# Knor itself
add_executable(knor)
target_sources(knor
  PRIVATE
    src/aigcircuit.cpp
    src/aiger.c
    src/abcminimization.cpp
    src/aigencoder.cpp
    src/bddtools.cpp
    src/bisim.cpp
    src/explicitgame.cpp
    src/knor.cpp
    src/simplehoa.c
    src/symgame.cpp
    ${BISON_hoaparser_OUTPUTS}
    ${FLEX_hoalexer_OUTPUTS}
)
target_compile_features(knor PUBLIC c_std_11 cxx_std_17)
target_compile_options(knor PRIVATE -Wall -Wextra -fno-strict-aliasing -Wno-deprecated -Wno-unused-parameter)
target_link_libraries(knor PRIVATE oink::oink sylvan::sylvan libabc)
target_link_libraries(knor PRIVATE Boost::boost)
if(CMAKE_HOST_SYSTEM MATCHES Linux)
  target_link_libraries(knor PRIVATE "-static")
endif()
target_include_directories(knor PRIVATE src)
