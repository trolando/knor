cmake_minimum_required(VERSION 3.14)

project(knor
    VERSION 1.0.0
    DESCRIPTION "Knor, a synthesis tool for HOA parity automata using BDDs"
    HOMEPAGE_URL "https://github.com/trolando/knor"
    LANGUAGES C CXX
)

add_subdirectory(libs)

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
        src/gameconstructor.cpp
        src/knor.cpp
        src/mealy.cpp
        src/paritize.cpp
        src/permissive.cpp
        src/runstats.cpp
        src/simplehoa.cpp
        src/symgame.cpp
)
target_compile_features(knor PUBLIC c_std_11 cxx_std_17)
target_compile_options(knor PRIVATE -Wall -Wextra -fno-strict-aliasing -Wno-deprecated -Wno-unused-parameter)
target_link_libraries(knor PRIVATE oink::oink sylvan::sylvan libabc)
target_link_libraries(knor PRIVATE Boost::boost)
if(CMAKE_HOST_SYSTEM MATCHES Linux)
    target_link_libraries(knor PRIVATE "-static")
endif()
target_include_directories(knor PRIVATE src)

# Independent checker of the controllers against the automata, for the tests
add_executable(aigcheck)
target_sources(aigcheck
    PRIVATE
        test/aigcheck.cpp
        src/aiger.c
        src/paritize.cpp
        src/simplehoa.cpp
)
target_compile_features(aigcheck PUBLIC c_std_11 cxx_std_17)
target_compile_options(aigcheck PRIVATE -Wall -Wextra -fno-strict-aliasing -Wno-deprecated -Wno-unused-parameter)
target_include_directories(aigcheck PRIVATE src)
find_package(Threads REQUIRED)
target_link_libraries(aigcheck PRIVATE Threads::Threads)

# Tests (ctest / make test)
enable_testing()
add_test(NAME verdict_regression
    COMMAND bash ${CMAKE_SOURCE_DIR}/test/verdict_regression.sh $<TARGET_FILE:knor>)
add_test(NAME parser
    COMMAND bash ${CMAKE_SOURCE_DIR}/test/parser.sh $<TARGET_FILE:knor>)
add_test(NAME acceptance
    COMMAND bash ${CMAKE_SOURCE_DIR}/test/acceptance.sh $<TARGET_FILE:knor>)
add_test(NAME controllers
    COMMAND bash ${CMAKE_SOURCE_DIR}/test/controllers.sh $<TARGET_FILE:knor> $<TARGET_FILE:aigcheck>)
add_test(NAME aigcheck
    COMMAND bash ${CMAKE_SOURCE_DIR}/test/aigcheck.sh $<TARGET_FILE:knor> $<TARGET_FILE:aigcheck>)

# Optional: check the controllers and aigcheck itself with Spot, which shares no code with knor
find_program(AUTFILT autfilt)
find_package(Python3 COMPONENTS Interpreter)
if(AUTFILT AND Python3_Interpreter_FOUND)
    add_test(NAME spot_crosscheck
        COMMAND bash ${CMAKE_SOURCE_DIR}/test/spot_crosscheck.sh $<TARGET_FILE:knor> $<TARGET_FILE:aigcheck>
                ${AUTFILT} ${Python3_EXECUTABLE})
else()
    message(STATUS "autfilt (Spot) or Python 3 not found: the test spot_crosscheck is skipped (set AUTFILT to enable it)")
endif()
