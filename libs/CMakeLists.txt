find_package(lace QUIET)
if(NOT lace_FOUND)
    add_subdirectory(lace)
endif()

find_package(sylvan QUIET)
if(NOT sylvan_FOUND)
    add_subdirectory(sylvan)
endif()

find_package(oink QUIET)
if(NOT oink_FOUND)
    # knor reads no compressed games through Oink; without these optional dependencies (shared
    # libraries on most systems) knor still links statically on Linux
    set(CMAKE_DISABLE_FIND_PACKAGE_ZLIB ON)
    set(CMAKE_DISABLE_FIND_PACKAGE_BZip2 ON)
    set(CMAKE_DISABLE_FIND_PACKAGE_LibLZMA ON)
    add_subdirectory(oink)
    unset(CMAKE_DISABLE_FIND_PACKAGE_ZLIB)
    unset(CMAKE_DISABLE_FIND_PACKAGE_BZip2)
    unset(CMAKE_DISABLE_FIND_PACKAGE_LibLZMA)
endif()

find_package(abc QUIET)
if(NOT abc_FOUND)
    set (READLINE_FOUND FALSE)
    add_definitions (-DABC_USE_STDINT_H)
    # not ABC's own test suite, which downloads googletest when it is configured
    set(ABC_SKIP_TESTS ON CACHE BOOL "Skip ABC's test suite" FORCE)
    add_subdirectory(abc)
    # ABC compiles its C++ without exceptions and passes its flags on to the targets that link it;
    # knor uses exceptions, so keep that one flag to ABC itself
    get_target_property(ABC_INTERFACE_OPTIONS libabc INTERFACE_COMPILE_OPTIONS)
    if(ABC_INTERFACE_OPTIONS)
        list(FILTER ABC_INTERFACE_OPTIONS EXCLUDE REGEX "-fno-exceptions")
        set_target_properties(libabc PROPERTIES INTERFACE_COMPILE_OPTIONS "${ABC_INTERFACE_OPTIONS}")
    endif()
endif()
