cmake_minimum_required(VERSION 3.13)
project(Geodesy VERSION 1.0.0 LANGUAGES CXX)

# C++17 required
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Header-only library
add_library(Geodesy INTERFACE)
add_library(dstl::Geodesy ALIAS Geodesy)

target_include_directories(Geodesy INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Depend on LinAlg
target_link_libraries(Geodesy INTERFACE dstl::linalg)

# Installation
install(TARGETS Geodesy
    EXPORT DSTLTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)

install(DIRECTORY include/
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
)

# Export is handled by parent DSTL CMakeLists.txt

# Tests
option(BUILD_GEODESY_TESTS "Build Geodesy tests" ON)
if(BUILD_GEODESY_TESTS AND BUILD_TESTING)
    find_package(GTest QUIET)
    if(GTest_FOUND)
        enable_testing()
        add_subdirectory(Tests)
    else()
        message(STATUS "GTest not found, skipping Geodesy tests")
    endif()
endif()
