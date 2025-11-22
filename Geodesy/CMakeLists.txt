cmake_minimum_required(VERSION 3.13)
project(Geodesy VERSION 1.0.0 LANGUAGES CXX)

# Determine if this is a standalone build or part of DSTL
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(GEODESY_STANDALONE ON)
else()
    set(GEODESY_STANDALONE OFF)
endif()

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

# Handle LinAlg dependency
if(GEODESY_STANDALONE)
    # Standalone build: LinAlg must be found or provided
    if(NOT TARGET dstl::linalg)
        # Try to find LinAlg in parent directory (if building alongside DSTL)
        if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../LinAlg/CMakeLists.txt")
            add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../LinAlg ${CMAKE_CURRENT_BINARY_DIR}/LinAlg)
        else()
            message(STATUS "LinAlg not found. Geodesy requires dstl::linalg.")
            message(STATUS "Please ensure LinAlg is available or build as part of DSTL.")
        endif()
    endif()
endif()

# Link to LinAlg if available
if(TARGET dstl::linalg)
    target_link_libraries(Geodesy INTERFACE dstl::linalg)
else()
    message(WARNING "dstl::linalg target not found. Geodesy may not build correctly.")
endif()

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
