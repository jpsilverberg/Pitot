# Minimum CMake version required
cmake_minimum_required(VERSION 3.12)

# Project name and version
project(DSTL VERSION 1.0.0 LANGUAGES CXX)

# Project-wide build knobs so every component shares the same defaults.
set(DSTL_CXX_STANDARD 17 CACHE STRING "C++ standard used when building DSTL")
set_property(CACHE DSTL_CXX_STANDARD PROPERTY STRINGS 14 17 20 23)

option(DSTL_ENABLE_WARNINGS "Enable common warning flags for every DSTL target" ON)

set(CMAKE_CXX_STANDARD ${DSTL_CXX_STANDARD})
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(DSTL_ENABLE_WARNINGS)
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# Default to a sensible build type for single-configuration generators
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, like: Debug or Release" FORCE)
endif()

# Optional: Set output directories for binaries and libraries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Enable testing if needed (default OFF so downstream not slowed by tests)
option(BUILD_TESTING "Build tests" OFF)
if(BUILD_TESTING)
  enable_testing()
endif()

# Helper: only add a subdirectory if it exists and contains a CMakeLists.txt
macro(maybe_add_subdirectory subdir)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/CMakeLists.txt")
    message(STATUS "Adding subdirectory: ${subdir}")
    add_subdirectory(${subdir})
  else()
    message(STATUS "Skipping missing subdirectory: ${subdir}")
  endif()
endmacro()

set(DSTL_COMPONENT_DIRS
  Mutex
  Logger
  DataPools
  Numbers
  LinAlg
  Quadrature
  DGSTFEM
)

foreach(component ${DSTL_COMPONENT_DIRS})
  maybe_add_subdirectory(${component})
endforeach()

add_library(dstl INTERFACE)
target_include_directories(dstl
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Bits/include>
    $<INSTALL_INTERFACE:include>
)
target_compile_features(dstl INTERFACE cxx_std_${DSTL_CXX_STANDARD})

install(DIRECTORY include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY Bits/include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY DataPools
  DESTINATION .
  PATTERN "Tests" EXCLUDE
  PATTERN "build" EXCLUDE
  PATTERN ".git" EXCLUDE
  PATTERN "CMakeLists.txt" EXCLUDE
)

install(DIRECTORY Logger
  DESTINATION .
  PATTERN "Tests" EXCLUDE
  PATTERN "build" EXCLUDE
  PATTERN ".git" EXCLUDE
  PATTERN "CMakeLists.txt" EXCLUDE
)

install(DIRECTORY Mutex
  DESTINATION .
  PATTERN "Tests" EXCLUDE
  PATTERN "build" EXCLUDE
  PATTERN ".git" EXCLUDE
  PATTERN "CMakeLists.txt" EXCLUDE
)

set(DSTL_CONFIG_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/DSTL)

configure_package_config_file(
  cmake/DSTLConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/DSTLConfig.cmake"
  INSTALL_DESTINATION ${DSTL_CONFIG_INSTALL_DIR}
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/DSTLConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(TARGETS dstl
  EXPORT DSTLTargets
)

install(EXPORT DSTLTargets
  FILE DSTLTargets.cmake
  NAMESPACE dstl::
  DESTINATION ${DSTL_CONFIG_INSTALL_DIR}
)

export(EXPORT DSTLTargets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/DSTLTargets.cmake"
  NAMESPACE dstl::
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/DSTLConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/DSTLConfigVersion.cmake"
  DESTINATION ${DSTL_CONFIG_INSTALL_DIR}
)
