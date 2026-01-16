cmake_minimum_required(VERSION 3.22.1)
project(Pitot VERSION 1.0.0 LANGUAGES CXX)

include(GNUInstallDirs)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# C++17 required (Pitot uses inline constexpr)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

get_filename_component(DSTL_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)

# Use common DSTL modules
list(APPEND CMAKE_MODULE_PATH "${DSTL_ROOT}/cmake")
include(DSTLTesting)

# Enable testing if needed
option(BUILD_TESTING "Build tests" ON)
if(BUILD_TESTING)
  enable_testing()
endif()

# Header-only library
add_library(dstl_pitot INTERFACE)

add_library(dstl::pitot ALIAS dstl_pitot)
add_library(dstl::Pitot ALIAS dstl_pitot)

target_include_directories(dstl_pitot
  INTERFACE
    $<BUILD_INTERFACE:${DSTL_ROOT}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(dstl_pitot INTERFACE cxx_std_17)

set_target_properties(dstl_pitot PROPERTIES EXPORT_NAME pitot)

#Installation
install(TARGETS dstl_pitot
  EXPORT DSTLTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY include/dstl/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dstl
  FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
)

# Tests
if(BUILD_TESTING)
  enable_testing()
  find_and_configure_gtest()
  
  # Add test directory if it exists
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/Tests")
    add_subdirectory(Tests)
  endif()
endif()
