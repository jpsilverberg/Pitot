cmake_minimum_required(VERSION 3.12)
project(pitot LANGUAGES CXX)

include(GNUInstallDirs)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(dstl_pitot INTERFACE)
add_library(dstl::pitot ALIAS dstl_pitot)

target_include_directories(dstl_pitot
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(dstl_pitot INTERFACE cxx_std_17)

set_target_properties(dstl_pitot PROPERTIES EXPORT_NAME pitot)

install(TARGETS dstl_pitot
  EXPORT DSTLTargets
)

install(DIRECTORY include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
)
