# Minimum CMake version required
cmake_minimum_required(VERSION 3.12)

# Project name and version
project(DSTL VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default to a sensible build type for single-configuration generators
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, like: Debug or Release" FORCE)
endif()

# Optional: Set output directories for binaries and libraries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Add include directories for DSTL
include_directories(${CMAKE_SOURCE_DIR}/include)


# Enable testing if needed
option(BUILD_TESTING "Build tests" ON)
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

maybe_add_subdirectory(DataMutex)
maybe_add_subdirectory(DataLog)
maybe_add_subdirectory(DataPools)

# Optional: Add tests and examples directories
# if(BUILD_TESTING)
#   add_subdirectory(tests)
# endif()
# add_subdirectory(examples)

# Installation rules (for headers, libraries, and CMake config)
# include(GNUInstallDirs)
# install(EXPORT MyProjectTargets
#   FILE MyProjectTargets.cmake
#   NAMESPACE MyProject::
#   DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MyProject
# )

# # Generate and install a CMake package configuration file
# include(CMakePackageConfigHelpers)
# write_basic_package_version_file(
#   "${CMAKE_BINARY_DIR}/MyProjectConfigVersion.cmake"
#   VERSION ${PROJECT_VERSION}
#   COMPATIBILITY SameMajorVersion
# )

# install(FILES
#   "${CMAKE_CURRENT_SOURCE_DIR}/cmake/MyProjectConfig.cmake"
#   "${CMAKE_BINARY_DIR}/MyProjectConfigVersion.cmake"
#   DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MyProject
# )
