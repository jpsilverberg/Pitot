# Example using CPM.cmake for more flexible package management
# CPM can handle both source and binary packages

cmake_minimum_required(VERSION 3.14)
project(MyApp)

# Download CPM.cmake if not already present
set(CPM_DOWNLOAD_VERSION 0.38.1)
if(NOT EXISTS "${CMAKE_BINARY_DIR}/cmake/CPM.cmake")
  file(DOWNLOAD
    "https://github.com/cpm-cmake/CPM.cmake/releases/download/v${CPM_DOWNLOAD_VERSION}/CPM.cmake"
    "${CMAKE_BINARY_DIR}/cmake/CPM.cmake"
  )
endif()
include(${CMAKE_BINARY_DIR}/cmake/CPM.cmake)

# Option A: Use a GitHub release archive (headers only, no tests)
CPMAddPackage(
  NAME dstl
  URL https://github.com/jpsilverberg/DSTL/releases/download/v1.0.0/dstl-1.0.0-headers.tar.gz
  # URL_HASH SHA256=<hash>
)

# Option B: Use git but exclude tests (fallback if no release exists)
# CPMAddPackage(
#   NAME dstl
#   GITHUB_REPOSITORY jpsilverberg/DSTL
#   GIT_TAG v1.0.0
#   OPTIONS
#     "BUILD_TESTING OFF"
# )

find_package(DSTL CONFIG REQUIRED COMPONENTS dmutex dpools)

add_executable(app main.cpp)
target_link_libraries(app PRIVATE dstl::dmutex dstl::dpools)
