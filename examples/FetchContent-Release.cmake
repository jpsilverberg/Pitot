# Example CMakeLists.txt showing how to consume DSTL from a GitHub Release

cmake_minimum_required(VERSION 3.14)
project(MyApp)

include(FetchContent)

# Option 1: Download pre-built release archive (headers + CMake configs only)
# Upload dstl-X.Y.Z-headers.tar.gz to GitHub Releases first
FetchContent_Declare(dstl
  URL https://github.com/jpsilverberg/DSTL/releases/download/v1.0.0/dstl-1.0.0-headers.tar.gz
  # Optional: verify integrity
  # URL_HASH SHA256=<hash-of-your-tarball>
)

# This extracts the archive and makes find_package work
FetchContent_MakeAvailable(dstl)

# Now find the package (it's in the extracted location)
find_package(DSTL CONFIG REQUIRED COMPONENTS dmutex dpools)

add_executable(app main.cpp)
target_link_libraries(app PRIVATE dstl::dmutex dstl::dpools)
