cmake_minimum_required(VERSION 3.15)
project(dstl_package_smoke LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(DSTL CONFIG REQUIRED COMPONENTS dmutex dlog)

add_executable(dstl_package_smoke main.cpp)
target_link_libraries(dstl_package_smoke
  PRIVATE
    dstl
    dstl::dmutex
    dstl::dlog
)
