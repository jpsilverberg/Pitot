# DSTLTesting.cmake
# Common testing setup for DSTL libraries

# Function to find and configure GTest
# Usage: find_and_configure_gtest()
function(find_and_configure_gtest)
  # Try CONFIG first (modern CMake), then MODULE mode
  find_package(GTest CONFIG QUIET)
  if(NOT GTest_FOUND)
    find_package(GTest QUIET)
  endif()
  
  # If still not found, fetch from GitHub
  if(NOT GTest_FOUND)
    message(STATUS "GTest not found, fetching via FetchContent...")
    include(FetchContent)
    FetchContent_Declare(
      googletest
      URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
    )
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
  else()
    message(STATUS "Using system-provided GTest: ${GTest_VERSION}")
  endif()
  
  # Include GoogleTest module for gtest_discover_tests
  include(GoogleTest)
endfunction()
