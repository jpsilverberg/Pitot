cmake_minimum_required(VERSION 3.13)
project(Numbers LANGUAGES CXX)

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

option(ENABLE_ASAN "Enable AddressSanitizer for tests" ON)

set(COMMON_WARNING_FLAGS -Wall -Wextra)
set(SANITIZE_FLAGS)
if(ENABLE_ASAN)
  list(APPEND SANITIZE_FLAGS -fsanitize=address)
endif()

get_filename_component(DSTL_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)

# When building DGSTFEM standalone ensure Numbers (and its dependencies) exist.
if(NOT TARGET dstl::numbers)
  if(EXISTS "${DSTL_ROOT}/DataMutex/CMakeLists.txt" AND EXISTS "${DSTL_ROOT}/DataLog/CMakeLists.txt")
    # Store the current BUILD_TESTING state
    set(NUMBERS_BUILD_TESTING ${BUILD_TESTING})
    
    # Temporarily disable testing for dependencies
    set(BUILD_TESTING OFF)
    
    if(NOT TARGET dstl::dmutex)
      add_subdirectory("${DSTL_ROOT}/DataMutex" "${CMAKE_CURRENT_BINARY_DIR}/deps/DataMutex")
    endif()
    add_subdirectory("${DSTL_ROOT}/DataLog" "${CMAKE_CURRENT_BINARY_DIR}/deps/DataLog")
    add_subdirectory("${DSTL_ROOT}/Numbers" "${CMAKE_CURRENT_BINARY_DIR}/deps/Numbers")
    
    # Restore testing setting for this project
    set(BUILD_TESTING ${NUMBERS_BUILD_TESTING})
  else()
    find_package(dstl CONFIG REQUIRED COMPONENTS numbers)
  endif()
endif()

add_library(dstl_dg0 STATIC
  Tests/DG0.cpp
)

add_library(dstl::dg0 ALIAS dstl_dg0)
target_include_directories(dstl_dg0
  PUBLIC
    $<BUILD_INTERFACE:${DSTL_ROOT}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(dstl_dg0 PUBLIC cxx_std_17)

target_compile_options(dstl_dg0 PUBLIC ${COMMON_WARNING_FLAGS})

target_link_libraries(dstl_dg0
  PUBLIC
    dstl::dlog
    dstl::numbers
)

set_target_properties(dstl_dg0 PROPERTIES EXPORT_NAME dg0)

if(SANITIZE_FLAGS)
  target_compile_options(dstl_dg0 PRIVATE ${SANITIZE_FLAGS})
  target_link_options(dstl_dg0 PRIVATE ${SANITIZE_FLAGS})
endif()

install(TARGETS dstl_dg0
  EXPORT DSTLTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY include/dstl/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dstl
)

install(FILES include/DG0.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

include(CTest)

if(BUILD_TESTING)
  enable_testing()

  find_package(GTest CONFIG QUIET)
  if(NOT GTest_FOUND)
    find_package(GTest QUIET)
  endif()

  if(NOT GTest_FOUND)
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

  include(GoogleTest)

  add_executable(
    DG0Tests
    Tests/DG0.cpp    
  )

  target_compile_options(DG0Tests PRIVATE ${COMMON_WARNING_FLAGS} ${SANITIZE_FLAGS})
  target_link_options(DG0Tests PRIVATE ${SANITIZE_FLAGS})
  target_link_libraries(
    DG0Tests
    PRIVATE
      dstl::dg0
      dstl::numbers
      GTest::gtest_main
      GTest::gtest
  )

  gtest_discover_tests(DG0Tests)

else()
  message(STATUS "BUILD_TESTING is OFF: skipping DataPools tests")
endif()
