cmake_minimum_required(VERSION 3.13)
project(LinAlg LANGUAGES CXX)

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

# Enable testing if needed
option(BUILD_TESTING "Build tests" ON)
if(BUILD_TESTING)
  enable_testing()
endif()

# When building Numbers standalone ensure its dstl dependencies are available.
if(NOT TARGET dstl::dlog)
  if(EXISTS "${DSTL_ROOT}/Mutex/CMakeLists.txt" AND EXISTS "${DSTL_ROOT}/Logger/CMakeLists.txt")
    message(STATUS "Using local Mutex and Logger dependencies")
    # Store the current BUILD_TESTING state
    set(NUMBERS_BUILD_TESTING ${BUILD_TESTING})
    
    # Temporarily disable testing for dependencies
    set(BUILD_TESTING OFF)
    
    if(NOT TARGET dstl::dmutex)
      add_subdirectory("${DSTL_ROOT}/Mutex" "${CMAKE_CURRENT_BINARY_DIR}/deps/Mutex")
    endif()
    add_subdirectory("${DSTL_ROOT}/Logger" "${CMAKE_CURRENT_BINARY_DIR}/deps/Logger")
    
    # Restore testing setting for this project
    set(BUILD_TESTING ${NUMBERS_BUILD_TESTING})
  else()
    find_package(dstl CONFIG REQUIRED COMPONENTS dlog)
  endif()
endif()

add_library(dstl_linalg STATIC
  src/LinAlg.cpp
)

add_library(dstl::linalg ALIAS dstl_linalg)

target_include_directories(dstl_linalg
  PUBLIC
    $<BUILD_INTERFACE:${DSTL_ROOT}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(dstl_linalg PUBLIC cxx_std_17)

target_compile_options(dstl_linalg PUBLIC ${COMMON_WARNING_FLAGS})

target_link_libraries(dstl_linalg
  PUBLIC
    dstl::dlog
)

set_target_properties(dstl_linalg PROPERTIES EXPORT_NAME linalg)

if(SANITIZE_FLAGS)
  target_compile_options(dstl_linalg PRIVATE ${SANITIZE_FLAGS})
  target_link_options(dstl_linalg PRIVATE ${SANITIZE_FLAGS})
endif()

install(TARGETS dstl_linalg
  EXPORT DSTLTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY include/dstl/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dstl
)

install(FILES include/LinAlg.h
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
    LinAlgTests
    Tests/LinAlgTests.cpp
    Tests/VecTests.cpp
    Tests/MatTests.cpp
  )

  target_compile_options(LinAlgTests PRIVATE ${COMMON_WARNING_FLAGS} ${SANITIZE_FLAGS})
  target_link_options(LinAlgTests PRIVATE ${SANITIZE_FLAGS})

  target_link_libraries(
    LinAlgTests
    PRIVATE
      dstl::linalg
      GTest::gtest_main
      GTest::gtest
  )

  gtest_discover_tests(LinAlgTests)

else()
  message(STATUS "BUILD_TESTING is OFF: skipping LinAlg tests")
endif()
