cmake_minimum_required(VERSION 3.13)

# Find GTest
find_package(GTest REQUIRED)

# Test executable
add_executable(GeodesyTests
    test_geodesy.cpp
    test_geodesy_reference.cpp
    test_geodesy_vincenty.cpp
    test_print_distances.cpp
    test_debug_failures.cpp
)

target_link_libraries(GeodesyTests
    PRIVATE
        dstl::Geodesy
        dstl::linalg
        GTest::gtest
        GTest::gtest_main
)

# Add test
add_test(NAME GeodesyTests COMMAND GeodesyTests)
