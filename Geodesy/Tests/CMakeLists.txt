cmake_minimum_required(VERSION 3.13)

# Find GTest
find_package(GTest REQUIRED)

# Test executable
add_executable(GeodesyTests
    test_geodesy.cpp
)

target_link_libraries(GeodesyTests
    PRIVATE
        dstl::Geodesy
        dstl::LinAlg
        GTest::GTest
        GTest::Main
)

# Add test
add_test(NAME GeodesyTests COMMAND GeodesyTests)
