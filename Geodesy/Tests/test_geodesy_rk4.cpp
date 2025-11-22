// test_geodesy_rk4.cpp - RK4 geodesic integration tests
// Tests 4th-order Runge-Kutta integration for geodesic paths with altitude profiles

#include <dstl/Geodesy.h>
#include <gtest/gtest.h>
#include "test_reference_data.h"
#include <cmath>
#include <functional>

using namespace dstl::geo;
using namespace dstl::geo::test_data;
using namespace dstl::linalg;

// ==================== Constant Altitude Tests ====================

TEST(RK4Test, ConstantAltitude_ShortDistance)
{
    auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 10000.0);
    
    // Move 100km North at constant altitude
    double distance = 100000.0;
    double bearing = 0.0;  // North
    int steps = 32;
    
    auto end_rk4 = start.move_geodesic_rk4(bearing, distance, nullptr, steps);
    auto end_vincenty = start.move_by_bearing_accurate(distance, bearing);
    
    // RK4 should match Vincenty for constant altitude
    EXPECT_NEAR(end_rk4.latitude_deg(), end_vincenty.latitude_deg(), 0.001)
        << "RK4 latitude should match Vincenty";
    EXPECT_NEAR(end_rk4.longitude_deg(), end_vincenty.longitude_deg(), 0.001)
        << "RK4 longitude should match Vincenty";
    EXPECT_NEAR(end_rk4.altitude(), start.altitude(), 1.0)
        << "Altitude should be preserved";
}

TEST(RK4Test, ConstantAltitude_MediumDistance)
{
    auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 5000.0);
    
    // Move 500km at 45° (NE)
    double distance = 500000.0;
    double bearing = PI / 4.0;
    
    auto end_rk4 = start.move_geodesic_rk4(bearing, distance, nullptr, 32);
    auto end_vincenty = start.move_by_bearing_accurate(distance, bearing);
    
    // Should be very close (within 10m over 500km)
    double position_error = end_rk4.distance_to(end_vincenty);
    EXPECT_LT(position_error, 10.0)
        << "RK4 should match Vincenty within 10m over 500km";
}

TEST(RK4Test, ConstantAltitude_LongDistance)
{
    auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
    
    // Move 5000km East
    double distance = 5000000.0;
    double bearing = PI / 2.0;  // East
    
    auto end_rk4 = start.move_geodesic_rk4(bearing, distance, nullptr, 64);
    auto end_vincenty = start.move_by_bearing_accurate(distance, bearing);
    
    // Should be close even over long distance (within 100m)
    double position_error = end_rk4.distance_to(end_vincenty);
    EXPECT_LT(position_error, 100.0)
        << "RK4 should match Vincenty within 100m over 5000km";
}

TEST(RK4Test, ConstantAltitude_AllDirections)
{
    auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 1000.0);
    double distance = 200000.0;  // 200km
    
    // Test all cardinal directions
    std::vector<std::pair<std::string, double>> directions = {
        {"North", 0.0},
        {"NE", PI / 4.0},
        {"East", PI / 2.0},
        {"SE", 3.0 * PI / 4.0},
        {"South", PI},
        {"SW", 5.0 * PI / 4.0},
        {"West", 3.0 * PI / 2.0},
        {"NW", 7.0 * PI / 4.0}
    };
    
    for (const auto& [name, bearing] : directions)
    {
        auto end_rk4 = origin.move_geodesic_rk4(bearing, distance, nullptr, 32);
        auto end_vincenty = origin.move_by_bearing_accurate(distance, bearing);
        
        double error = end_rk4.distance_to(end_vincenty);
        EXPECT_LT(error, 10.0)
            << "RK4 error for " << name << " direction: " << error << "m";
    }
}

// ==================== Altitude Profile Tests ====================

TEST(RK4Test, AltitudeProfile_LinearClimb)
{
    auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
    
    // Linear climb: 0m to 10,000m over 500km
    double distance = 500000.0;
    double bearing = 0.0;  // North
    auto profile = [](double d) { return d * 0.02; };  // 10,000m / 500,000m = 0.02
    
    auto end = start.move_geodesic_rk4(bearing, distance, profile, 32);
    
    // Should end at 10,000m altitude
    EXPECT_NEAR(end.altitude(), 10000.0, 10.0)
        << "Final altitude should be 10,000m";
    
    // Should have moved North
    EXPECT_GT(end.latitude_deg(), start.latitude_deg())
        << "Should have moved North";
}

TEST(RK4Test, AltitudeProfile_LinearDescent)
{
    auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 10000.0);
    
    // Linear descent: 10,000m to 0m over 500km
    double distance = 500000.0;
    double bearing = PI / 2.0;  // East
    auto profile = [](double d) { return 10000.0 - d * 0.02; };
    
    auto end = start.move_geodesic_rk4(bearing, distance, profile, 32);
    
    // Should end at 0m altitude
    EXPECT_NEAR(end.altitude(), 0.0, 10.0)
        << "Final altitude should be 0m";
    
    // Should have moved East
    EXPECT_GT(end.longitude_deg(), start.longitude_deg())
        << "Should have moved East";
}

TEST(RK4Test, AltitudeProfile_Parabolic)
{
    auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
    
    // Parabolic: climb to 10,000m at midpoint, return to 0m
    double distance = 500000.0;
    double bearing = PI / 4.0;  // NE
    auto profile = [](double d) {
        double x = d / 500000.0;  // Normalize to [0,1]
        return 10000.0 * 4.0 * x * (1.0 - x);  // Peak at x=0.5
    };
    
    auto end = start.move_geodesic_rk4(bearing, distance, profile, 64);
    
    // Should end near starting altitude
    EXPECT_NEAR(end.altitude(), 0.0, 50.0)
        << "Final altitude should return to ~0m";
    
    // Should have moved diagonally
    EXPECT_GT(end.latitude_deg(), start.latitude_deg());
    EXPECT_GT(end.longitude_deg(), start.longitude_deg());
}

TEST(RK4Test, AltitudeProfile_Constant)
{
    auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 5000.0);
    
    // Constant altitude profile
    double distance = 300000.0;
    double bearing = PI;  // South
    auto profile = [](double) { return 5000.0; };
    
    auto end_with_profile = start.move_geodesic_rk4(bearing, distance, profile, 32);
    auto end_without_profile = start.move_geodesic_rk4(bearing, distance, nullptr, 32);
    
    // Both should give same result
    EXPECT_NEAR(end_with_profile.latitude_deg(), end_without_profile.latitude_deg(), 0.001);
    EXPECT_NEAR(end_with_profile.longitude_deg(), end_without_profile.longitude_deg(), 0.001);
    EXPECT_NEAR(end_with_profile.altitude(), 5000.0, 1.0);
}

// ==================== Step Count Sensitivity ====================

TEST(RK4Test, StepCount_Convergence)
{
    auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
    double distance = 1000000.0;  // 1000km
    double bearing = 0.0;
    
    // Test with different step counts
    auto end_8 = start.move_geodesic_rk4(bearing, distance, nullptr, 8);
    auto end_16 = start.move_geodesic_rk4(bearing, distance, nullptr, 16);
    auto end_32 = start.move_geodesic_rk4(bearing, distance, nullptr, 32);
    auto end_64 = start.move_geodesic_rk4(bearing, distance, nullptr, 64);
    
    // More steps should converge to same answer
    double error_16_32 = end_16.distance_to(end_32);
    double error_32_64 = end_32.distance_to(end_64);
    
    // Error should decrease with more steps
    EXPECT_LT(error_32_64, error_16_32)
        << "Error should decrease with more steps";
    
    // With 64 steps, should be very accurate
    EXPECT_LT(error_32_64, 10.0)
        << "32 vs 64 steps should differ by < 10m";
}

TEST(RK4Test, StepCount_MinimumSteps)
{
    auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
    double distance = 100000.0;
    double bearing = 0.0;
    
    // Test with very few steps (should still work)
    auto end_4 = start.move_geodesic_rk4(bearing, distance, nullptr, 4);
    auto end_32 = start.move_geodesic_rk4(bearing, distance, nullptr, 32);
    
    // Even with 4 steps, should be reasonably close
    double error = end_4.distance_to(end_32);
    EXPECT_LT(error, 100.0)
        << "Even 4 steps should give < 100m error over 100km";
}

TEST(RK4Test, StepCount_VeryShortDistance)
{
    auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
    double distance = 100.0;  // 100m - very short
    double bearing = 0.0;
    
    // For very short distances, should fall back to direct method
    auto end = start.move_geodesic_rk4(bearing, distance, nullptr, 32);
    
    // Should still give reasonable result
    double actual_dist = start.geodesic_distance_to(end);
    EXPECT_NEAR(actual_dist, distance, 1.0);
}

// ==================== Accuracy vs Vincenty ====================

TEST(RK4Test, Accuracy_VsVincenty_ShortRange)
{
    auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
    double distance = 50000.0;  // 50km
    double bearing = PI / 3.0;  // 60°
    
    auto end_rk4 = start.move_geodesic_rk4(bearing, distance, nullptr, 32);
    auto end_vincenty = start.move_by_bearing_accurate(distance, bearing);
    
    double error = end_rk4.distance_to(end_vincenty);
    EXPECT_LT(error, 1.0)
        << "RK4 should match Vincenty within 1m over 50km";
}

TEST(RK4Test, Accuracy_VsVincenty_LongRange)
{
    auto start = ECEFCoordinate::from_geodetic_deg(0.0, 0.0, 0.0);
    double distance = 10000000.0;  // 10,000km
    double bearing = PI / 4.0;
    
    auto end_rk4 = start.move_geodesic_rk4(bearing, distance, nullptr, 128);
    auto end_vincenty = start.move_by_bearing_accurate(distance, bearing);
    
    double error = end_rk4.distance_to(end_vincenty);
    EXPECT_LT(error, 1000.0)
        << "RK4 should match Vincenty within 1km over 10,000km";
}

TEST(RK4Test, Accuracy_NearPole)
{
    auto start = ECEFCoordinate::from_geodetic_deg(85.0, 0.0, 0.0);
    double distance = 500000.0;  // 500km
    double bearing = PI / 2.0;  // East
    
    auto end_rk4 = start.move_geodesic_rk4(bearing, distance, nullptr, 64);
    auto end_vincenty = start.move_by_bearing_accurate(distance, bearing);
    
    double error = end_rk4.distance_to(end_vincenty);
    EXPECT_LT(error, 100.0)
        << "RK4 should work near pole within 100m";
}

// ==================== Convenience Functions ====================

TEST(RK4Test, ConvenienceFunction_ConstantAltitude)
{
    auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 5000.0);
    double distance = 200000.0;
    double bearing = 0.0;
    
    auto end1 = start.move_geodesic_rk4_constant_altitude(bearing, distance, 32);
    auto end2 = start.move_geodesic_rk4(bearing, distance, nullptr, 32);
    
    // Both should give same result
    EXPECT_NEAR(end1.latitude_deg(), end2.latitude_deg(), 0.0001);
    EXPECT_NEAR(end1.longitude_deg(), end2.longitude_deg(), 0.0001);
    EXPECT_NEAR(end1.altitude(), end2.altitude(), 0.1);
}

TEST(RK4Test, FlyConstantAltitude)
{
    auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 10000.0);
    
    // Fly 270nm (500km) at 90° (East) at 35,000 ft
    double track_deg = 90.0;
    double distance_nm = 270.0;
    double altitude_ft = 35000.0;
    
    auto end = start.fly_constant_altitude(track_deg, distance_nm, altitude_ft, 64);
    
    // Should end at cruise altitude
    EXPECT_NEAR(end.altitude(), altitude_ft * FT_TO_M, 50.0)
        << "Should reach cruise altitude";
    
    // Should have moved East
    EXPECT_GT(end.longitude_deg(), start.longitude_deg());
}

TEST(RK4Test, FlyAltitudeProfile)
{
    auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
    
    // Fly with climb profile: climb to 30,000 ft over 400nm
    double track_deg = 45.0;  // NE
    double distance_nm = 400.0;
    auto profile = [](double d_nm) { return d_nm * 75.0; };  // 30,000 ft / 400 nm = 75 ft/nm
    
    auto end = start.fly_altitude_profile(track_deg, distance_nm, profile, 32);
    
    // Should have climbed to ~30,000 ft
    EXPECT_GT(end.altitude(), 8000.0)  // > 8000m (~26,000 ft)
        << "Should have climbed";
    
    // Should have moved diagonally
    EXPECT_GT(end.latitude_deg(), start.latitude_deg());
    EXPECT_GT(end.longitude_deg(), start.longitude_deg());
}

// ==================== Edge Cases ====================

TEST(RK4Test, EdgeCase_ZeroDistance)
{
    auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 1000.0);
    
    auto end = start.move_geodesic_rk4(0.0, 0.0, nullptr, 32);
    
    // Should stay at same position
    EXPECT_NEAR(end.latitude_deg(), start.latitude_deg(), 0.0001);
    EXPECT_NEAR(end.longitude_deg(), start.longitude_deg(), 0.0001);
    EXPECT_NEAR(end.altitude(), start.altitude(), 0.1);
}

TEST(RK4Test, EdgeCase_VeryLongDistance)
{
    auto start = ECEFCoordinate::from_geodetic_deg(0.0, 0.0, 0.0);
    
    // Nearly quarter of Earth's circumference
    double distance = 9000000.0;  // 9,000km
    double bearing = PI / 2.0;  // East
    
    auto end = start.move_geodesic_rk4(bearing, distance, nullptr, 128);
    
    // Should have moved significantly
    EXPECT_GT(end.longitude_deg(), 70.0)
        << "Should have moved far East";
}

TEST(RK4Test, EdgeCase_AltitudeProfileNegative)
{
    auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 5000.0);
    
    // Descend below sea level
    double distance = 100000.0;
    auto profile = [](double d) { return 5000.0 - d * 0.06; };  // End at -1000m
    
    auto end = start.move_geodesic_rk4(0.0, distance, profile, 32);
    
    // Should handle negative altitude
    EXPECT_LT(end.altitude(), 0.0)
        << "Should handle negative altitude";
}

// ==================== Consistency Tests ====================

TEST(RK4Test, Consistency_RoundTrip)
{
    auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 5000.0);
    double distance = 500000.0;
    double bearing = PI / 3.0;
    
    // Go forward
    auto mid = start.move_geodesic_rk4(bearing, distance, nullptr, 32);
    
    // Go back
    double bearing_back = mid.bearing_to(start);
    double dist_back = mid.geodesic_distance_to(start);
    auto end = mid.move_geodesic_rk4(bearing_back, dist_back, nullptr, 32);
    
    // Should return close to start
    double error = start.distance_to(end);
    EXPECT_LT(error, 100.0)
        << "Round-trip should return within 100m";
}

TEST(RK4Test, Consistency_MultipleHops)
{
    auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
    double bearing = 0.0;
    
    // One 1000km hop
    auto end_single = start.move_geodesic_rk4(bearing, 1000000.0, nullptr, 64);
    
    // Two 500km hops
    auto mid = start.move_geodesic_rk4(bearing, 500000.0, nullptr, 32);
    auto end_double = mid.move_geodesic_rk4(bearing, 500000.0, nullptr, 32);
    
    // Should be close
    double error = end_single.distance_to(end_double);
    EXPECT_LT(error, 50.0)
        << "Multiple hops should give consistent result";
}

TEST(RK4Test, Consistency_WithAltitudeProfile)
{
    auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
    double distance = 500000.0;
    double bearing = 0.0;
    
    // Profile that ends at same altitude
    auto profile = [](double d) {
        double x = d / 500000.0;
        return 5000.0 * std::sin(2.0 * PI * x);  // Sine wave, ends at 0
    };
    
    auto end = start.move_geodesic_rk4(bearing, distance, profile, 64);
    
    // Should end near starting altitude
    EXPECT_NEAR(end.altitude(), 0.0, 100.0)
        << "Sine wave profile should return to starting altitude";
}

// ==================== Performance Characteristics ====================

TEST(RK4Test, Performance_ReasonableSteps)
{
    auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
    
    // Should complete quickly even with many steps
    auto end = start.move_geodesic_rk4(0.0, 5000000.0, nullptr, 256);
    
    // Just verify it completes
    EXPECT_NE(end.latitude_deg(), start.latitude_deg());
}

TEST(RK4Test, Performance_ComplexProfile)
{
    auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
    
    // Complex altitude profile
    auto profile = [](double d) {
        return 5000.0 * std::sin(d / 50000.0) + 3000.0 * std::cos(d / 30000.0) + 10000.0;
    };
    
    // Should complete quickly
    auto end = start.move_geodesic_rk4(0.0, 1000000.0, profile, 128);
    
    // Just verify it completes
    EXPECT_GT(end.altitude(), 0.0);
}
