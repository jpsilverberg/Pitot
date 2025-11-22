// test_geodesy_100percent.cpp - Final tests to reach 100% coverage
// Tests for any remaining uncovered functions and edge cases

#include <dstl/Geodesy.h>
#include <gtest/gtest.h>
#include "test_reference_data.h"
#include <cmath>

using namespace dstl::geo;
using namespace dstl::geo::test_data;
using namespace dstl::linalg;

// ==================== Remaining RK4 Functions ====================

TEST(Coverage100Test, DisplaceConstantAltitudeRK4)
{
    auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 10000.0);
    
    // Displace using RK4 with constant altitude
    double east_m = 100000.0;   // 100km East
    double north_m = 100000.0;  // 100km North
    
    auto target = origin.displace_constant_altitude_rk4(east_m, north_m, 64);
    
    // Should have moved diagonally
    EXPECT_GT(target.latitude_deg(), origin.latitude_deg());
    EXPECT_GT(target.longitude_deg(), origin.longitude_deg());
    
    // Altitude should be approximately preserved
    EXPECT_NEAR(target.altitude(), origin.altitude(), 50.0);
}

TEST(Coverage100Test, DisplaceConstantAltitudeRK4_ZeroDisplacement)
{
    auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 5000.0);
    
    auto target = origin.displace_constant_altitude_rk4(0.0, 0.0, 32);
    
    // Should stay at same position
    EXPECT_NEAR(target.latitude_deg(), origin.latitude_deg(), 0.0001);
    EXPECT_NEAR(target.longitude_deg(), origin.longitude_deg(), 0.0001);
    EXPECT_NEAR(target.altitude(), origin.altitude(), 0.1);
}

TEST(Coverage100Test, DisplaceConstantAltitudeRK4_VsSpherical)
{
    auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
    
    double east_m = 50000.0;
    double north_m = 50000.0;
    Vec3 enu{east_m, north_m, 0.0};
    
    auto rk4_result = origin.displace_constant_altitude_rk4(east_m, north_m, 32);
    auto spherical_result = origin.apply_enu_displacement_spherical(enu, 0.0, 3);
    
    // Should be reasonably close
    double diff = rk4_result.distance_to(spherical_result);
    EXPECT_LT(diff, 200.0)
        << "RK4 and spherical should be similar for medium distances";
}

// ==================== NED Constant Altitude ====================

TEST(Coverage100Test, NEDDisplacement_ConstantAltitude)
{
    auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 5000.0);
    
    // Move 10km N, 5km E, climb 500m
    Vec3 ned{10000.0, 5000.0, 0.0};
    double alt_change = 500.0;
    
    auto target = origin.apply_ned_displacement_constant_altitude(ned, alt_change);
    
    // Should have moved
    EXPECT_GT(target.latitude_deg(), origin.latitude_deg());
    EXPECT_GT(target.longitude_deg(), origin.longitude_deg());
    
    // Altitude change should be approximately correct
    double actual_alt_change = target.altitude() - origin.altitude();
    EXPECT_NEAR(actual_alt_change, alt_change, 20.0);
}

TEST(Coverage100Test, NEDDisplacement_ConstantAltitude_ZeroChange)
{
    auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 3000.0);
    
    Vec3 ned{5000.0, 5000.0, 0.0};
    auto target = origin.apply_ned_displacement_constant_altitude(ned, 0.0);
    
    // Altitude should be approximately preserved
    EXPECT_NEAR(target.altitude(), origin.altitude(), 10.0);
}

TEST(Coverage100Test, NEDDisplacement_ConstantAltitude_Descent)
{
    auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 10000.0);
    
    // Descend 3000m
    Vec3 ned{20000.0, 10000.0, 0.0};
    double alt_change = -3000.0;
    
    auto target = origin.apply_ned_displacement_constant_altitude(ned, alt_change);
    
    // Should have descended
    EXPECT_LT(target.altitude(), origin.altitude());
    
    double actual_alt_change = target.altitude() - origin.altitude();
    EXPECT_NEAR(actual_alt_change, alt_change, 30.0);
}

// ==================== Operator Overloads ====================

TEST(Coverage100Test, OperatorMinus_VectorOffset)
{
    auto coord = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 1000.0);
    Vec3 offset{1000.0, 2000.0, 3000.0};
    
    auto result = coord - offset;
    
    // Should have subtracted offset
    EXPECT_NEAR(result.x(), coord.x() - 1000.0, 0.001);
    EXPECT_NEAR(result.y(), coord.y() - 2000.0, 0.001);
    EXPECT_NEAR(result.z(), coord.z() - 3000.0, 0.001);
}

TEST(Coverage100Test, OperatorMinusEquals_VectorOffset)
{
    auto coord = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 1000.0);
    double orig_x = coord.x();
    double orig_y = coord.y();
    double orig_z = coord.z();
    
    Vec3 offset{1000.0, 2000.0, 3000.0};
    coord -= offset;
    
    // Should have subtracted offset in place
    EXPECT_NEAR(coord.x(), orig_x - 1000.0, 0.001);
    EXPECT_NEAR(coord.y(), orig_y - 2000.0, 0.001);
    EXPECT_NEAR(coord.z(), orig_z - 3000.0, 0.001);
}

TEST(Coverage100Test, OperatorPlus_Commutative)
{
    auto coord = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 1000.0);
    Vec3 offset{1000.0, 2000.0, 3000.0};
    
    auto result1 = coord + offset;
    auto result2 = offset + coord;
    
    // Should be commutative
    EXPECT_NEAR(result1.x(), result2.x(), 0.001);
    EXPECT_NEAR(result1.y(), result2.y(), 0.001);
    EXPECT_NEAR(result1.z(), result2.z(), 0.001);
}

// ==================== Default Constructor ====================

TEST(Coverage100Test, DefaultConstructor)
{
    ECEFCoordinate coord;
    
    // Default should be at origin
    EXPECT_DOUBLE_EQ(coord.x(), 0.0);
    EXPECT_DOUBLE_EQ(coord.y(), 0.0);
    EXPECT_DOUBLE_EQ(coord.z(), 0.0);
}

TEST(Coverage100Test, DefaultConstructor_GeodeticAccess)
{
    ECEFCoordinate coord;  // At Earth center
    
    // Accessing geodetic coordinates of Earth center
    // This is a degenerate case but should not crash
    double lat = coord.latitude();
    double lon = coord.longitude();
    double alt = coord.altitude();
    
    // Should return some values (exact values depend on implementation)
    // Just verify it doesn't crash and returns finite values
    EXPECT_TRUE(std::isfinite(lat));
    EXPECT_TRUE(std::isfinite(lon));
    EXPECT_TRUE(std::isfinite(alt));
}

// ==================== Edge Cases for Frames ====================

TEST(Coverage100Test, ENUFrame_AtPole)
{
    // Test ENU frame at North Pole
    auto coord = ECEFCoordinate::from_geodetic_deg(90.0, 0.0, 0.0);
    auto frame = coord.local_frame();
    
    // Frame vectors should be unit vectors
    EXPECT_NEAR(frame.east.norm(), 1.0, 1e-10);
    EXPECT_NEAR(frame.north.norm(), 1.0, 1e-10);
    EXPECT_NEAR(frame.up.norm(), 1.0, 1e-10);
    
    // Should be orthogonal
    EXPECT_NEAR(frame.east.dot(frame.north), 0.0, 1e-10);
    EXPECT_NEAR(frame.east.dot(frame.up), 0.0, 1e-10);
    EXPECT_NEAR(frame.north.dot(frame.up), 0.0, 1e-10);
}

TEST(Coverage100Test, NEDFrame_AtPole)
{
    // Test NED frame at South Pole
    auto coord = ECEFCoordinate::from_geodetic_deg(-90.0, 0.0, 0.0);
    auto frame = coord.local_ned_frame();
    
    // Frame vectors should be unit vectors
    EXPECT_NEAR(frame.north.norm(), 1.0, 1e-10);
    EXPECT_NEAR(frame.east.norm(), 1.0, 1e-10);
    EXPECT_NEAR(frame.down.norm(), 1.0, 1e-10);
    
    // Should be orthogonal
    EXPECT_NEAR(frame.north.dot(frame.east), 0.0, 1e-10);
    EXPECT_NEAR(frame.north.dot(frame.down), 0.0, 1e-10);
    EXPECT_NEAR(frame.east.dot(frame.down), 0.0, 1e-10);
}

TEST(Coverage100Test, ENUFrame_AtEquatorDateLine)
{
    // Test at equator on date line
    auto coord = ECEFCoordinate::from_geodetic_deg(0.0, 180.0, 0.0);
    auto frame = coord.local_frame();
    
    // Should still have valid orthonormal frame
    EXPECT_NEAR(frame.east.norm(), 1.0, 1e-10);
    EXPECT_NEAR(frame.north.norm(), 1.0, 1e-10);
    EXPECT_NEAR(frame.up.norm(), 1.0, 1e-10);
}

// ==================== Move by Bearing Edge Cases ====================

TEST(Coverage100Test, MoveByBearing_ZeroDistance)
{
    auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
    
    auto result = origin.move_by_bearing(0.0, 0.0);
    
    // Should stay at same position
    EXPECT_NEAR(result.latitude_deg(), origin.latitude_deg(), 0.0001);
    EXPECT_NEAR(result.longitude_deg(), origin.longitude_deg(), 0.0001);
}

TEST(Coverage100Test, MoveByBearing_AllBearings)
{
    auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
    double distance = 10000.0;  // 10km
    
    // Test every 45 degrees
    for (int i = 0; i < 8; ++i)
    {
        double bearing = i * PI / 4.0;
        auto result = origin.move_by_bearing(distance, bearing);
        
        // Should have moved approximately the right distance
        double actual_dist = origin.geodesic_distance_to(result);
        EXPECT_NEAR(actual_dist, distance, 100.0)
            << "Distance mismatch for bearing " << (bearing * RAD_TO_DEG) << "°";
    }
}

TEST(Coverage100Test, MoveByBearingAccurate_MaxIterations)
{
    auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
    
    // Test with different iteration counts
    auto result_5 = origin.move_by_bearing_accurate(100000.0, 0.0, 5);
    auto result_20 = origin.move_by_bearing_accurate(100000.0, 0.0, 20);
    auto result_50 = origin.move_by_bearing_accurate(100000.0, 0.0, 50);
    
    // More iterations should converge to same answer
    double diff_5_20 = result_5.distance_to(result_20);
    double diff_20_50 = result_20.distance_to(result_50);
    
    // Error should decrease with more iterations (or stay at zero if already converged)
    EXPECT_LE(diff_20_50, diff_5_20)
        << "More iterations should give better or equal convergence";
    
    // With 20 iterations, should be well converged
    EXPECT_LT(diff_20_50, 1.0)
        << "20 vs 50 iterations should differ by < 1m";
}

// ==================== Geodesic Distance Edge Cases ====================

TEST(Coverage100Test, GeodesicDistance_MaxIterations)
{
    auto coord1 = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
    auto coord2 = ECEFCoordinate::from_geodetic_deg(40.0, -118.0, 0.0);
    
    // Test with different iteration counts
    double dist_5 = coord1.geodesic_distance_to(coord2, 5);
    double dist_20 = coord1.geodesic_distance_to(coord2, 20);
    double dist_50 = coord1.geodesic_distance_to(coord2, 50);
    
    // Should converge with more iterations
    double error_5_20 = std::abs(dist_5 - dist_20);
    double error_20_50 = std::abs(dist_20 - dist_50);
    
    EXPECT_LE(error_20_50, error_5_20)
        << "More iterations should converge better or equal";
    
    // Should converge
    EXPECT_NEAR(dist_20, dist_50, 0.1)
        << "20 vs 50 iterations should give same distance";
}

TEST(Coverage100Test, GeodesicDistance_SamePoint)
{
    auto coord = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
    
    double dist = coord.geodesic_distance_to(coord);
    
    EXPECT_NEAR(dist, 0.0, 0.001)
        << "Distance to self should be zero";
}

// ==================== Trajectory Edge Cases ====================

TEST(Coverage100Test, Trajectory_EmptyLegs)
{
    auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 1000.0);
    
    std::vector<Vec3> empty_legs;
    auto result = origin.apply_enu_trajectory(empty_legs);
    
    // Should stay at origin
    EXPECT_NEAR(result.latitude_deg(), origin.latitude_deg(), 0.0001);
    EXPECT_NEAR(result.longitude_deg(), origin.longitude_deg(), 0.0001);
}

TEST(Coverage100Test, Trajectory_SingleLeg)
{
    auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 1000.0);
    
    std::vector<Vec3> single_leg = {{10000.0, 10000.0, 100.0}};
    
    auto result_trajectory = origin.apply_enu_trajectory(single_leg);
    auto result_direct = origin.apply_enu_displacement(single_leg[0]);
    
    // Should be same as direct displacement
    EXPECT_NEAR(result_trajectory.latitude_deg(), result_direct.latitude_deg(), 0.001);
    EXPECT_NEAR(result_trajectory.longitude_deg(), result_direct.longitude_deg(), 0.001);
}

TEST(Coverage100Test, Trajectory_NED_EmptyLegs)
{
    auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 1000.0);
    
    std::vector<Vec3> empty_legs;
    auto result = origin.apply_ned_trajectory(empty_legs);
    
    // Should stay at origin
    EXPECT_NEAR(result.latitude_deg(), origin.latitude_deg(), 0.0001);
    EXPECT_NEAR(result.longitude_deg(), origin.longitude_deg(), 0.0001);
}

// ==================== RK4 Edge Cases ====================

TEST(Coverage100Test, RK4_VeryFewSteps)
{
    auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
    
    // Test with minimum steps (should fall back or handle gracefully)
    auto result_1 = origin.move_geodesic_rk4(0.0, 10000.0, nullptr, 1);
    auto result_2 = origin.move_geodesic_rk4(0.0, 10000.0, nullptr, 2);
    
    // Should still give reasonable results
    double dist_1 = origin.geodesic_distance_to(result_1);
    double dist_2 = origin.geodesic_distance_to(result_2);
    
    EXPECT_GT(dist_1, 9000.0);
    EXPECT_LT(dist_1, 11000.0);
    EXPECT_GT(dist_2, 9000.0);
    EXPECT_LT(dist_2, 11000.0);
}

TEST(Coverage100Test, RK4_NegativeAltitudeProfile)
{
    auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
    
    // Profile that goes deep underwater
    auto profile = [](double d) { return -d * 0.1; };  // Descend 100m per km
    
    auto result = origin.move_geodesic_rk4(0.0, 50000.0, profile, 32);
    
    // Should end deep underwater
    EXPECT_LT(result.altitude(), -4000.0);
}

// ==================== Free Functions ====================

TEST(Coverage100Test, FreeFunction_Midpoint)
{
    auto coord1 = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 100.0);
    auto coord2 = ECEFCoordinate::from_geodetic_deg(38.0, -121.0, 200.0);
    
    auto mid = midpoint(coord1, coord2);
    
    // Should be between the two points
    EXPECT_GT(mid.latitude_deg(), coord1.latitude_deg());
    EXPECT_LT(mid.latitude_deg(), coord2.latitude_deg());
    EXPECT_GT(mid.longitude_deg(), coord1.longitude_deg());
    EXPECT_LT(mid.longitude_deg(), coord2.longitude_deg());
}

TEST(Coverage100Test, FreeFunction_Interpolate)
{
    auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 100.0);
    auto end = ECEFCoordinate::from_geodetic_deg(38.0, -121.0, 200.0);
    
    // Test various interpolation factors
    auto at_0 = interpolate(start, end, 0.0);
    auto at_25 = interpolate(start, end, 0.25);
    auto at_50 = interpolate(start, end, 0.5);
    auto at_75 = interpolate(start, end, 0.75);
    auto at_100 = interpolate(start, end, 1.0);
    
    // Endpoints should match
    EXPECT_NEAR(at_0.latitude_deg(), start.latitude_deg(), 0.001);
    EXPECT_NEAR(at_100.latitude_deg(), end.latitude_deg(), 0.001);
    
    // Midpoint should be between
    EXPECT_GT(at_50.latitude_deg(), start.latitude_deg());
    EXPECT_LT(at_50.latitude_deg(), end.latitude_deg());
}

TEST(Coverage100Test, FreeFunction_Centroid)
{
    std::vector<ECEFCoordinate> coords = {
        ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 100.0),
        ECEFCoordinate::from_geodetic_deg(37.1, -122.0, 100.0),
        ECEFCoordinate::from_geodetic_deg(37.0, -121.9, 100.0)
    };
    
    auto center = centroid(coords);
    
    // Should be near the geometric center
    EXPECT_NEAR(center.latitude_deg(), 37.033, 0.01);
    EXPECT_NEAR(center.longitude_deg(), -121.967, 0.01);
}

TEST(Coverage100Test, FreeFunction_Centroid_SinglePoint)
{
    std::vector<ECEFCoordinate> coords = {
        ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 100.0)
    };
    
    auto center = centroid(coords);
    
    // Should be the same point
    EXPECT_NEAR(center.latitude_deg(), 37.0, 0.001);
    EXPECT_NEAR(center.longitude_deg(), -122.0, 0.001);
}

// ==================== String Conversions ====================

TEST(Coverage100Test, ToString_Radians)
{
    auto coord = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 10.0);
    
    std::string str = coord.to_string(6, false);  // Radians
    
    // Should contain "rad"
    EXPECT_NE(str.find("rad"), std::string::npos);
}

TEST(Coverage100Test, ToString_DifferentPrecision)
{
    auto coord = ECEFCoordinate::from_geodetic_deg(37.123456789, -122.987654321, 10.0);
    
    std::string str_2 = coord.to_string(2, true);
    std::string str_8 = coord.to_string(8, true);
    
    // Higher precision should have longer string
    EXPECT_GT(str_8.length(), str_2.length());
}

TEST(Coverage100Test, ToStringECEF_DifferentPrecision)
{
    auto coord = ECEFCoordinate(1234567.89, 2345678.91, 3456789.12);
    
    std::string str_1 = coord.to_string_ecef(1);
    std::string str_5 = coord.to_string_ecef(5);
    
    // Higher precision should have longer string
    EXPECT_GT(str_5.length(), str_1.length());
}

// ==================== Final Edge Cases ====================

TEST(Coverage100Test, VeryHighLatitude)
{
    // Test at 89.999° (very close to pole)
    auto coord = ECEFCoordinate::from_geodetic_deg(89.999, 0.0, 0.0);
    
    // Should handle gracefully
    EXPECT_NEAR(coord.latitude_deg(), 89.999, 0.001);
    
    // Should be able to compute radii
    double M = coord.meridian_radius();
    double N = coord.prime_vertical_radius();
    
    EXPECT_GT(M, 6000000.0);
    EXPECT_GT(N, 6000000.0);
}

TEST(Coverage100Test, LongitudeWrapAround)
{
    // Test longitude wrap-around
    auto coord1 = ECEFCoordinate::from_geodetic_deg(0.0, 179.9, 0.0);
    auto coord2 = ECEFCoordinate::from_geodetic_deg(0.0, -179.9, 0.0);
    
    // Should be close together (0.2° apart)
    double dist = coord1.geodesic_distance_to(coord2);
    EXPECT_LT(dist, 25000.0);  // < 25km
}

TEST(Coverage100Test, CacheInvalidation_AfterModification)
{
    auto coord = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 1000.0);
    
    // Access geodetic to populate cache
    double lat_before = coord.latitude();
    EXPECT_TRUE(coord.has_geodetic_cache());
    
    // Modify ECEF
    coord += Vec3{10000.0, 0.0, 0.0};
    
    // Cache should still exist (it's mutable)
    // But values should be different
    double lat_after = coord.latitude();
    EXPECT_NE(lat_before, lat_after);
}
