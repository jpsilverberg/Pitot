// test_geodesy_complete.cpp - Complete coverage tests
// Tests for spherical methods, English units, and remaining functions

#include <dstl/Geodesy.h>
#include <gtest/gtest.h>
#include "test_reference_data.h"
#include <cmath>

using namespace dstl::geo;
using namespace dstl::geo::test_data;
using namespace dstl::linalg;

// ==================== English Units - Construction ====================

TEST(CompleteTest, Construction_FromGeodeticFeet)
{
    // Test construction with altitude in feet
    double lat_deg = 37.0;
    double lon_deg = -122.0;
    double alt_ft = 32808.4;  // ~10,000m
    
    auto coord = ECEFCoordinate::from_geodetic_deg_ft(lat_deg, lon_deg, alt_ft);
    
    // Verify conversion
    EXPECT_NEAR(coord.latitude_deg(), lat_deg, 0.0001);
    EXPECT_NEAR(coord.longitude_deg(), lon_deg, 0.0001);
    EXPECT_NEAR(coord.altitude(), alt_ft * FT_TO_M, 1.0);
}

TEST(CompleteTest, Construction_FeetRoundTrip)
{
    double alt_ft = 39000.0;  // Typical cruise altitude
    
    auto coord = ECEFCoordinate::from_geodetic_deg_ft(37.0, -122.0, alt_ft);
    double alt_back_ft = coord.altitude_ft();
    
    EXPECT_NEAR(alt_back_ft, alt_ft, 0.1)
        << "Feet round-trip should be accurate";
}

// ==================== English Units - Altitude ====================

TEST(CompleteTest, AltitudeFeet_Conversion)
{
    auto coord = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 10000.0);
    
    double alt_ft = coord.altitude_ft();
    double expected_ft = 10000.0 * M_TO_FT;
    
    EXPECT_NEAR(alt_ft, expected_ft, 0.1)
        << "Altitude in feet should match conversion";
}

TEST(CompleteTest, AltitudeFeet_SeaLevel)
{
    auto coord = ECEFCoordinate::from_geodetic_deg(0.0, 0.0, 0.0);
    
    EXPECT_NEAR(coord.altitude_ft(), 0.0, 0.1)
        << "Sea level should be 0 feet";
}

TEST(CompleteTest, AltitudeFeet_Negative)
{
    // Dead Sea is ~430m below sea level
    auto coord = ECEFCoordinate::from_geodetic_deg(31.5, 35.5, -430.0);
    
    double alt_ft = coord.altitude_ft();
    EXPECT_LT(alt_ft, 0.0)
        << "Below sea level should be negative feet";
    EXPECT_NEAR(alt_ft, -430.0 * M_TO_FT, 1.0);
}

// ==================== English Units - Displacements ====================

TEST(CompleteTest, ENUDisplacement_ConstantAltitudeFeet)
{
    auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 10000.0);
    
    // Move 10km E, 10km N, climb 1000 ft
    Vec3 enu{10000.0, 10000.0, 0.0};
    double alt_change_ft = 1000.0;
    
    auto target = origin.apply_enu_displacement_constant_altitude_ft(enu, alt_change_ft);
    
    // Verify altitude change
    double alt_change_m = target.altitude() - origin.altitude();
    EXPECT_NEAR(alt_change_m, alt_change_ft * FT_TO_M, 1.0)
        << "Altitude change should match feet conversion";
}

TEST(CompleteTest, ENUDisplacement_FeetZeroAltitudeChange)
{
    auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 5000.0);
    
    Vec3 enu{5000.0, 5000.0, 0.0};
    auto target = origin.apply_enu_displacement_constant_altitude_ft(enu, 0.0);
    
    // Altitude should be approximately preserved
    EXPECT_NEAR(target.altitude(), origin.altitude(), 10.0);
}

// ==================== Spherical Displacement Methods ====================

TEST(CompleteTest, SphericalDisplacement_ShortDistance)
{
    auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 1000.0);
    
    // Move 10km E, 10km N on spherical surface
    Vec3 enu{10000.0, 10000.0, 0.0};
    
    auto target = origin.apply_enu_displacement_spherical(enu, 0.0, 3);
    
    // Should have moved
    EXPECT_GT(target.latitude_deg(), origin.latitude_deg());
    EXPECT_GT(target.longitude_deg(), origin.longitude_deg());
    
    // Altitude should be approximately preserved
    EXPECT_NEAR(target.altitude(), origin.altitude(), 10.0);
}

TEST(CompleteTest, SphericalDisplacement_WithAltitudeChange)
{
    auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 5000.0);
    
    // Move horizontally and climb 2000m
    Vec3 enu{20000.0, 0.0, 0.0};  // 20km East
    double alt_change = 2000.0;
    
    auto target = origin.apply_enu_displacement_spherical(enu, alt_change, 3);
    
    // Should have moved East
    EXPECT_GT(target.longitude_deg(), origin.longitude_deg());
    
    // Should have climbed
    double actual_alt_change = target.altitude() - origin.altitude();
    EXPECT_NEAR(actual_alt_change, alt_change, 50.0)
        << "Altitude change should be approximately correct";
}

TEST(CompleteTest, SphericalDisplacement_IterativeRefinement)
{
    auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 10000.0);
    
    // Test with different iteration counts
    Vec3 enu{50000.0, 50000.0, 0.0};  // 50km diagonal
    
    auto result_1 = origin.apply_enu_displacement_spherical(enu, 0.0, 1);
    auto result_3 = origin.apply_enu_displacement_spherical(enu, 0.0, 3);
    auto result_5 = origin.apply_enu_displacement_spherical(enu, 0.0, 5);
    
    // More iterations should converge (or at least not diverge significantly)
    double error_1_3 = result_1.distance_to(result_3);
    double error_3_5 = result_3.distance_to(result_5);
    
    // With more iterations, error should be small
    EXPECT_LT(error_3_5, 100.0)
        << "With 5 iterations, should be well converged";
}

TEST(CompleteTest, SphericalDisplacement_VsFlatEarth)
{
    auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
    
    // For short distances, spherical and flat-earth should be similar
    Vec3 enu{5000.0, 5000.0, 0.0};  // 5km diagonal
    
    auto flat = origin.apply_enu_displacement(enu);
    auto spherical = origin.apply_enu_displacement_spherical(enu, 0.0, 3);
    
    // Should be reasonably close for short distances
    double diff = flat.distance_to(spherical);
    EXPECT_LT(diff, 50.0)
        << "Flat-earth and spherical should be similar for short distances";
}

TEST(CompleteTest, SphericalDisplacement_LongDistance)
{
    auto origin = ECEFCoordinate::from_geodetic_deg(0.0, 0.0, 0.0);
    
    // Long distance where spherical is better than flat-earth
    Vec3 enu{500000.0, 0.0, 0.0};  // 500km East
    
    auto spherical = origin.apply_enu_displacement_spherical(enu, 0.0, 5);
    
    // Should have moved significantly
    EXPECT_GT(spherical.longitude_deg(), 4.0)
        << "Should have moved >4° East";
}

TEST(CompleteTest, SphericalDisplacement_Feet)
{
    auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 10000.0);
    
    // Move with altitude change in feet
    Vec3 enu{10000.0, 10000.0, 0.0};
    double alt_change_ft = 5000.0;
    
    auto target = origin.apply_enu_displacement_spherical_ft(enu, alt_change_ft, 3);
    
    // Verify altitude change
    double alt_change_m = target.altitude() - origin.altitude();
    EXPECT_NEAR(alt_change_m, alt_change_ft * FT_TO_M, 50.0);
}

// ==================== Unit Conversions ====================

TEST(CompleteTest, UnitConversions_FeetToMeters)
{
    auto conversions = get_unit_conversions();
    
    for (const auto& conv : conversions)
    {
        double meters = conv.feet * FT_TO_M;
        EXPECT_NEAR(meters, conv.meters, 0.01)
            << "Feet to meters conversion";
    }
}

TEST(CompleteTest, UnitConversions_NauticalMiles)
{
    auto conversions = get_unit_conversions();
    
    for (const auto& conv : conversions)
    {
        double meters = conv.nautical_miles * NM_TO_M;
        EXPECT_NEAR(meters, conv.meters, 0.1)
            << "Nautical miles to meters conversion";
    }
}

TEST(CompleteTest, UnitConversions_Knots)
{
    auto conversions = get_unit_conversions();
    
    for (const auto& conv : conversions)
    {
        double ms = conv.knots * KT_TO_MS;
        EXPECT_NEAR(ms, conv.meters_per_sec, 0.01)
            << "Knots to m/s conversion";
    }
}

// ==================== Trajectory Frame Updates ====================

TEST(CompleteTest, Trajectory_FrameUpdateInterval)
{
    auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 1000.0);
    
    // Create trajectory with multiple legs
    std::vector<Vec3> legs;
    for (int i = 0; i < 10; ++i)
    {
        legs.push_back({5000.0, 5000.0, 100.0});  // 10 legs of 5km each
    }
    
    // Test with different update intervals
    auto result_no_update = origin.apply_enu_trajectory(legs, 0.0);  // Single frame
    auto result_10km = origin.apply_enu_trajectory(legs, 10000.0);   // Update every 10km
    auto result_5km = origin.apply_enu_trajectory(legs, 5000.0);     // Update every 5km
    
    // All should reach similar final position
    double diff_no_10 = result_no_update.distance_to(result_10km);
    double diff_10_5 = result_10km.distance_to(result_5km);
    
    // More frequent updates should be more accurate
    EXPECT_LT(diff_10_5, 100.0)
        << "Different update intervals should give similar results";
}

TEST(CompleteTest, Trajectory_SingleFrameVsMultiple)
{
    auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 1000.0);
    
    std::vector<Vec3> legs = {
        {10000.0, 0.0, 0.0},
        {0.0, 10000.0, 0.0},
        {-5000.0, 5000.0, 0.0}
    };
    
    // Single frame (update_interval <= 0)
    auto single = origin.apply_enu_trajectory(legs, 0.0);
    
    // Multiple frames
    auto multiple = origin.apply_enu_trajectory(legs, 5000.0);
    
    // Should be reasonably close
    double diff = single.distance_to(multiple);
    EXPECT_LT(diff, 100.0)
        << "Single vs multiple frames should be close";
}

TEST(CompleteTest, Trajectory_NEDFrameUpdates)
{
    auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 1000.0);
    
    std::vector<Vec3> legs = {
        {20000.0, 10000.0, -500.0},  // 20km N, 10km E, climb 500m
        {10000.0, 20000.0, 200.0}    // 10km N, 20km E, descend 200m
    };
    
    auto result = origin.apply_ned_trajectory(legs, 10000.0);
    
    // Should have moved significantly
    double dist = origin.geodesic_distance_to(result);
    EXPECT_GT(dist, 40000.0)
        << "Should have moved >40km";
}

// ==================== Additional Coverage ====================

TEST(CompleteTest, Degrees_RadiansConversion)
{
    // Test degree/radian conversions
    double deg = 45.0;
    double rad = deg * DEG_TO_RAD;
    double deg_back = rad * RAD_TO_DEG;
    
    EXPECT_NEAR(deg_back, deg, 1e-10);
}

TEST(CompleteTest, Constants_Consistency)
{
    // Verify unit conversion constants are consistent
    EXPECT_NEAR(FT_TO_M * M_TO_FT, 1.0, 1e-10);
    EXPECT_NEAR(NM_TO_M * M_TO_NM, 1.0, 1e-10);
    EXPECT_NEAR(KT_TO_MS * MS_TO_KT, 1.0, 1e-10);
}

TEST(CompleteTest, WGS84_Constants)
{
    // Verify WGS84 constants are consistent
    double b_computed = WGS84::a * (1.0 - WGS84::f);
    EXPECT_NEAR(WGS84::b, b_computed, 1e-6);
    
    double e2_computed = WGS84::f * (2.0 - WGS84::f);
    EXPECT_NEAR(WGS84::e2, e2_computed, 1e-12);
}

// ==================== Stress Tests ====================

TEST(CompleteTest, Stress_ManySmallSteps)
{
    auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
    
    // Apply many small displacements
    auto current = origin;
    for (int i = 0; i < 100; ++i)
    {
        Vec3 enu{100.0, 100.0, 10.0};  // 100m steps
        current = current.apply_enu_displacement(enu);
    }
    
    // Should have moved ~14km total
    double dist = origin.geodesic_distance_to(current);
    EXPECT_GT(dist, 13000.0);
    EXPECT_LT(dist, 15000.0);
}

TEST(CompleteTest, Stress_HighAltitude)
{
    // Test at LEO orbit altitude
    auto coord = ECEFCoordinate::from_geodetic_deg(0.0, 0.0, 400000.0);
    
    // Should handle gracefully
    EXPECT_NEAR(coord.latitude_deg(), 0.0, 0.0001);
    EXPECT_NEAR(coord.longitude_deg(), 0.0, 0.0001);
    EXPECT_NEAR(coord.altitude(), 400000.0, 1.0);
}

TEST(CompleteTest, Stress_DeepUnderwater)
{
    // Mariana Trench depth
    auto coord = ECEFCoordinate::from_geodetic_deg(11.35, 142.2, -10994.0);
    
    // Should handle negative altitude
    EXPECT_LT(coord.altitude(), 0.0);
    EXPECT_NEAR(coord.altitude(), -10994.0, 1.0);
}

TEST(CompleteTest, Stress_RapidAltitudeChanges)
{
    auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
    
    // Rapid climb and descent
    std::vector<Vec3> legs = {
        {1000.0, 0.0, 5000.0},   // Climb 5000m
        {1000.0, 0.0, -5000.0},  // Descend 5000m
        {1000.0, 0.0, 3000.0},   // Climb 3000m
        {1000.0, 0.0, -3000.0}   // Descend 3000m
    };
    
    auto final = origin.apply_enu_trajectory(legs);
    
    // Should end near starting altitude
    EXPECT_NEAR(final.altitude(), origin.altitude(), 100.0);
}

// ==================== Integration Tests ====================

TEST(CompleteTest, Integration_CompleteFlightProfile)
{
    // Simulate a complete flight: takeoff, climb, cruise, descend, land
    auto origin = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 4.0);  // SFO
    
    // Fly to LAX (33.9416, -118.4085)
    double track_deg = 137.5;  // SE
    double distance_nm = 293.0;  // ~543 km
    
    // Altitude profile: climb to 35,000 ft, cruise, descend
    auto profile = [](double d_nm) {
        if (d_nm < 50.0) return d_nm * 700.0;           // Climb: 0 to 35,000 ft
        if (d_nm < 243.0) return 35000.0;              // Cruise
        return 35000.0 - (d_nm - 243.0) * 700.0;       // Descend
    };
    
    auto destination = origin.fly_altitude_profile(track_deg, distance_nm, profile, 64);
    
    // Should be in the general vicinity of LAX (within ~20km)
    EXPECT_NEAR(destination.latitude_deg(), 33.9416, 0.2);
    EXPECT_NEAR(destination.longitude_deg(), -118.4085, 0.2);
    
    // Should be near ground level
    EXPECT_LT(destination.altitude(), 200.0);
}

TEST(CompleteTest, Integration_AroundTheWorld)
{
    auto start = ECEFCoordinate::from_geodetic_deg(0.0, 0.0, 0.0);
    
    // Go around the world in 4 steps (90° each)
    auto east1 = start.move_by_bearing_accurate(10018754.0, PI / 2.0);  // 90° East
    auto east2 = east1.move_by_bearing_accurate(10018754.0, PI / 2.0);  // 180° total
    auto east3 = east2.move_by_bearing_accurate(10018754.0, PI / 2.0);  // 270° total
    auto back = east3.move_by_bearing_accurate(10018754.0, PI / 2.0);   // 360° total
    
    // Should be back near start
    double error = start.distance_to(back);
    EXPECT_LT(error, 1000.0)
        << "Around the world should return within 1km";
}

TEST(CompleteTest, Integration_PoleToEquator)
{
    auto north_pole = ECEFCoordinate::from_geodetic_deg(90.0, 0.0, 0.0);
    auto equator = ECEFCoordinate::from_geodetic_deg(0.0, 0.0, 0.0);
    
    // Distance should be approximately quarter of meridian
    // (WGS84 ellipsoid is shorter than sphere at poles)
    double dist = north_pole.geodesic_distance_to(equator);
    double quarter_meridian = PI * WGS84::a / 2.0;
    
    // Allow larger tolerance due to ellipsoid vs sphere difference
    EXPECT_NEAR(dist, quarter_meridian, 20000.0)
        << "Pole to equator should be approximately quarter meridian";
}
