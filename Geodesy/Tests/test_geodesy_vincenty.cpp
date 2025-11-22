// test_geodesy_vincenty.cpp - Vincenty formula accuracy tests
// Tests geodesic distance and direct formula against known solutions

#include <dstl/Geodesy.h>
#include <gtest/gtest.h>
#include "test_reference_data.h"
#include <cmath>

using namespace dstl::geo;
using namespace dstl::geo::test_data;
using namespace dstl::linalg;

// ==================== Known Airport Distances ====================

TEST(VincentyTest, KnownAirportDistances_AllPairs)
{
   // Test all known airport distance pairs from reference data
   auto distances = get_known_distances();
   auto airports = get_airports();
   
   for (const auto& ref : distances)
   {
      // Find the airports
      const Airport *from_ap = nullptr, *to_ap = nullptr;
      for (const auto& ap : airports)
      {
         if (ap.code == ref.from_code) from_ap = &ap;
         if (ap.code == ref.to_code) to_ap = &ap;
      }
      
      if (!from_ap || !to_ap) continue;
      
      auto from = ECEFCoordinate::from_geodetic_deg(from_ap->lat_deg, from_ap->lon_deg, from_ap->elevation_m);
      auto to = ECEFCoordinate::from_geodetic_deg(to_ap->lat_deg, to_ap->lon_deg, to_ap->elevation_m);
      
      double distance = from.geodesic_distance_to(to);
      
      // Use distance-dependent tolerance: 0.01% of distance, minimum 10m
      double tolerance = std::max(10.0, ref.distance_m * 0.0001);
      
      EXPECT_NEAR(distance, ref.distance_m, tolerance)
          << ref.from_code << " to " << ref.to_code << " distance incorrect";
   }
}

// ==================== Vincenty's Original Test Cases ====================

TEST(VincentyTest, VincentyOriginal_Equatorial)
{
   // Equatorial path: (0°, 0°) to (0°, 90°)
   auto start = ECEFCoordinate::from_geodetic_deg(0.0, 0.0, 0.0);
   auto end = ECEFCoordinate::from_geodetic_deg(0.0, 90.0, 0.0);
   
   double distance = start.geodesic_distance_to(end);
   
   // Quarter of Earth's circumference at equator
   double expected = PI * WGS84::a / 2.0;  // ~10,018,754 m
   
   EXPECT_NEAR(distance, expected, 100.0)
       << "Equatorial quarter-circle distance incorrect";
}

TEST(VincentyTest, VincentyOriginal_Meridional)
{
   // Meridional path: (0°, 0°) to (45°, 0°)
   auto start = ECEFCoordinate::from_geodetic_deg(0.0, 0.0, 0.0);
   auto end = ECEFCoordinate::from_geodetic_deg(45.0, 0.0, 0.0);
   
   double distance = start.geodesic_distance_to(end);
   
   // Expected: ~4,984,944 m (from Vincenty 1975)
   EXPECT_NEAR(distance, 4984944.0, 100.0)
       << "Meridional path distance incorrect";
}

TEST(VincentyTest, VincentyOriginal_Diagonal)
{
   // Diagonal path: (30°, 30°) to (60°, 60°)
   auto start = ECEFCoordinate::from_geodetic_deg(30.0, 30.0, 0.0);
   auto end = ECEFCoordinate::from_geodetic_deg(60.0, 60.0, 0.0);
   
   double distance = start.geodesic_distance_to(end);
   
   // Expected: ~4,015,703 m
   EXPECT_NEAR(distance, 4015703.0, 100.0)
       << "Diagonal path distance incorrect";
}

// ==================== Vincenty Direct Formula Tests ====================

TEST(VincentyTest, DirectFormula_ShortDistance)
{
   auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
   
   // Move 10km due North
   double distance = 10000.0;
   double bearing = 0.0;  // North
   
   auto end = start.move_by_bearing_accurate(distance, bearing);
   
   // Verify distance back
   double dist_back = start.geodesic_distance_to(end);
   EXPECT_NEAR(dist_back, distance, 1.0)  // Within 1m
       << "Direct formula round-trip failed for 10km";
   
   // Verify bearing back (handle wrap-around at 0°/360°)
   double bearing_back = start.bearing_to(end);
   // Normalize both to [0, 2π)
   double bearing_norm = std::fmod(bearing + 2.0 * PI, 2.0 * PI);
   double bearing_back_norm = std::fmod(bearing_back + 2.0 * PI, 2.0 * PI);
   // Check if they're close, accounting for wrap-around
   double diff = std::abs(bearing_norm - bearing_back_norm);
   if (diff > PI) diff = 2.0 * PI - diff;
   EXPECT_LT(diff, 0.001)  // Within ~0.06°
       << "Bearing round-trip failed";
}

TEST(VincentyTest, DirectFormula_MediumDistance)
{
   auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
   
   // Move 500km at 45° (NE)
   double distance = 500000.0;
   double bearing = PI / 4.0;  // 45° NE
   
   auto end = start.move_by_bearing_accurate(distance, bearing);
   
   // Verify distance back
   double dist_back = start.geodesic_distance_to(end);
   EXPECT_NEAR(dist_back, distance, 10.0)  // Within 10m
       << "Direct formula round-trip failed for 500km";
}

TEST(VincentyTest, DirectFormula_LongDistance)
{
   auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
   
   // Move 5000km at 90° (East)
   double distance = 5000000.0;
   double bearing = PI / 2.0;  // East
   
   auto end = start.move_by_bearing_accurate(distance, bearing);
   
   // Verify distance back
   double dist_back = start.geodesic_distance_to(end);
   EXPECT_NEAR(dist_back, distance, 100.0)  // Within 100m
       << "Direct formula round-trip failed for 5000km";
}

TEST(VincentyTest, DirectFormula_AllCardinalDirections)
{
   auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
   double distance = 100000.0;  // 100km
   
   // Test all cardinal directions
   std::vector<std::pair<std::string, double>> directions = {
       {"North", 0.0},
       {"East", PI / 2.0},
       {"South", PI},
       {"West", 3.0 * PI / 2.0}
   };
   
   for (const auto& [name, bearing] : directions)
   {
      auto target = origin.move_by_bearing_accurate(distance, bearing);
      double dist_back = origin.geodesic_distance_to(target);
      
      EXPECT_NEAR(dist_back, distance, 10.0)
          << "Direct formula failed for " << name << " direction";
   }
}

// ==================== Convergence Behavior Tests ====================

TEST(VincentyTest, Convergence_NormalCase)
{
   // Test that Vincenty converges quickly for normal cases
   auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
   auto end = ECEFCoordinate::from_geodetic_deg(40.0, -118.0, 0.0);
   
   // Should converge in < 10 iterations
   double distance = start.geodesic_distance_to(end);
   
   EXPECT_GT(distance, 0.0);
   EXPECT_LT(distance, 20000000.0);  // Sanity check
}

TEST(VincentyTest, Convergence_NearPole)
{
   // Test convergence near pole
   auto start = ECEFCoordinate::from_geodetic_deg(89.0, 0.0, 0.0);
   auto end = ECEFCoordinate::from_geodetic_deg(89.0, 90.0, 0.0);
   
   double distance = start.geodesic_distance_to(end);
   
   EXPECT_GT(distance, 0.0);
   EXPECT_LT(distance, 200000.0);  // Should be < 200km
}

TEST(VincentyTest, Convergence_AcrossDateLine)
{
   // Test across international date line
   auto start = ECEFCoordinate::from_geodetic_deg(0.0, 179.0, 0.0);
   auto end = ECEFCoordinate::from_geodetic_deg(0.0, -179.0, 0.0);
   
   double distance = start.geodesic_distance_to(end);
   
   // Should be ~222 km (2° at equator)
   EXPECT_NEAR(distance, 222000.0, 1000.0);
}

// ==================== Symmetry Tests ====================

TEST(VincentyTest, Symmetry_DistanceCommutative)
{
   auto point1 = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
   auto point2 = ECEFCoordinate::from_geodetic_deg(40.0, -118.0, 0.0);
   
   double dist12 = point1.geodesic_distance_to(point2);
   double dist21 = point2.geodesic_distance_to(point1);
   
   EXPECT_NEAR(dist12, dist21, 0.001)
       << "Geodesic distance should be symmetric";
}

TEST(VincentyTest, Symmetry_BearingReciprocal)
{
   auto point1 = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
   auto point2 = ECEFCoordinate::from_geodetic_deg(40.0, -118.0, 0.0);
   
   double bearing12 = point1.bearing_to(point2);
   double bearing21 = point2.bearing_to(point1);
   
   // Reciprocal bearing should differ by π (180°)
   // Normalize the difference
   double diff = bearing12 - bearing21;
   // Normalize to [-π, π]
   while (diff > PI) diff -= 2.0 * PI;
   while (diff < -PI) diff += 2.0 * PI;
   diff = std::abs(diff);
   
   // Should be close to π
   EXPECT_NEAR(diff, PI, 0.05)  // Within ~3° (bearings change slightly over distance)
       << "Bearings should be approximately reciprocal";
}

// ==================== Altitude Independence Tests ====================

TEST(VincentyTest, AltitudeIndependence_GeodesicDistance)
{
   // Geodesic distance should be independent of altitude
   // (it's measured on the ellipsoid surface)
   
   auto point1_sea = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
   auto point2_sea = ECEFCoordinate::from_geodetic_deg(40.0, -118.0, 0.0);
   
   auto point1_high = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 10000.0);
   auto point2_high = ECEFCoordinate::from_geodetic_deg(40.0, -118.0, 10000.0);
   
   double dist_sea = point1_sea.geodesic_distance_to(point2_sea);
   double dist_high = point1_high.geodesic_distance_to(point2_high);
   
   // Should be nearly identical (within numerical precision)
   EXPECT_NEAR(dist_sea, dist_high, 1.0)
       << "Geodesic distance should be altitude-independent";
}

TEST(VincentyTest, AltitudePreservation_DirectFormula)
{
   auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 5000.0);
   
   // Move 100km North
   auto end = start.move_by_bearing_accurate(100000.0, 0.0);
   
   // Altitude should be preserved
   EXPECT_NEAR(end.altitude(), start.altitude(), 1.0)
       << "Direct formula should preserve altitude";
}

// ==================== Edge Cases ====================

TEST(VincentyTest, EdgeCase_ZeroDistance)
{
   auto point = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
   
   double distance = point.geodesic_distance_to(point);
   
   EXPECT_NEAR(distance, 0.0, 0.001)
       << "Distance to self should be zero";
}

TEST(VincentyTest, EdgeCase_VeryShortDistance)
{
   auto point1 = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
   auto point2 = ECEFCoordinate::from_geodetic_deg(37.0001, -122.0, 0.0);
   
   double distance = point1.geodesic_distance_to(point2);
   
   // ~11 meters (0.0001° at this latitude)
   EXPECT_NEAR(distance, 11.0, 1.0);
}

TEST(VincentyTest, EdgeCase_EquatorCrossing)
{
   auto north = ECEFCoordinate::from_geodetic_deg(10.0, 0.0, 0.0);
   auto south = ECEFCoordinate::from_geodetic_deg(-10.0, 0.0, 0.0);
   
   double distance = north.geodesic_distance_to(south);
   
   // 20° along meridian
   EXPECT_GT(distance, 2200000.0);  // > 2200 km
   EXPECT_LT(distance, 2300000.0);  // < 2300 km
}

TEST(VincentyTest, EdgeCase_PrimeMeridianCrossing)
{
   auto west = ECEFCoordinate::from_geodetic_deg(45.0, -10.0, 0.0);
   auto east = ECEFCoordinate::from_geodetic_deg(45.0, 10.0, 0.0);
   
   double distance = west.geodesic_distance_to(east);
   
   // 20° at 45° latitude
   EXPECT_GT(distance, 1500000.0);  // > 1500 km
   EXPECT_LT(distance, 1600000.0);  // < 1600 km
}

// ==================== Accuracy vs Simple Methods ====================

TEST(VincentyTest, Accuracy_VsChordDistance)
{
   auto point1 = ECEFCoordinate::from_geodetic_deg(0.0, 0.0, 0.0);
   auto point2 = ECEFCoordinate::from_geodetic_deg(0.0, 90.0, 0.0);
   
   double geodesic = point1.geodesic_distance_to(point2);
   double chord = point1.distance_to(point2);
   
   // Geodesic should be longer than chord
   EXPECT_GT(geodesic, chord);
   
   // For quarter circle, difference should be significant
   double diff_percent = (geodesic - chord) / geodesic * 100.0;
   EXPECT_GT(diff_percent, 9.5)  // > 9.5% difference (actual is ~9.97%)
       << "Geodesic should be significantly longer than chord for long distances";
}

TEST(VincentyTest, Accuracy_ShortDistanceConvergence)
{
   // For very short distances, geodesic and chord should converge
   auto point1 = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
   auto point2 = ECEFCoordinate::from_geodetic_deg(37.001, -122.0, 0.0);
   
   double geodesic = point1.geodesic_distance_to(point2);
   double chord = point1.distance_to(point2);
   
   // Should differ by < 0.01%
   double diff_percent = std::abs(geodesic - chord) / geodesic * 100.0;
   EXPECT_LT(diff_percent, 0.01)
       << "Geodesic and chord should be nearly equal for short distances";
}

// ==================== Consistency Tests ====================

TEST(VincentyTest, Consistency_InverseDirectRoundTrip)
{
   auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
   auto end = ECEFCoordinate::from_geodetic_deg(40.0, -118.0, 0.0);
   
   // Inverse: compute distance and bearing
   double distance = start.geodesic_distance_to(end);
   double bearing = start.bearing_to(end);
   
   // Direct: use distance and bearing to compute endpoint
   auto end_computed = start.move_by_bearing_accurate(distance, bearing);
   
   // Should arrive at same point
   EXPECT_NEAR(end_computed.latitude_deg(), end.latitude_deg(), 0.001)
       << "Inverse-direct round-trip latitude mismatch";
   EXPECT_NEAR(end_computed.longitude_deg(), end.longitude_deg(), 0.001)
       << "Inverse-direct round-trip longitude mismatch";
}

TEST(VincentyTest, Consistency_MultipleHops)
{
   auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
   
   // Move 100km North in one hop
   auto end_single = start.move_by_bearing_accurate(100000.0, 0.0);
   
   // Move 100km North in two 50km hops
   auto mid = start.move_by_bearing_accurate(50000.0, 0.0);
   auto end_double = mid.move_by_bearing_accurate(50000.0, 0.0);
   
   // Should arrive at approximately same point
   EXPECT_NEAR(end_single.latitude_deg(), end_double.latitude_deg(), 0.01)
       << "Multiple hops should give consistent result";
   EXPECT_NEAR(end_single.longitude_deg(), end_double.longitude_deg(), 0.01)
       << "Multiple hops should give consistent result";
}

// ==================== Performance Characteristics ====================

TEST(VincentyTest, Performance_ReasonableIterations)
{
   // Verify that Vincenty converges in reasonable iterations
   // (This is more of a sanity check than a performance test)
   
   auto point1 = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
   auto point2 = ECEFCoordinate::from_geodetic_deg(51.5, -0.1, 0.0);
   
   // Should complete without hanging
   double distance = point1.geodesic_distance_to(point2);
   
   EXPECT_GT(distance, 8000000.0);  // > 8000 km
   EXPECT_LT(distance, 9000000.0);  // < 9000 km
}

TEST(VincentyTest, Performance_DirectFormulaSpeed)
{
   auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
   
   // Should complete quickly even for long distances
   auto end = start.move_by_bearing_accurate(10000000.0, PI / 4.0);
   
   // Just verify it completes and gives reasonable result
   EXPECT_NE(end.latitude_deg(), start.latitude_deg());
   EXPECT_NE(end.longitude_deg(), start.longitude_deg());
}
