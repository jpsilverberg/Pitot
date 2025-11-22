// test_geodesy_reference.cpp - Reference point validation tests
// Tests against mathematically known WGS84 values and reference coordinates

#include <dstl/Geodesy.h>
#include <gtest/gtest.h>
#include "test_reference_data.h"
#include <cmath>

using namespace dstl::geo;
using namespace dstl::geo::test_data;
using namespace dstl::linalg;

// ==================== WGS84 Constants Validation ====================

TEST(ReferenceTest, WGS84_Constants)
{
   // Verify WGS84 constants match reference values
   EXPECT_DOUBLE_EQ(WGS84::a, WGS84Reference::a);
   EXPECT_NEAR(WGS84::f, WGS84Reference::f, EPSILON_MICRO);
   EXPECT_NEAR(WGS84::b, WGS84Reference::b, EPSILON_MM);
   EXPECT_NEAR(WGS84::e2, WGS84Reference::e2, EPSILON_MICRO);
   EXPECT_NEAR(WGS84::ep2, WGS84Reference::ep2, EPSILON_MICRO);
}

TEST(ReferenceTest, WGS84_DerivedValues)
{
   // Verify derived relationships
   double b_computed = WGS84::a * (1.0 - WGS84::f);
   EXPECT_NEAR(WGS84::b, b_computed, EPSILON_MM);
   
   double e2_computed = WGS84::f * (2.0 - WGS84::f);
   EXPECT_NEAR(WGS84::e2, e2_computed, EPSILON_MICRO);
   
   double ep2_computed = WGS84::e2 / (1.0 - WGS84::e2);
   EXPECT_NEAR(WGS84::ep2, ep2_computed, EPSILON_MICRO);
}

// ==================== Special Points ECEF Validation ====================

TEST(ReferenceTest, EquatorPrimeMeridian_ECEF)
{
   auto coord = ECEFCoordinate::from_geodetic_deg(0.0, 0.0, 0.0);
   
   // At (0°, 0°, 0m), ECEF should be (a, 0, 0)
   EXPECT_NEAR(coord.x(), WGS84::a, EPSILON_MM);
   EXPECT_NEAR(coord.y(), 0.0, EPSILON_MM);
   EXPECT_NEAR(coord.z(), 0.0, EPSILON_MM);
}

TEST(ReferenceTest, NorthPole_ECEF)
{
   auto coord = ECEFCoordinate::from_geodetic_deg(90.0, 0.0, 0.0);
   
   // At North Pole, ECEF should be (0, 0, b)
   EXPECT_NEAR(coord.x(), 0.0, EPSILON_MM);
   EXPECT_NEAR(coord.y(), 0.0, EPSILON_MM);
   EXPECT_NEAR(coord.z(), WGS84::b, EPSILON_MM);
}

TEST(ReferenceTest, SouthPole_ECEF)
{
   auto coord = ECEFCoordinate::from_geodetic_deg(-90.0, 0.0, 0.0);
   
   // At South Pole, ECEF should be (0, 0, -b)
   EXPECT_NEAR(coord.x(), 0.0, EPSILON_MM);
   EXPECT_NEAR(coord.y(), 0.0, EPSILON_MM);
   EXPECT_NEAR(coord.z(), -WGS84::b, EPSILON_MM);
}

TEST(ReferenceTest, Equator90E_ECEF)
{
   auto coord = ECEFCoordinate::from_geodetic_deg(0.0, 90.0, 0.0);
   
   // At (0°, 90°E, 0m), ECEF should be (0, a, 0)
   EXPECT_NEAR(coord.x(), 0.0, EPSILON_MM);
   EXPECT_NEAR(coord.y(), WGS84::a, EPSILON_MM);
   EXPECT_NEAR(coord.z(), 0.0, EPSILON_MM);
}

TEST(ReferenceTest, Equator180E_ECEF)
{
   auto coord = ECEFCoordinate::from_geodetic_deg(0.0, 180.0, 0.0);
   
   // At (0°, 180°E, 0m), ECEF should be (-a, 0, 0)
   EXPECT_NEAR(coord.x(), -WGS84::a, EPSILON_MM);
   EXPECT_NEAR(coord.y(), 0.0, EPSILON_MM);
   EXPECT_NEAR(coord.z(), 0.0, EPSILON_MM);
}

TEST(ReferenceTest, Equator90W_ECEF)
{
   auto coord = ECEFCoordinate::from_geodetic_deg(0.0, -90.0, 0.0);
   
   // At (0°, 90°W, 0m), ECEF should be (0, -a, 0)
   EXPECT_NEAR(coord.x(), 0.0, EPSILON_MM);
   EXPECT_NEAR(coord.y(), -WGS84::a, EPSILON_MM);
   EXPECT_NEAR(coord.z(), 0.0, EPSILON_MM);
}

TEST(ReferenceTest, AllSpecialPoints_ECEF)
{
   auto points = get_special_points();
   
   for (const auto& pt : points)
   {
      auto coord = ECEFCoordinate::from_geodetic_deg(pt.lat_deg, pt.lon_deg, pt.alt_m);
      
      EXPECT_NEAR(coord.x(), pt.ecef_x, EPSILON_M)
          << "X mismatch for " << pt.name << " (source: " << pt.source << ")";
      EXPECT_NEAR(coord.y(), pt.ecef_y, EPSILON_M)
          << "Y mismatch for " << pt.name << " (source: " << pt.source << ")";
      EXPECT_NEAR(coord.z(), pt.ecef_z, EPSILON_M)
          << "Z mismatch for " << pt.name << " (source: " << pt.source << ")";
   }
}

// ==================== Special Points Geodetic Round-Trip ====================

TEST(ReferenceTest, SpecialPoints_GeodeticRoundTrip)
{
   auto points = get_special_points();
   
   for (const auto& pt : points)
   {
      auto coord = ECEFCoordinate::from_geodetic_deg(pt.lat_deg, pt.lon_deg, pt.alt_m);
      
      // Convert back to geodetic
      double lat_back = coord.latitude_deg();
      double lon_back = coord.longitude_deg();
      double alt_back = coord.altitude();
      
      EXPECT_NEAR(lat_back, pt.lat_deg, EPSILON_DEG)
          << "Latitude round-trip failed for " << pt.name;
      EXPECT_NEAR(lon_back, pt.lon_deg, EPSILON_DEG)
          << "Longitude round-trip failed for " << pt.name;
      EXPECT_NEAR(alt_back, pt.alt_m, EPSILON_CM)
          << "Altitude round-trip failed for " << pt.name;
   }
}

TEST(ReferenceTest, SpecialPoints_ECEFRoundTrip)
{
   auto points = get_special_points();
   
   for (const auto& pt : points)
   {
      // Start from ECEF
      auto coord = ECEFCoordinate(pt.ecef_x, pt.ecef_y, pt.ecef_z);
      
      // Convert to geodetic and back to ECEF
      auto coord2 = ECEFCoordinate::from_geodetic(
          coord.latitude(), coord.longitude(), coord.altitude());
      
      EXPECT_NEAR(coord2.x(), pt.ecef_x, EPSILON_CM)
          << "X round-trip failed for " << pt.name;
      EXPECT_NEAR(coord2.y(), pt.ecef_y, EPSILON_CM)
          << "Y round-trip failed for " << pt.name;
      EXPECT_NEAR(coord2.z(), pt.ecef_z, EPSILON_CM)
          << "Z round-trip failed for " << pt.name;
   }
}

// ==================== Radii of Curvature Validation ====================

TEST(ReferenceTest, MeridianRadius_Equator)
{
   auto coord = ECEFCoordinate::from_geodetic_deg(0.0, 0.0, 0.0);
   double M = coord.meridian_radius();
   
   // At equator: M = a(1-e²)/(1-e²sin²φ)^(3/2) = a(1-e²)
   double M_expected = WGS84::a * (1.0 - WGS84::e2);
   
   EXPECT_NEAR(M, M_expected, EPSILON_MM)
       << "Meridian radius at equator should be " << M_expected;
}

TEST(ReferenceTest, PrimeVerticalRadius_Equator)
{
   auto coord = ECEFCoordinate::from_geodetic_deg(0.0, 0.0, 0.0);
   double N = coord.prime_vertical_radius();
   
   // At equator: N = a/√(1-e²sin²φ) = a
   EXPECT_NEAR(N, WGS84::a, EPSILON_MM)
       << "Prime vertical radius at equator should be a";
}

TEST(ReferenceTest, MeridianRadius_Pole)
{
   auto coord = ECEFCoordinate::from_geodetic_deg(90.0, 0.0, 0.0);
   double M = coord.meridian_radius();
   
   // At pole: M = a/√(1-e²) = a²/b
   double M_expected = WGS84::a * WGS84::a / WGS84::b;
   
   EXPECT_NEAR(M, M_expected, EPSILON_M)
       << "Meridian radius at pole should be " << M_expected;
}

TEST(ReferenceTest, PrimeVerticalRadius_Pole)
{
   auto coord = ECEFCoordinate::from_geodetic_deg(90.0, 0.0, 0.0);
   double N = coord.prime_vertical_radius();
   
   // At pole: N = M (both radii are equal)
   double M = coord.meridian_radius();
   
   EXPECT_NEAR(N, M, EPSILON_MM)
       << "At pole, N should equal M";
}

TEST(ReferenceTest, RadiiOfCurvature_AllLatitudes)
{
   auto radii = get_radii_of_curvature();
   
   for (const auto& r : radii)
   {
      auto coord = ECEFCoordinate::from_geodetic_deg(r.lat_deg, 0.0, 0.0);
      
      double M = coord.meridian_radius();
      double N = coord.prime_vertical_radius();
      
      EXPECT_NEAR(M, r.meridian_radius_m, EPSILON_M)
          << "Meridian radius mismatch at " << r.lat_deg << "°";
      EXPECT_NEAR(N, r.prime_vertical_radius_m, EPSILON_M)
          << "Prime vertical radius mismatch at " << r.lat_deg << "°";
   }
}

TEST(ReferenceTest, RadiiOfCurvature_Monotonicity)
{
   // M increases from equator to pole
   // N increases from equator to pole
   
   auto M_eq = ECEFCoordinate::from_geodetic_deg(0.0, 0.0, 0.0).meridian_radius();
   auto M_30 = ECEFCoordinate::from_geodetic_deg(30.0, 0.0, 0.0).meridian_radius();
   auto M_60 = ECEFCoordinate::from_geodetic_deg(60.0, 0.0, 0.0).meridian_radius();
   auto M_90 = ECEFCoordinate::from_geodetic_deg(90.0, 0.0, 0.0).meridian_radius();
   
   EXPECT_LT(M_eq, M_30);
   EXPECT_LT(M_30, M_60);
   EXPECT_LT(M_60, M_90);
   
   auto N_eq = ECEFCoordinate::from_geodetic_deg(0.0, 0.0, 0.0).prime_vertical_radius();
   auto N_30 = ECEFCoordinate::from_geodetic_deg(30.0, 0.0, 0.0).prime_vertical_radius();
   auto N_60 = ECEFCoordinate::from_geodetic_deg(60.0, 0.0, 0.0).prime_vertical_radius();
   auto N_90 = ECEFCoordinate::from_geodetic_deg(90.0, 0.0, 0.0).prime_vertical_radius();
   
   EXPECT_LT(N_eq, N_30);
   EXPECT_LT(N_30, N_60);
   EXPECT_LT(N_60, N_90);
}

// ==================== Airport Coordinates Validation ====================

TEST(ReferenceTest, Airports_GeodeticRoundTrip)
{
   auto airports = get_airports();
   
   for (const auto& airport : airports)
   {
      auto coord = ECEFCoordinate::from_geodetic_deg(
          airport.lat_deg, airport.lon_deg, airport.elevation_m);
      
      // Convert back
      double lat_back = coord.latitude_deg();
      double lon_back = coord.longitude_deg();
      double alt_back = coord.altitude();
      
      EXPECT_NEAR(lat_back, airport.lat_deg, EPSILON_DEG)
          << "Latitude round-trip failed for " << airport.name;
      EXPECT_NEAR(lon_back, airport.lon_deg, EPSILON_DEG)
          << "Longitude round-trip failed for " << airport.name;
      EXPECT_NEAR(alt_back, airport.elevation_m, EPSILON_CM)
          << "Altitude round-trip failed for " << airport.name;
   }
}

TEST(ReferenceTest, Airports_ECEFMagnitude)
{
   auto airports = get_airports();
   
   for (const auto& airport : airports)
   {
      auto coord = ECEFCoordinate::from_geodetic_deg(
          airport.lat_deg, airport.lon_deg, airport.elevation_m);
      
      // ECEF magnitude should be approximately Earth radius + altitude
      double ecef_mag = coord.ecef().norm();
      double expected_mag = WGS84::a + airport.elevation_m;
      
      // Allow 1% tolerance due to ellipsoid shape
      EXPECT_NEAR(ecef_mag, expected_mag, expected_mag * 0.01)
          << "ECEF magnitude unexpected for " << airport.name;
   }
}

// ==================== High-Precision Conversions ====================

TEST(ReferenceTest, HighPrecision_45N_0E)
{
   // Test point from GeographicLib: 45°N, 0°E, 0m
   // ECEF: (4517590.878, 0.0, 4487348.409)
   
   auto coord = ECEFCoordinate::from_geodetic_deg(45.0, 0.0, 0.0);
   
   EXPECT_NEAR(coord.x(), 4517590.878, EPSILON_M);
   EXPECT_NEAR(coord.y(), 0.0, EPSILON_M);
   EXPECT_NEAR(coord.z(), 4487348.409, EPSILON_M);
   
   // Round-trip
   auto coord2 = ECEFCoordinate(4517590.878, 0.0, 4487348.409);
   EXPECT_NEAR(coord2.latitude_deg(), 45.0, EPSILON_DEG);
   EXPECT_NEAR(coord2.longitude_deg(), 0.0, EPSILON_DEG);
   EXPECT_NEAR(coord2.altitude(), 0.0, EPSILON_CM);
}

TEST(ReferenceTest, HighPrecision_45S_0E)
{
   // Test point: 45°S, 0°E, 0m
   // ECEF: (4517590.878, 0.0, -4487348.409)
   
   auto coord = ECEFCoordinate::from_geodetic_deg(-45.0, 0.0, 0.0);
   
   EXPECT_NEAR(coord.x(), 4517590.878, EPSILON_M);
   EXPECT_NEAR(coord.y(), 0.0, EPSILON_M);
   EXPECT_NEAR(coord.z(), -4487348.409, EPSILON_M);
   
   // Round-trip
   auto coord2 = ECEFCoordinate(4517590.878, 0.0, -4487348.409);
   EXPECT_NEAR(coord2.latitude_deg(), -45.0, EPSILON_DEG);
   EXPECT_NEAR(coord2.longitude_deg(), 0.0, EPSILON_DEG);
   EXPECT_NEAR(coord2.altitude(), 0.0, EPSILON_CM);
}

// ==================== Altitude Effects ====================

TEST(ReferenceTest, Altitude_ECEF_Scaling)
{
   // Test that altitude increases ECEF magnitude proportionally
   
   auto coord_0 = ECEFCoordinate::from_geodetic_deg(45.0, 0.0, 0.0);
   auto coord_1000 = ECEFCoordinate::from_geodetic_deg(45.0, 0.0, 1000.0);
   auto coord_10000 = ECEFCoordinate::from_geodetic_deg(45.0, 0.0, 10000.0);
   
   double mag_0 = coord_0.ecef().norm();
   double mag_1000 = coord_1000.ecef().norm();
   double mag_10000 = coord_10000.ecef().norm();
   
   // Magnitude should increase by approximately the altitude
   EXPECT_NEAR(mag_1000 - mag_0, 1000.0, EPSILON_M);
   EXPECT_NEAR(mag_10000 - mag_0, 10000.0, EPSILON_M);
}

TEST(ReferenceTest, Altitude_Direction_Preserved)
{
   // Test that altitude changes preserve horizontal direction
   
   auto coord_0 = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
   auto coord_1000 = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 1000.0);
   
   // Directions should be nearly parallel
   auto dir_0 = coord_0.ecef().normalized();
   auto dir_1000 = coord_1000.ecef().normalized();
   
   double dot = dir_0.dot(dir_1000);
   EXPECT_NEAR(dot, 1.0, 1e-6)
       << "Altitude change should preserve direction";
}

// ==================== Symmetry Tests ====================

TEST(ReferenceTest, Symmetry_NorthSouth)
{
   // Test that North and South hemispheres are symmetric
   
   auto north = ECEFCoordinate::from_geodetic_deg(45.0, 30.0, 100.0);
   auto south = ECEFCoordinate::from_geodetic_deg(-45.0, 30.0, 100.0);
   
   // X and Y should be equal, Z should be opposite
   EXPECT_NEAR(north.x(), south.x(), EPSILON_MM);
   EXPECT_NEAR(north.y(), south.y(), EPSILON_MM);
   EXPECT_NEAR(north.z(), -south.z(), EPSILON_MM);
}

TEST(ReferenceTest, Symmetry_EastWest)
{
   // Test that East and West hemispheres are symmetric
   
   auto east = ECEFCoordinate::from_geodetic_deg(45.0, 30.0, 100.0);
   auto west = ECEFCoordinate::from_geodetic_deg(45.0, -30.0, 100.0);
   
   // X and Z should be equal, Y should be opposite
   EXPECT_NEAR(east.x(), west.x(), EPSILON_MM);
   EXPECT_NEAR(east.y(), -west.y(), EPSILON_MM);
   EXPECT_NEAR(east.z(), west.z(), EPSILON_MM);
}

TEST(ReferenceTest, Symmetry_Antipodal)
{
   // Test antipodal points (opposite sides of Earth)
   
   auto point1 = ECEFCoordinate::from_geodetic_deg(30.0, 45.0, 0.0);
   auto point2 = ECEFCoordinate::from_geodetic_deg(-30.0, -135.0, 0.0);
   
   // Antipodal points should have opposite ECEF vectors
   auto ecef1 = point1.ecef();
   auto ecef2 = point2.ecef();
   
   EXPECT_NEAR(ecef1.x, -ecef2.x, EPSILON_M);
   EXPECT_NEAR(ecef1.y, -ecef2.y, EPSILON_M);
   EXPECT_NEAR(ecef1.z, -ecef2.z, EPSILON_M);
}

// ==================== Consistency Tests ====================

TEST(ReferenceTest, Consistency_MultipleConstructions)
{
   // Test that different construction methods give same result
   
   double lat_rad = 37.7749 * DEG_TO_RAD;
   double lon_rad = -122.4194 * DEG_TO_RAD;
   double alt = 10.0;
   
   auto coord1 = ECEFCoordinate::from_geodetic(lat_rad, lon_rad, alt);
   auto coord2 = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 10.0);
   
   EXPECT_NEAR(coord1.x(), coord2.x(), EPSILON_MM);
   EXPECT_NEAR(coord1.y(), coord2.y(), EPSILON_MM);
   EXPECT_NEAR(coord1.z(), coord2.z(), EPSILON_MM);
}

TEST(ReferenceTest, Consistency_CachedValues)
{
   // Test that cached trig values are consistent
   
   auto coord = ECEFCoordinate::from_geodetic_deg(45.0, 30.0, 100.0);
   
   double lat = coord.latitude();
   double lon = coord.longitude();
   
   // Cached trig should match computed trig
   EXPECT_NEAR(coord.sin_lat(), std::sin(lat), EPSILON_MICRO);
   EXPECT_NEAR(coord.cos_lat(), std::cos(lat), EPSILON_MICRO);
   EXPECT_NEAR(coord.sin_lon(), std::sin(lon), EPSILON_MICRO);
   EXPECT_NEAR(coord.cos_lon(), std::cos(lon), EPSILON_MICRO);
}

TEST(ReferenceTest, Consistency_TrigIdentities)
{
   // Test that trig values satisfy fundamental identities
   
   auto coord = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 100.0);
   
   double sl = coord.sin_lat();
   double cl = coord.cos_lat();
   double sn = coord.sin_lon();
   double cn = coord.cos_lon();
   
   // sin²θ + cos²θ = 1
   EXPECT_NEAR(sl*sl + cl*cl, 1.0, EPSILON_MICRO);
   EXPECT_NEAR(sn*sn + cn*cn, 1.0, EPSILON_MICRO);
}
