#include <dstl/Geodesy.h>
#include <gtest/gtest.h>
#include <cmath>
#include <unordered_map>

using namespace dstl::geo;
using namespace dstl::linalg;

constexpr double EPSILON = 1e-6;
constexpr double EPSILON_METERS = 0.01; // 1cm tolerance

// Helper to compare doubles
bool approx_equal(double a, double b, double eps = EPSILON)
{
   return std::abs(a - b) < eps;
}

// ==================== Construction Tests ====================

TEST(GeodesyTest, ECEFConstruction)
{
   ECEFCoordinate coord(1000000, 2000000, 3000000);
   EXPECT_DOUBLE_EQ(coord.x(), 1000000);
   EXPECT_DOUBLE_EQ(coord.y(), 2000000);
   EXPECT_DOUBLE_EQ(coord.z(), 3000000);
}

TEST(GeodesyTest, GeodeticConstruction)
{
   // Test point: 45°N, 120°W, 100m altitude
   double lat = 45.0 * M_PI / 180.0;
   double lon = -120.0 * M_PI / 180.0;
   double alt = 100.0;

   auto coord = ECEFCoordinate::from_geodetic(lat, lon, alt);

   // Verify we can get ECEF coordinates
   EXPECT_GT(std::abs(coord.x()), 1e6); // Should be millions of meters
   EXPECT_GT(std::abs(coord.y()), 1e6);
   EXPECT_GT(std::abs(coord.z()), 1e6);
}

TEST(GeodesyTest, GeodeticConstructionDegrees)
{
   auto coord = ECEFCoordinate::from_geodetic_deg(45.0, -120.0, 100.0);

   EXPECT_GT(std::abs(coord.x()), 1e6);
   EXPECT_GT(std::abs(coord.y()), 1e6);
   EXPECT_GT(std::abs(coord.z()), 1e6);
}

// ==================== Round-Trip Tests ====================

TEST(GeodesyTest, GeodeticRoundTrip)
{
   // Test various points around the globe
   std::vector<std::tuple<double, double, double>> test_points = {
       {0.0, 0.0, 0.0},                                  // Equator, Prime Meridian
       {45.0, -120.0, 100.0},                            // Mid-latitude
       {-33.8688, 151.2093, 50.0},                       // Sydney
       {37.7749, -122.4194, 10.0},                       // San Francisco
       {89.0, 0.0, 0.0},                                 // Near North Pole
       {-89.0, 0.0, 0.0},                                // Near South Pole
   };

   for (const auto &[lat_deg, lon_deg, alt] : test_points)
   {
      auto coord = ECEFCoordinate::from_geodetic_deg(lat_deg, lon_deg, alt);

      // Convert back to geodetic
      double lat_back = coord.latitude_deg();
      double lon_back = coord.longitude_deg();
      double alt_back = coord.altitude();

      EXPECT_TRUE(approx_equal(lat_deg, lat_back, 1e-6))
          << "Latitude mismatch for (" << lat_deg << ", " << lon_deg << ", " << alt << ")";
      EXPECT_TRUE(approx_equal(lon_deg, lon_back, 1e-6))
          << "Longitude mismatch for (" << lat_deg << ", " << lon_deg << ", " << alt << ")";
      EXPECT_TRUE(approx_equal(alt, alt_back, EPSILON_METERS))
          << "Altitude mismatch for (" << lat_deg << ", " << lon_deg << ", " << alt << ")";
   }
}

// ==================== Caching Tests ====================

TEST(GeodesyTest, GeodeticCaching)
{
   auto coord = ECEFCoordinate::from_geodetic_deg(45.0, -120.0, 100.0);

   // First access should compute
   EXPECT_TRUE(coord.has_geodetic_cache());

   // Subsequent accesses should use cache (we can't directly test this,
   // but we verify consistency)
   double lat1 = coord.latitude();
   double lat2 = coord.latitude();
   EXPECT_DOUBLE_EQ(lat1, lat2);
}

TEST(GeodesyTest, CacheInvalidation)
{
   auto coord = ECEFCoordinate::from_geodetic_deg(45.0, -120.0, 100.0);

   double lat_before = coord.latitude();

   // Modify ECEF coordinates
   coord += Vec3{1000, 0, 0};

   // Geodetic should change
   double lat_after = coord.latitude();
   EXPECT_NE(lat_before, lat_after);
}

TEST(GeodesyTest, TrigCaching)
{
   auto coord = ECEFCoordinate::from_geodetic_deg(45.0, -120.0, 100.0);

   // Access trig values
   double sl = coord.sin_lat();
   double cl = coord.cos_lat();
   double sn = coord.sin_lon();
   double cn = coord.cos_lon();

   // Verify they're consistent with lat/lon
   double lat = coord.latitude();
   double lon = coord.longitude();

   EXPECT_TRUE(approx_equal(sl, std::sin(lat), 1e-10));
   EXPECT_TRUE(approx_equal(cl, std::cos(lat), 1e-10));
   EXPECT_TRUE(approx_equal(sn, std::sin(lon), 1e-10));
   EXPECT_TRUE(approx_equal(cn, std::cos(lon), 1e-10));
}

// ==================== Vector Operations Tests ====================

TEST(GeodesyTest, VectorTo)
{
   auto coord1 = ECEFCoordinate::from_geodetic_deg(0.0, 0.0, 0.0);
   auto coord2 = ECEFCoordinate::from_geodetic_deg(0.0, 1.0, 0.0);

   Vec3 vec = coord1.vector_to(coord2);

   // Should be non-zero
   EXPECT_GT(vec.norm(), 1000.0); // At least 1km
}

TEST(GeodesyTest, DistanceTo)
{
   auto coord1 = ECEFCoordinate::from_geodetic_deg(0.0, 0.0, 0.0);
   auto coord2 = ECEFCoordinate::from_geodetic_deg(0.0, 1.0, 0.0);

   double dist = coord1.distance_to(coord2);

   // 1 degree at equator ≈ 111 km
   EXPECT_GT(dist, 100000.0);  // > 100 km
   EXPECT_LT(dist, 120000.0);  // < 120 km
}

TEST(GeodesyTest, ECEFOffset)
{
   auto coord1 = ECEFCoordinate::from_geodetic_deg(45.0, -120.0, 100.0);
   Vec3 offset{1000, 2000, 3000};

   auto coord2 = coord1 + offset;

   // Verify offset was applied
   EXPECT_DOUBLE_EQ(coord2.x(), coord1.x() + 1000);
   EXPECT_DOUBLE_EQ(coord2.y(), coord1.y() + 2000);
   EXPECT_DOUBLE_EQ(coord2.z(), coord1.z() + 3000);
}

TEST(GeodesyTest, ECEFOffsetInPlace)
{
   auto coord = ECEFCoordinate::from_geodetic_deg(45.0, -120.0, 100.0);
   double x_orig = coord.x();

   coord += Vec3{1000, 0, 0};

   EXPECT_DOUBLE_EQ(coord.x(), x_orig + 1000);
}

// ==================== ENU Frame Tests ====================

TEST(GeodesyTest, ENUFrameCreation)
{
   auto coord = ECEFCoordinate::from_geodetic_deg(45.0, -120.0, 100.0);
   auto frame = coord.local_frame();

   // Verify frame vectors are unit vectors
   EXPECT_TRUE(approx_equal(frame.east.norm(), 1.0, 1e-10));
   EXPECT_TRUE(approx_equal(frame.north.norm(), 1.0, 1e-10));
   EXPECT_TRUE(approx_equal(frame.up.norm(), 1.0, 1e-10));

   // Verify frame vectors are orthogonal
   EXPECT_TRUE(approx_equal(frame.east.dot(frame.north), 0.0, 1e-10));
   EXPECT_TRUE(approx_equal(frame.east.dot(frame.up), 0.0, 1e-10));
   EXPECT_TRUE(approx_equal(frame.north.dot(frame.up), 0.0, 1e-10));
}

TEST(GeodesyTest, ENURoundTrip)
{
   auto coord = ECEFCoordinate::from_geodetic_deg(45.0, -120.0, 100.0);
   auto frame = coord.local_frame();

   Vec3 enu_offset{100, 200, 50}; // 100m E, 200m N, 50m U

   // Convert to ECEF and back
   Vec3 ecef_offset = frame.to_ecef(enu_offset);
   Vec3 enu_back = frame.to_enu(ecef_offset);

   EXPECT_TRUE(approx_equal(enu_offset.x, enu_back.x, 1e-6));
   EXPECT_TRUE(approx_equal(enu_offset.y, enu_back.y, 1e-6));
   EXPECT_TRUE(approx_equal(enu_offset.z, enu_back.z, 1e-6));
}

TEST(GeodesyTest, ENUCoordinateConversion)
{
   auto origin = ECEFCoordinate::from_geodetic_deg(45.0, -120.0, 100.0);
   auto frame = origin.local_frame();

   // Create a point 1km north
   Vec3 enu_offset{0, 1000, 0};
   auto target = frame.enu_to_coord(enu_offset);

   // Convert back to ENU
   Vec3 enu_back = frame.coord_to_enu(target);

   EXPECT_TRUE(approx_equal(enu_offset.x, enu_back.x, EPSILON_METERS));
   EXPECT_TRUE(approx_equal(enu_offset.y, enu_back.y, EPSILON_METERS));
   EXPECT_TRUE(approx_equal(enu_offset.z, enu_back.z, EPSILON_METERS));
}

// ==================== Special Cases ====================

TEST(GeodesyTest, EquatorPrimeMeridian)
{
   auto coord = ECEFCoordinate::from_geodetic_deg(0.0, 0.0, 0.0);

   EXPECT_TRUE(approx_equal(coord.latitude_deg(), 0.0, 1e-6));
   EXPECT_TRUE(approx_equal(coord.longitude_deg(), 0.0, 1e-6));
   EXPECT_TRUE(approx_equal(coord.altitude(), 0.0, EPSILON_METERS));

   // At equator on prime meridian, X should be ≈ Earth radius
   EXPECT_TRUE(approx_equal(coord.x(), WGS84::a, 1.0));
   EXPECT_TRUE(approx_equal(coord.y(), 0.0, 1.0));
   EXPECT_TRUE(approx_equal(coord.z(), 0.0, 1.0));
}

TEST(GeodesyTest, NorthPole)
{
   auto coord = ECEFCoordinate::from_geodetic_deg(90.0, 0.0, 0.0);

   EXPECT_TRUE(approx_equal(coord.latitude_deg(), 90.0, 1e-6));
   EXPECT_TRUE(approx_equal(coord.altitude(), 0.0, EPSILON_METERS));

   // At pole, X and Y should be near zero, Z should be ≈ polar radius
   EXPECT_TRUE(approx_equal(coord.x(), 0.0, 1.0));
   EXPECT_TRUE(approx_equal(coord.y(), 0.0, 1.0));
   EXPECT_TRUE(approx_equal(coord.z(), WGS84::b, 1.0));
}

TEST(GeodesyTest, SouthPole)
{
   auto coord = ECEFCoordinate::from_geodetic_deg(-90.0, 0.0, 0.0);

   EXPECT_TRUE(approx_equal(coord.latitude_deg(), -90.0, 1e-6));
   EXPECT_TRUE(approx_equal(coord.altitude(), 0.0, EPSILON_METERS));

   EXPECT_TRUE(approx_equal(coord.x(), 0.0, 1.0));
   EXPECT_TRUE(approx_equal(coord.y(), 0.0, 1.0));
   EXPECT_TRUE(approx_equal(coord.z(), -WGS84::b, 1.0));
}

// ==================== Performance Characteristics ====================

TEST(GeodesyTest, NoTrigForECEFOperations)
{
   // This test verifies that ECEF operations don't trigger geodetic computation
   auto coord1 = ECEFCoordinate(1000000, 2000000, 3000000);
   auto coord2 = ECEFCoordinate(1001000, 2001000, 3001000);

   // These operations should not compute geodetic
   EXPECT_FALSE(coord1.has_geodetic_cache());
   EXPECT_FALSE(coord2.has_geodetic_cache());

   double dist = coord1.distance_to(coord2);
   EXPECT_GT(dist, 0);

   // Still no geodetic computation
   EXPECT_FALSE(coord1.has_geodetic_cache());
   EXPECT_FALSE(coord2.has_geodetic_cache());

   Vec3 vec = coord1.vector_to(coord2);
   EXPECT_GT(vec.norm(), 0);

   // Still no geodetic computation
   EXPECT_FALSE(coord1.has_geodetic_cache());
   EXPECT_FALSE(coord2.has_geodetic_cache());
}

TEST(GeodesyTest, PrecomputeGeodetic)
{
   auto coord = ECEFCoordinate(1000000, 2000000, 3000000);

   EXPECT_FALSE(coord.has_geodetic_cache());

   coord.precompute_geodetic();

   EXPECT_TRUE(coord.has_geodetic_cache());
}

// ==================== Geodesic Distance Tests ====================

TEST(GeodesyTest, GeodesicDistance)
{
   // Test SF to LA (known distance ~559 km)
   auto sf = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 0);
   auto la = ECEFCoordinate::from_geodetic_deg(34.0522, -118.2437, 0);

   double geodesic_dist = sf.geodesic_distance_to(la);
   double chord_dist = sf.distance_to(la);

   // Geodesic should be longer than chord
   EXPECT_GT(geodesic_dist, chord_dist);

   // Should be approximately 559 km
   EXPECT_GT(geodesic_dist, 550000.0);
   EXPECT_LT(geodesic_dist, 570000.0);
}

TEST(GeodesyTest, GeodesicDistanceShort)
{
   // For short distances, geodesic and chord should be similar
   auto coord1 = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 0);
   auto coord2 = ECEFCoordinate::from_geodetic_deg(37.7750, -122.4194, 0);

   double geodesic_dist = coord1.geodesic_distance_to(coord2);
   double chord_dist = coord1.distance_to(coord2);

   // Should differ by less than 0.1%
   EXPECT_TRUE(approx_equal(geodesic_dist, chord_dist, chord_dist * 0.001));
}

TEST(GeodesyTest, GeodesicDistanceAntipodal)
{
   // Test near-antipodal points (opposite sides of Earth)
   auto coord1 = ECEFCoordinate::from_geodetic_deg(0, 0, 0);
   auto coord2 = ECEFCoordinate::from_geodetic_deg(0, 179, 0);

   double geodesic_dist = coord1.geodesic_distance_to(coord2);

   // 179° along equator (not 180°)
   double expected_dist = (179.0 / 180.0) * PI * WGS84::a;
   EXPECT_TRUE(approx_equal(geodesic_dist, expected_dist, 1000.0));
}

// ==================== Bearing Tests ====================

TEST(GeodesyTest, BearingNorth)
{
   auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0);
   auto north = ECEFCoordinate::from_geodetic_deg(38.0, -122.0, 0);

   double bearing = origin.bearing_to(north);

   // Should be approximately 0 (North) - handle 360° wrap-around
   // Normalize bearing to [0, 2π) and check if close to 0 or 2π
   double bearing_norm = std::fmod(bearing + 2.0 * PI, 2.0 * PI);
   bool is_north = (bearing_norm < 0.01) || (bearing_norm > 2.0 * PI - 0.01);
   EXPECT_TRUE(is_north) << "Bearing was " << bearing << " rad (" << (bearing * RAD_TO_DEG) << "°)";
}

TEST(GeodesyTest, BearingEast)
{
   auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0);
   auto east = ECEFCoordinate::from_geodetic_deg(37.0, -121.0, 0);

   double bearing = origin.bearing_to(east);

   // Should be approximately π/2 (East)
   EXPECT_TRUE(approx_equal(bearing, PI / 2.0, 0.01));
}

TEST(GeodesyTest, BearingSouth)
{
   auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0);
   auto south = ECEFCoordinate::from_geodetic_deg(36.0, -122.0, 0);

   double bearing = origin.bearing_to(south);

   // Should be approximately π (South)
   EXPECT_TRUE(approx_equal(bearing, PI, 0.01));
}

TEST(GeodesyTest, BearingWest)
{
   auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0);
   auto west = ECEFCoordinate::from_geodetic_deg(37.0, -123.0, 0);

   double bearing = origin.bearing_to(west);

   // Should be approximately 3π/2 (West)
   EXPECT_TRUE(approx_equal(bearing, 3.0 * PI / 2.0, 0.01));
}

// ==================== Movement Tests ====================

TEST(GeodesyTest, MoveByBearing)
{
   auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0);

   // Move 1km North
   auto north = origin.move_by_bearing(1000.0, 0.0);

   // Latitude should increase
   EXPECT_GT(north.latitude_deg(), origin.latitude_deg());

   // Longitude should be approximately the same
   EXPECT_TRUE(approx_equal(north.longitude_deg(), origin.longitude_deg(), 0.001));
}

TEST(GeodesyTest, MoveByBearingRoundTrip)
{
   auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0);

   // Move 10km at 45° bearing
   double distance = 10000.0;
   double bearing = PI / 4.0;  // 45° (NE)

   auto target = origin.move_by_bearing(distance, bearing);

   // Compute distance and bearing back
   double dist_back = origin.geodesic_distance_to(target);
   double bearing_back = origin.bearing_to(target);

   // Distance should match (within 1%)
   EXPECT_TRUE(approx_equal(distance, dist_back, distance * 0.01));

   // Bearing should match (within 1°)
   EXPECT_TRUE(approx_equal(bearing, bearing_back, 0.02));
}

// ==================== NED Frame Tests ====================

TEST(GeodesyTest, NEDFrameCreation)
{
   auto coord = ECEFCoordinate::from_geodetic_deg(45.0, -120.0, 100.0);
   auto frame = coord.local_ned_frame();

   // Verify frame vectors are unit vectors
   EXPECT_TRUE(approx_equal(frame.north.norm(), 1.0, 1e-10));
   EXPECT_TRUE(approx_equal(frame.east.norm(), 1.0, 1e-10));
   EXPECT_TRUE(approx_equal(frame.down.norm(), 1.0, 1e-10));

   // Verify frame vectors are orthogonal
   EXPECT_TRUE(approx_equal(frame.north.dot(frame.east), 0.0, 1e-10));
   EXPECT_TRUE(approx_equal(frame.north.dot(frame.down), 0.0, 1e-10));
   EXPECT_TRUE(approx_equal(frame.east.dot(frame.down), 0.0, 1e-10));
}

TEST(GeodesyTest, NEDRoundTrip)
{
   auto coord = ECEFCoordinate::from_geodetic_deg(45.0, -120.0, 100.0);
   auto frame = coord.local_ned_frame();

   Vec3 ned_offset{100, 200, -50}; // 100m N, 200m E, 50m Up

   // Convert to ECEF and back
   Vec3 ecef_offset = frame.to_ecef(ned_offset);
   Vec3 ned_back = frame.to_ned(ecef_offset);

   EXPECT_TRUE(approx_equal(ned_offset.x, ned_back.x, 1e-6));
   EXPECT_TRUE(approx_equal(ned_offset.y, ned_back.y, 1e-6));
   EXPECT_TRUE(approx_equal(ned_offset.z, ned_back.z, 1e-6));
}

// ==================== String Conversion Tests ====================

TEST(GeodesyTest, ToString)
{
   auto coord = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 10.0);

   std::string str = coord.to_string();

   // Should contain key information
   EXPECT_NE(str.find("37.7"), std::string::npos);
   EXPECT_NE(str.find("122.4"), std::string::npos);
   EXPECT_NE(str.find("10"), std::string::npos);
}

TEST(GeodesyTest, ToStringECEF)
{
   auto coord = ECEFCoordinate(1000000, 2000000, 3000000);

   std::string str = coord.to_string_ecef();

   // Should contain ECEF coordinates
   EXPECT_NE(str.find("1000000"), std::string::npos);
   EXPECT_NE(str.find("2000000"), std::string::npos);
   EXPECT_NE(str.find("3000000"), std::string::npos);
}

// ==================== Equality Tests ====================

TEST(GeodesyTest, EqualityOperator)
{
   auto coord1 = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 10.0);
   auto coord2 = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 10.0);
   auto coord3 = ECEFCoordinate::from_geodetic_deg(37.7750, -122.4194, 10.0);

   EXPECT_TRUE(coord1 == coord2);
   EXPECT_FALSE(coord1 == coord3);
   EXPECT_TRUE(coord1 != coord3);
}

TEST(GeodesyTest, ApproxEqual)
{
   auto coord1 = ECEFCoordinate(1000000.0, 2000000.0, 3000000.0);
   auto coord2 = ECEFCoordinate(1000000.5, 2000000.5, 3000000.5);
   auto coord3 = ECEFCoordinate(1000010.0, 2000000.0, 3000000.0);

   // Within 1m tolerance
   EXPECT_TRUE(coord1.approx_equal(coord2, 1.0));

   // Not within 1m tolerance (exactly 10m apart)
   EXPECT_FALSE(coord1.approx_equal(coord3, 1.0));

   // Within 10m tolerance (need slightly more than 10m since distance is exactly 10m)
   EXPECT_TRUE(coord1.approx_equal(coord3, 10.1));
}

// ==================== Hash Tests ====================

TEST(GeodesyTest, HashFunction)
{
   auto coord1 = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 10.0);
   auto coord2 = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 10.0);
   auto coord3 = ECEFCoordinate::from_geodetic_deg(37.7750, -122.4194, 10.0);

   std::hash<ECEFCoordinate> hasher;

   // Same coordinates should have same hash
   EXPECT_EQ(hasher(coord1), hasher(coord2));

   // Different coordinates should (probably) have different hash
   EXPECT_NE(hasher(coord1), hasher(coord3));
}

TEST(GeodesyTest, UnorderedMapUsage)
{
   // Test that ECEFCoordinate can be used in unordered_map
   std::unordered_map<ECEFCoordinate, std::string> coord_map;

   auto sf = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 10.0);
   auto la = ECEFCoordinate::from_geodetic_deg(34.0522, -118.2437, 50.0);

   coord_map[sf] = "San Francisco";
   coord_map[la] = "Los Angeles";

   EXPECT_EQ(coord_map[sf], "San Francisco");
   EXPECT_EQ(coord_map[la], "Los Angeles");
   EXPECT_EQ(coord_map.size(), 2);
}

// ==================== Flat-Earth Displacement Tests ====================

TEST(GeodesyTest, ApplyENUDisplacement)
{
   auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 1000.0);

   // Move 10km East, 5km North, climb 500m
   Vec3 enu_disp{10000.0, 5000.0, 500.0};
   auto target = origin.apply_enu_displacement(enu_disp);

   // Verify displacement was applied
   EXPECT_GT(target.longitude_deg(), origin.longitude_deg());  // Moved East
   EXPECT_GT(target.latitude_deg(), origin.latitude_deg());    // Moved North
   EXPECT_GT(target.altitude(), origin.altitude());            // Climbed

   // Verify approximate distance
   double dist = origin.geodesic_distance_to(target);
   double expected_dist = std::sqrt(10000.0 * 10000.0 + 5000.0 * 5000.0);
   EXPECT_TRUE(approx_equal(dist, expected_dist, 100.0));  // Within 100m
}

TEST(GeodesyTest, ApplyNEDDisplacement)
{
   auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 1000.0);

   // Move 10km North, 5km East, descend 500m (positive down)
   Vec3 ned_disp{10000.0, 5000.0, 500.0};
   auto target = origin.apply_ned_displacement(ned_disp);

   // Verify displacement was applied
   EXPECT_GT(target.latitude_deg(), origin.latitude_deg());   // Moved North
   EXPECT_GT(target.longitude_deg(), origin.longitude_deg()); // Moved East
   EXPECT_LT(target.altitude(), origin.altitude());           // Descended

   // Verify altitude change (flat-earth approximation has ~10m error over 10km)
   double alt_change = target.altitude() - origin.altitude();
   EXPECT_TRUE(approx_equal(alt_change, -500.0, 15.0));  // Within 15m
}

TEST(GeodesyTest, FlatEarthRoundTrip)
{
   auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 1000.0);
   auto frame = origin.local_frame();

   // Apply displacement and convert back to ENU
   Vec3 enu_disp{1000.0, 2000.0, 100.0};
   auto target = origin.apply_enu_displacement(enu_disp);
   Vec3 enu_back = frame.coord_to_enu(target);

   // Should match original displacement (within tolerance)
   EXPECT_TRUE(approx_equal(enu_disp.x, enu_back.x, 0.1));
   EXPECT_TRUE(approx_equal(enu_disp.y, enu_back.y, 0.1));
   EXPECT_TRUE(approx_equal(enu_disp.z, enu_back.z, 0.1));
}

TEST(GeodesyTest, ENUTrajectory)
{
   auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 1000.0);

   // Multi-leg trajectory
   std::vector<Vec3> legs = {
       {10000.0, 10000.0, 500.0},   // 10km E, 10km N, climb 500m
       {5000.0, -5000.0, -200.0},   // 5km E, 5km S, descend 200m
       {-3000.0, 8000.0, 100.0}     // 3km W, 8km N, climb 100m
   };

   auto final_pos = origin.apply_enu_trajectory(legs);

   // Verify we moved from origin
   EXPECT_NE(final_pos.latitude_deg(), origin.latitude_deg());
   EXPECT_NE(final_pos.longitude_deg(), origin.longitude_deg());

   // Verify altitude change (flat-earth approximation accumulates error)
   double total_alt_change = 500.0 - 200.0 + 100.0;
   double actual_alt_change = final_pos.altitude() - origin.altitude();
   EXPECT_TRUE(approx_equal(actual_alt_change, total_alt_change, 30.0));  // Within 30m
}

TEST(GeodesyTest, NEDTrajectory)
{
   auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 1000.0);

   // Multi-leg trajectory in NED
   std::vector<Vec3> legs = {
       {10000.0, 5000.0, -500.0},   // 10km N, 5km E, climb 500m
       {-5000.0, 10000.0, 200.0},   // 5km S, 10km E, descend 200m
   };

   auto final_pos = origin.apply_ned_trajectory(legs);

   // Verify we moved
   double dist = origin.geodesic_distance_to(final_pos);
   EXPECT_GT(dist, 10000.0);  // At least 10km

   // Verify altitude change (flat-earth approximation has error)
   double total_alt_change = 500.0 - 200.0;
   double actual_alt_change = final_pos.altitude() - origin.altitude();
   EXPECT_TRUE(approx_equal(actual_alt_change, total_alt_change, 25.0));  // Within 25m
}

TEST(GeodesyTest, FlatEarthAccuracy)
{
   auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 1000.0);

   // Test accuracy at different distances
   std::vector<double> distances = {1000.0, 10000.0, 50000.0, 100000.0};

   for (double dist : distances)
   {
      // Apply flat-earth displacement (due East)
      Vec3 enu_disp{dist, 0.0, 0.0};
      auto target = origin.apply_enu_displacement(enu_disp);

      // Measure actual geodesic distance
      double geodesic_dist = origin.geodesic_distance_to(target);
      double error = std::abs(geodesic_dist - dist);

      // Error should be small for short distances
      if (dist <= 10000.0)
      {
         EXPECT_LT(error, 10.0);  // < 10m error for 10km
      }
      else if (dist <= 50000.0)
      {
         EXPECT_LT(error, 100.0);  // < 100m error for 50km
      }
   }
}

// ==================== Geometric Utilities Tests ====================

TEST(GeodesyTest, SurfaceNormal)
{
   auto coord = ECEFCoordinate::from_geodetic_deg(45.0, -120.0, 0.0);
   auto normal = coord.surface_normal();

   // Normal should be unit vector
   EXPECT_TRUE(approx_equal(normal.norm(), 1.0, 1e-10));

   // Normal should point away from Earth center
   // For a point on the surface, normal should be roughly aligned with ECEF position
   auto ecef_normalized = coord.ecef().normalized();
   double alignment = normal.dot(ecef_normalized);
   EXPECT_GT(alignment, 0.99);  // Should be nearly parallel
}

TEST(GeodesyTest, Midpoint)
{
   auto sf = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 10.0);
   auto la = ECEFCoordinate::from_geodetic_deg(34.0522, -118.2437, 50.0);

   auto mid = sf.midpoint(la);

   // Midpoint should be roughly equidistant from both
   double dist_sf = sf.geodesic_distance_to(mid);
   double dist_la = la.geodesic_distance_to(mid);
   double total_dist = sf.geodesic_distance_to(la);

   EXPECT_TRUE(approx_equal(dist_sf, total_dist / 2.0, 1000.0));  // Within 1km
   EXPECT_TRUE(approx_equal(dist_la, total_dist / 2.0, 1000.0));

   // Altitude should be average
   double expected_alt = (sf.altitude() + la.altitude()) / 2.0;
   EXPECT_TRUE(approx_equal(mid.altitude(), expected_alt, 10.0));
}

TEST(GeodesyTest, MidpointFreeFunction)
{
   auto coord1 = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 100.0);
   auto coord2 = ECEFCoordinate::from_geodetic_deg(38.0, -121.0, 200.0);

   auto mid1 = coord1.midpoint(coord2);
   auto mid2 = midpoint(coord1, coord2);

   // Both methods should give same result
   EXPECT_TRUE(mid1 == mid2);
}

TEST(GeodesyTest, Interpolate)
{
   auto start = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 100.0);
   auto end = ECEFCoordinate::from_geodetic_deg(38.0, -121.0, 200.0);

   // Test boundary conditions
   auto at_start = interpolate(start, end, 0.0);
   auto at_end = interpolate(start, end, 1.0);
   auto at_mid = interpolate(start, end, 0.5);

   EXPECT_TRUE(at_start.approx_equal(start, 1.0));
   EXPECT_TRUE(at_end.approx_equal(end, 1.0));

   // Midpoint should be equidistant
   double dist_to_mid = start.geodesic_distance_to(at_mid);
   double total_dist = start.geodesic_distance_to(end);
   EXPECT_TRUE(approx_equal(dist_to_mid, total_dist / 2.0, 1000.0));
}

TEST(GeodesyTest, Centroid)
{
   std::vector<ECEFCoordinate> coords = {
       ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 100.0),
       ECEFCoordinate::from_geodetic_deg(37.1, -122.0, 100.0),
       ECEFCoordinate::from_geodetic_deg(37.0, -121.9, 100.0),
       ECEFCoordinate::from_geodetic_deg(37.1, -121.9, 100.0)
   };

   auto center = centroid(coords);

   // Centroid should be near the geometric center
   EXPECT_TRUE(approx_equal(center.latitude_deg(), 37.05, 0.01));
   EXPECT_TRUE(approx_equal(center.longitude_deg(), -121.95, 0.01));
   EXPECT_TRUE(approx_equal(center.altitude(), 100.0, 10.0));
}

TEST(GeodesyTest, CommutativeAddition)
{
   auto coord = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 100.0);
   Vec3 offset{1000, 2000, 500};

   auto result1 = coord + offset;
   auto result2 = offset + coord;

   // Should be commutative
   EXPECT_TRUE(result1 == result2);
}
