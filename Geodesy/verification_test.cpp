// Verification Test - Flat-Earth Displacement Accuracy
// Test: Fly 100 NM NE from SF at 30,000 ft, verify altitude and distance preservation

#include <dstl/Geodesy.h>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace dstl::geo;
using namespace dstl::linalg;

int main() {
    std::cout << std::fixed << std::setprecision(6);
    
    // Constants
    constexpr double NM_TO_M = 1852.0;      // Nautical miles to meters
    constexpr double FT_TO_M = 0.3048;      // Feet to meters
    
    std::cout << "=================================================================\n";
    std::cout << "  Verification Test: Flat-Earth Displacement Accuracy\n";
    std::cout << "=================================================================\n\n";
    
    // Starting position: Overhead San Francisco at 30,000 ft
    double sf_lat = 37.7749;
    double sf_lon = -122.4194;
    double altitude_ft = 30000.0;
    double altitude_m = altitude_ft * FT_TO_M;
    
    auto start = ECEFCoordinate::from_geodetic_deg(sf_lat, sf_lon, altitude_m);
    
    std::cout << "Starting Position:\n";
    std::cout << "  Location: San Francisco\n";
    std::cout << "  Latitude: " << start.latitude_deg() << "°\n";
    std::cout << "  Longitude: " << start.longitude_deg() << "°\n";
    std::cout << "  Altitude: " << start.altitude() / FT_TO_M << " ft (" 
              << start.altitude() << " m)\n";
    std::cout << "  ECEF: (" << start.x() << ", " << start.y() << ", " << start.z() << ")\n";
    
    // Maneuver: Fly 100 NM Northeast in local plane
    std::cout << "\n-----------------------------------------------------------------\n";
    std::cout << "Maneuver: Fly 100 NM Northeast in local flat-earth plane\n";
    std::cout << "-----------------------------------------------------------------\n\n";
    
    double distance_nm = 100.0;
    double distance_m = distance_nm * NM_TO_M;
    double bearing_deg = 45.0;  // Northeast
    
    // In local ENU frame: NE = equal East and North components
    double east_m = distance_m / std::sqrt(2.0);
    double north_m = distance_m / std::sqrt(2.0);
    double altitude_change_m = 0.0;  // Maintain altitude
    
    std::cout << "Flat-Earth Displacement (ENU frame):\n";
    std::cout << "  East: " << east_m << " m (" << east_m / NM_TO_M << " NM)\n";
    std::cout << "  North: " << north_m << " m (" << north_m / NM_TO_M << " NM)\n";
    std::cout << "  Up: " << altitude_change_m << " m\n";
    std::cout << "  Horizontal distance: " << distance_m << " m (" << distance_nm << " NM)\n";
    std::cout << "  Bearing: " << bearing_deg << "° from North\n";
    
    // Apply flat-earth displacement with constant altitude
    Vec3 enu_displacement{east_m, north_m, altitude_change_m};
    auto end = start.apply_enu_displacement_constant_altitude(enu_displacement, altitude_change_m);
    
    std::cout << "\n-----------------------------------------------------------------\n";
    std::cout << "Final Position:\n";
    std::cout << "-----------------------------------------------------------------\n\n";
    std::cout << "  Latitude: " << end.latitude_deg() << "°\n";
    std::cout << "  Longitude: " << end.longitude_deg() << "°\n";
    std::cout << "  Altitude: " << end.altitude() / FT_TO_M << " ft (" 
              << end.altitude() << " m)\n";
    std::cout << "  ECEF: (" << end.x() << ", " << end.y() << ", " << end.z() << ")\n";
    
    // Verification
    std::cout << "\n=================================================================\n";
    std::cout << "  VERIFICATION RESULTS\n";
    std::cout << "=================================================================\n\n";
    
    // 1. Altitude preservation
    double altitude_change_ft = (end.altitude() - start.altitude()) / FT_TO_M;
    double altitude_error_ft = std::abs(altitude_change_ft);
    
    std::cout << "1. Altitude Preservation:\n";
    std::cout << "   Start altitude: " << start.altitude() / FT_TO_M << " ft\n";
    std::cout << "   End altitude: " << end.altitude() / FT_TO_M << " ft\n";
    std::cout << "   Change: " << altitude_change_ft << " ft\n";
    std::cout << "   Error: " << altitude_error_ft << " ft\n";
    std::cout << "   Status: " << (altitude_error_ft < 1.0 ? "✓ PASS" : "✗ FAIL") 
              << " (< 1 ft tolerance)\n\n";
    
    // 2. Horizontal distance (chord distance)
    double chord_distance_m = start.distance_to(end);
    double chord_distance_nm = chord_distance_m / NM_TO_M;
    double chord_error_nm = std::abs(chord_distance_nm - distance_nm);
    double chord_error_percent = (chord_error_nm / distance_nm) * 100.0;
    
    std::cout << "2. Horizontal Distance (Chord):\n";
    std::cout << "   Planned: " << distance_nm << " NM\n";
    std::cout << "   Actual: " << chord_distance_nm << " NM\n";
    std::cout << "   Error: " << chord_error_nm << " NM (" << chord_error_percent << "%)\n";
    std::cout << "   Status: " << (chord_error_nm < 0.1 ? "✓ PASS" : "✗ FAIL") 
              << " (< 0.1 NM tolerance)\n\n";
    
    // 3. Geodesic distance (great circle at surface)
    // Note: Geodesic measures surface distance, we're at altitude, so expect slight difference
    double geodesic_distance_m = start.geodesic_distance_to(end);
    double geodesic_distance_nm = geodesic_distance_m / NM_TO_M;
    double geodesic_error_nm = std::abs(geodesic_distance_nm - distance_nm);
    double geodesic_error_percent = (geodesic_error_nm / distance_nm) * 100.0;
    
    std::cout << "3. Geodesic Distance (Great Circle at surface):\n";
    std::cout << "   Planned (at altitude): " << distance_nm << " NM\n";
    std::cout << "   Actual (at surface): " << geodesic_distance_nm << " NM\n";
    std::cout << "   Difference: " << geodesic_error_nm << " NM (" << geodesic_error_percent << "%)\n";
    std::cout << "   Note: Geodesic measures surface distance, we're at " << altitude_ft << " ft\n";
    std::cout << "   Status: " << (geodesic_error_nm < 0.2 ? "✓ PASS" : "✗ FAIL") 
              << " (< 0.2 NM tolerance for altitude effect)\n\n";
    
    // 4. Bearing verification
    double actual_bearing_rad = start.bearing_to(end);
    double actual_bearing_deg = actual_bearing_rad * RAD_TO_DEG;
    double bearing_error_deg = std::abs(actual_bearing_deg - bearing_deg);
    
    std::cout << "4. Bearing:\n";
    std::cout << "   Planned: " << bearing_deg << "° from North\n";
    std::cout << "   Actual: " << actual_bearing_deg << "° from North\n";
    std::cout << "   Error: " << bearing_error_deg << "°\n";
    std::cout << "   Status: " << (bearing_error_deg < 0.1 ? "✓ PASS" : "✗ FAIL") 
              << " (< 0.1° tolerance)\n\n";
    
    // 5. Position change
    double lat_change = end.latitude_deg() - start.latitude_deg();
    double lon_change = end.longitude_deg() - start.longitude_deg();
    
    std::cout << "5. Position Change:\n";
    std::cout << "   Latitude: " << lat_change << "° (" 
              << lat_change * 60.0 << " arc-minutes)\n";
    std::cout << "   Longitude: " << lon_change << "° (" 
              << lon_change * 60.0 << " arc-minutes)\n\n";
    
    // Overall summary
    std::cout << "=================================================================\n";
    std::cout << "  OVERALL SUMMARY\n";
    std::cout << "=================================================================\n\n";
    
    bool altitude_ok = altitude_error_ft < 1.0;
    bool chord_ok = chord_error_nm < 0.1;
    bool geodesic_ok = geodesic_error_nm < 0.2;  // Relaxed for altitude effect
    bool bearing_ok = bearing_error_deg < 0.1;
    bool all_pass = altitude_ok && chord_ok && geodesic_ok && bearing_ok;
    
    std::cout << "Test Results:\n";
    std::cout << "  Altitude preservation: " << (altitude_ok ? "✓ PASS" : "✗ FAIL") << "\n";
    std::cout << "  Chord distance: " << (chord_ok ? "✓ PASS" : "✗ FAIL") << "\n";
    std::cout << "  Geodesic distance: " << (geodesic_ok ? "✓ PASS" : "✗ FAIL") << "\n";
    std::cout << "  Bearing: " << (bearing_ok ? "✓ PASS" : "✗ FAIL") << "\n";
    std::cout << "\n";
    
    if (all_pass) {
        std::cout << "✓✓✓ ALL TESTS PASSED ✓✓✓\n\n";
        std::cout << "Conclusion:\n";
        std::cout << "  • Altitude is preserved exactly (< 1 ft error)\n";
        std::cout << "  • Horizontal distance matches flat-earth calculation\n";
        std::cout << "  • Great circle distance is accurate\n";
        std::cout << "  • Bearing is maintained correctly\n";
        std::cout << "  • Flat-earth displacement works perfectly!\n";
    } else {
        std::cout << "✗✗✗ SOME TESTS FAILED ✗✗✗\n";
    }
    
    std::cout << "\n=================================================================\n\n";
    
    return all_pass ? 0 : 1;
}
