// Comprehensive Example - All Geodesy Features
// Demonstrates the complete API of the Geodesy module

#include <dstl/Geodesy.h>
#include <iostream>
#include <iomanip>
#include <vector>

using namespace dstl::geo;
using namespace dstl::linalg;

void print_separator(const std::string& title) {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  " << title << "\n";
    std::cout << std::string(70, '=') << "\n\n";
}

int main() {
    std::cout << std::fixed << std::setprecision(6);
    
    print_separator("DSTL Geodesy - Comprehensive Feature Demonstration");
    
    // ========================================================================
    // 1. Construction and Access
    // ========================================================================
    print_separator("1. Construction and Access");
    
    // From geodetic (degrees)
    auto sf = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 10.0);
    std::cout << "San Francisco: " << sf << "\n";
    std::cout << "  ECEF: " << sf.to_string_ecef() << "\n";
    std::cout << "  Cached trig: sin(lat)=" << sf.sin_lat() 
              << ", cos(lat)=" << sf.cos_lat() << "\n";
    
    // From ECEF (fast)
    auto from_ecef = ECEFCoordinate(sf.ecef());
    std::cout << "\nFrom ECEF: " << from_ecef << "\n";
    
    // ========================================================================
    // 2. Distance and Bearing
    // ========================================================================
    print_separator("2. Distance and Bearing");
    
    auto la = ECEFCoordinate::from_geodetic_deg(34.0522, -118.2437, 50.0);
    
    double chord = sf.distance_to(la);
    double geodesic = sf.geodesic_distance_to(la);
    double bearing = sf.bearing_to(la) * RAD_TO_DEG;
    
    std::cout << "San Francisco to Los Angeles:\n";
    std::cout << "  Chord distance: " << chord / 1000.0 << " km\n";
    std::cout << "  Geodesic distance: " << geodesic / 1000.0 << " km\n";
    std::cout << "  Bearing: " << bearing << "° from North\n";
    std::cout << "  Difference: " << (geodesic - chord) / 1000.0 << " km\n";
    
    // ========================================================================
    // 3. Movement and Navigation
    // ========================================================================
    print_separator("3. Movement and Navigation");
    
    // Move by bearing
    auto target1 = sf.move_by_bearing(50000.0, 45.0 * DEG_TO_RAD);
    std::cout << "Move 50km NE from SF: " << target1 << "\n";
    
    // Accurate movement for long distances
    auto target2 = sf.move_by_bearing_accurate(500000.0, 90.0 * DEG_TO_RAD);
    std::cout << "Move 500km East from SF: " << target2 << "\n";
    
    // ========================================================================
    // 4. Local Frames (ENU and NED)
    // ========================================================================
    print_separator("4. Local Frames");
    
    auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 1000.0);
    
    // ENU frame
    auto enu_frame = origin.local_frame();
    Vec3 enu_offset{10000, 20000, 500};  // 10km E, 20km N, 500m Up
    auto enu_target = enu_frame.enu_to_coord(enu_offset);
    std::cout << "ENU displacement (10km E, 20km N, 500m Up):\n";
    std::cout << "  Result: " << enu_target << "\n";
    
    // NED frame
    auto ned_frame = origin.local_ned_frame();
    Vec3 ned_offset{20000, 10000, -500};  // 20km N, 10km E, 500m Up
    auto ned_target = ned_frame.ned_to_coord(ned_offset);
    std::cout << "\nNED displacement (20km N, 10km E, 500m Up):\n";
    std::cout << "  Result: " << ned_target << "\n";
    
    // ========================================================================
    // 5. Flat-Earth Dynamics
    // ========================================================================
    print_separator("5. Flat-Earth Dynamics");
    
    constexpr double NM_TO_M = 1852.0;
    constexpr double FT_TO_M = 0.3048;
    
    // Single displacement
    Vec3 flight_disp{50.0 * NM_TO_M / std::sqrt(2.0),  // 50 NM NE
                     50.0 * NM_TO_M / std::sqrt(2.0),
                     10000.0 * FT_TO_M};                // Climb 10,000 ft
    auto flight_pos = origin.apply_enu_displacement(flight_disp);
    std::cout << "Flight maneuver (50 NM NE, climb 10,000 ft):\n";
    std::cout << "  Final: " << flight_pos << "\n";
    
    // Multi-leg trajectory
    std::vector<Vec3> legs = {
        {30000, 30000, 1000},   // 30km E, 30km N, climb 1km
        {20000, -10000, -500},  // 20km E, 10km S, descend 500m
    };
    auto traj_end = origin.apply_enu_trajectory(legs);
    std::cout << "\nMulti-leg trajectory:\n";
    std::cout << "  Final: " << traj_end << "\n";
    
    // ========================================================================
    // 6. Geometric Utilities
    // ========================================================================
    print_separator("6. Geometric Utilities");
    
    // Surface normal
    auto normal = sf.surface_normal();
    std::cout << "Surface normal at SF: (" << normal.x << ", " 
              << normal.y << ", " << normal.z << ")\n";
    std::cout << "  Magnitude: " << normal.norm() << " (should be 1.0)\n";
    
    // Midpoint
    auto mid = midpoint(sf, la);
    std::cout << "\nMidpoint between SF and LA: " << mid << "\n";
    double dist_sf = sf.geodesic_distance_to(mid);
    double dist_la = la.geodesic_distance_to(mid);
    std::cout << "  Distance from SF: " << dist_sf / 1000.0 << " km\n";
    std::cout << "  Distance from LA: " << dist_la / 1000.0 << " km\n";
    
    // Interpolation
    auto quarter = interpolate(sf, la, 0.25);
    auto three_quarter = interpolate(sf, la, 0.75);
    std::cout << "\n25% along SF to LA: " << quarter << "\n";
    std::cout << "75% along SF to LA: " << three_quarter << "\n";
    
    // Centroid
    std::vector<ECEFCoordinate> cities = {
        ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 0),  // SF
        ECEFCoordinate::from_geodetic_deg(34.0522, -118.2437, 0),  // LA
        ECEFCoordinate::from_geodetic_deg(32.7157, -117.1611, 0),  // San Diego
        ECEFCoordinate::from_geodetic_deg(36.7783, -119.4179, 0),  // Fresno
    };
    auto center = centroid(cities);
    std::cout << "\nCentroid of 4 California cities: " << center << "\n";
    
    // ========================================================================
    // 7. Comparison and Equality
    // ========================================================================
    print_separator("7. Comparison and Equality");
    
    auto coord1 = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 100.0);
    auto coord2 = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 100.0);
    auto coord3 = ECEFCoordinate::from_geodetic_deg(37.0001, -122.0, 100.0);
    
    std::cout << "coord1 == coord2: " << (coord1 == coord2 ? "true" : "false") << "\n";
    std::cout << "coord1 == coord3: " << (coord1 == coord3 ? "true" : "false") << "\n";
    std::cout << "coord1.approx_equal(coord3, 100.0): " 
              << (coord1.approx_equal(coord3, 100.0) ? "true" : "false") << "\n";
    
    // ========================================================================
    // 8. Vector Operations
    // ========================================================================
    print_separator("8. Vector Operations");
    
    Vec3 offset{1000, 2000, 500};
    auto result1 = coord1 + offset;
    auto result2 = offset + coord1;  // Commutative
    
    std::cout << "coord + offset: " << result1 << "\n";
    std::cout << "offset + coord: " << result2 << "\n";
    std::cout << "Are they equal? " << (result1 == result2 ? "yes" : "no") << "\n";
    
    Vec3 vec = coord1.vector_to(coord3);
    std::cout << "\nVector from coord1 to coord3: (" << vec.x << ", " 
              << vec.y << ", " << vec.z << ")\n";
    
    // ========================================================================
    // 9. Container Support
    // ========================================================================
    print_separator("9. Container Support");
    
    std::unordered_map<ECEFCoordinate, std::string> locations;
    locations[sf] = "San Francisco";
    locations[la] = "Los Angeles";
    locations[center] = "California Center";
    
    std::cout << "Locations in map:\n";
    for (const auto& [coord, name] : locations) {
        std::cout << "  " << name << ": " << coord.latitude_deg() 
                  << "°, " << coord.longitude_deg() << "°\n";
    }
    
    // ========================================================================
    // 10. Performance Demonstration
    // ========================================================================
    print_separator("10. Performance Demonstration");
    
    std::cout << "Creating 10,000 coordinates and computing distances...\n";
    
    auto ref = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0);
    std::vector<ECEFCoordinate> points;
    
    for (int i = 0; i < 10000; ++i) {
        double lat = 37.0 + (i % 100) * 0.01;
        double lon = -122.0 + (i / 100) * 0.01;
        points.push_back(ECEFCoordinate::from_geodetic_deg(lat, lon, 0));
    }
    
    // Compute all distances (fast - no trig per distance!)
    double total_dist = 0;
    for (const auto& point : points) {
        total_dist += ref.distance_to(point);
    }
    
    std::cout << "  Computed 10,000 distances\n";
    std::cout << "  Average distance: " << total_dist / points.size() / 1000.0 << " km\n";
    std::cout << "  All done with zero trigonometry per distance!\n";
    
    // ========================================================================
    print_separator("All Features Demonstrated Successfully!");
    
    std::cout << "\nKey Takeaways:\n";
    std::cout << "  ✓ ECEF primary storage for fast vector operations\n";
    std::cout << "  ✓ Lazy geodetic computation with caching\n";
    std::cout << "  ✓ Cached trigonometry for repeated use\n";
    std::cout << "  ✓ Sub-millimeter accuracy with Vincenty formulas\n";
    std::cout << "  ✓ Real-time performance for all operations\n";
    std::cout << "  ✓ Complete flat-earth dynamics support\n";
    std::cout << "  ✓ Rich geometric utilities\n";
    std::cout << "  ✓ Container support with hashing\n";
    std::cout << "\nReady for production use in aerospace, robotics, and navigation!\n";
    
    return 0;
}
