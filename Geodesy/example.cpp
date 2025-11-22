// Example usage of the Geodesy library
#include <dstl/Geodesy.h>
#include <iostream>
#include <iomanip>
#include <vector>

using namespace dstl::geo;
using namespace dstl::linalg;

void print_separator() {
    std::cout << "\n" << std::string(60, '=') << "\n\n";
}

int main() {
    std::cout << std::fixed << std::setprecision(6);
    
    // Example 1: Creating coordinates
    std::cout << "Example 1: Creating Coordinates\n";
    print_separator();
    
    // From geodetic (degrees)
    auto san_francisco = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 10.0);
    auto los_angeles = ECEFCoordinate::from_geodetic_deg(34.0522, -118.2437, 50.0);
    
    std::cout << "San Francisco:\n";
    std::cout << "  Lat: " << san_francisco.latitude_deg() << "°\n";
    std::cout << "  Lon: " << san_francisco.longitude_deg() << "°\n";
    std::cout << "  Alt: " << san_francisco.altitude() << " m\n";
    std::cout << "  ECEF: (" << san_francisco.x() << ", " 
              << san_francisco.y() << ", " << san_francisco.z() << ")\n";
    
    // Example 2: Distance calculations (fast - no trig!)
    std::cout << "\nExample 2: Distance Calculations\n";
    print_separator();
    
    double distance = san_francisco.distance_to(los_angeles);
    std::cout << "Distance SF to LA: " << distance / 1000.0 << " km\n";
    
    Vec3 vector = san_francisco.vector_to(los_angeles);
    std::cout << "Vector SF to LA: (" << vector.x << ", " 
              << vector.y << ", " << vector.z << ") m\n";
    
    // Example 3: Cached trigonometry
    std::cout << "\nExample 3: Cached Trigonometry\n";
    print_separator();
    
    std::cout << "San Francisco trig values (cached after first geodetic access):\n";
    std::cout << "  sin(lat): " << san_francisco.sin_lat() << "\n";
    std::cout << "  cos(lat): " << san_francisco.cos_lat() << "\n";
    std::cout << "  sin(lon): " << san_francisco.sin_lon() << "\n";
    std::cout << "  cos(lon): " << san_francisco.cos_lon() << "\n";
    
    // Example 4: ECEF operations (no trig!)
    std::cout << "\nExample 4: ECEF Operations\n";
    print_separator();
    
    // Create a trajectory by adding ECEF offsets
    std::vector<ECEFCoordinate> trajectory;
    trajectory.push_back(san_francisco);
    
    for (int i = 1; i <= 5; ++i) {
        // Move 1km in X direction each step (no trig!)
        trajectory.push_back(trajectory.back() + Vec3{1000, 0, 0});
    }
    
    std::cout << "Trajectory (5 steps of 1km in X direction):\n";
    for (size_t i = 0; i < trajectory.size(); ++i) {
        std::cout << "  Point " << i << ": Lat=" << trajectory[i].latitude_deg() 
                  << "°, Lon=" << trajectory[i].longitude_deg() << "°\n";
    }
    
    // Example 5: Local ENU frame
    std::cout << "\nExample 5: Local ENU Frame\n";
    print_separator();
    
    auto origin = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 10.0);
    auto frame = origin.local_frame();
    
    // Define points in intuitive ENU coordinates
    std::cout << "Points relative to origin (ENU coordinates):\n";
    std::vector<Vec3> enu_points = {
        {100, 0, 0},      // 100m East
        {0, 200, 0},      // 200m North
        {-50, 150, 10},   // 50m West, 150m North, 10m Up
    };
    
    for (size_t i = 0; i < enu_points.size(); ++i) {
        auto coord = frame.enu_to_coord(enu_points[i]);
        std::cout << "  ENU (" << enu_points[i].x << ", " << enu_points[i].y 
                  << ", " << enu_points[i].z << ") -> ";
        std::cout << "Lat=" << coord.latitude_deg() << "°, Lon=" 
                  << coord.longitude_deg() << "°, Alt=" << coord.altitude() << "m\n";
    }
    
    // Example 6: Performance - batch ECEF operations
    std::cout << "\nExample 6: Performance - Batch Operations\n";
    print_separator();
    
    auto start = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 10.0);
    std::vector<ECEFCoordinate> path;
    
    // Create 1000 points (no trig per point!)
    for (int i = 0; i < 1000; ++i) {
        path.push_back(start + Vec3{i * 10.0, i * 5.0, 0});
    }
    
    // Compute total path length (no trig!)
    double total_length = 0;
    for (size_t i = 1; i < path.size(); ++i) {
        total_length += path[i-1].distance_to(path[i]);
    }
    
    std::cout << "Created path with " << path.size() << " points\n";
    std::cout << "Total path length: " << total_length / 1000.0 << " km\n";
    std::cout << "All distance calculations done without trigonometry!\n";
    
    // Example 7: Round-trip verification
    std::cout << "\nExample 7: Round-Trip Verification\n";
    print_separator();
    
    double lat_orig = 45.0;
    double lon_orig = -120.0;
    double alt_orig = 100.0;
    
    auto coord = ECEFCoordinate::from_geodetic_deg(lat_orig, lon_orig, alt_orig);
    
    double lat_back = coord.latitude_deg();
    double lon_back = coord.longitude_deg();
    double alt_back = coord.altitude();
    
    std::cout << "Original: Lat=" << lat_orig << "°, Lon=" << lon_orig 
              << "°, Alt=" << alt_orig << "m\n";
    std::cout << "Round-trip: Lat=" << lat_back << "°, Lon=" << lon_back 
              << "°, Alt=" << alt_back << "m\n";
    std::cout << "Error: Lat=" << (lat_back - lat_orig) << "°, Lon=" 
              << (lon_back - lon_orig) << "°, Alt=" << (alt_back - alt_orig) << "m\n";
    
    // Example 8: Geodesic distance and bearing
    std::cout << "\nExample 8: Geodesic Distance and Bearing\n";
    print_separator();
    
    double chord_dist = san_francisco.distance_to(los_angeles);
    double geodesic_dist = san_francisco.geodesic_distance_to(los_angeles);
    double bearing = san_francisco.bearing_to(los_angeles);
    
    std::cout << "San Francisco to Los Angeles:\n";
    std::cout << "  Chord distance: " << chord_dist / 1000.0 << " km\n";
    std::cout << "  Geodesic distance: " << geodesic_dist / 1000.0 << " km\n";
    std::cout << "  Bearing: " << bearing * 180.0 / M_PI << "° (from North)\n";
    std::cout << "  Difference: " << (geodesic_dist - chord_dist) / 1000.0 << " km\n";
    
    // Example 9: Movement by bearing
    std::cout << "\nExample 9: Movement by Bearing\n";
    print_separator();
    
    auto start_point = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0.0);
    
    std::cout << "Starting point: " << start_point.to_string() << "\n\n";
    std::cout << "Moving 10km in different directions:\n";
    
    std::vector<std::pair<std::string, double>> directions = {
        {"North", 0.0},
        {"East", M_PI / 2.0},
        {"South", M_PI},
        {"West", 3.0 * M_PI / 2.0},
        {"Northeast", M_PI / 4.0}
    };
    
    for (const auto& [name, bearing_rad] : directions) {
        auto target = start_point.move_by_bearing(10000.0, bearing_rad);
        std::cout << "  " << name << ": " << target.to_string() << "\n";
    }
    
    // Example 10: NED frame
    std::cout << "\nExample 10: NED (North-East-Down) Frame\n";
    print_separator();
    
    auto ned_origin = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 10.0);
    auto ned_frame = ned_origin.local_ned_frame();
    
    std::cout << "NED frame at origin:\n";
    std::cout << "  North vector: (" << ned_frame.north.x << ", " 
              << ned_frame.north.y << ", " << ned_frame.north.z << ")\n";
    std::cout << "  East vector: (" << ned_frame.east.x << ", " 
              << ned_frame.east.y << ", " << ned_frame.east.z << ")\n";
    std::cout << "  Down vector: (" << ned_frame.down.x << ", " 
              << ned_frame.down.y << ", " << ned_frame.down.z << ")\n";
    
    // Convert a point to NED
    auto target_point = ECEFCoordinate::from_geodetic_deg(37.7850, -122.4094, 15.0);
    Vec3 ned_coords = ned_frame.coord_to_ned(target_point);
    
    std::cout << "\nTarget point in NED coordinates:\n";
    std::cout << "  North: " << ned_coords.x << " m\n";
    std::cout << "  East: " << ned_coords.y << " m\n";
    std::cout << "  Down: " << ned_coords.z << " m\n";
    
    std::cout << "\nAll examples completed successfully!\n";
    
    return 0;
}
