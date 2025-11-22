// Flight Dynamics Example - Flat Earth to WGS84 Mapping
// Demonstrates computing maneuvers in local flat-earth frame and mapping back to WGS84

#include <dstl/Geodesy.h>
#include <iostream>
#include <iomanip>
#include <vector>

using namespace dstl::geo;
using namespace dstl::linalg;

// Conversion constants
constexpr double NM_TO_M = 1852.0;      // Nautical miles to meters
constexpr double FT_TO_M = 0.3048;      // Feet to meters
constexpr double KT_TO_MS = 0.514444;   // Knots to m/s

void print_separator() {
    std::cout << "\n" << std::string(70, '=') << "\n\n";
}

void print_position(const std::string& label, const ECEFCoordinate& pos) {
    std::cout << label << ":\n";
    std::cout << "  Lat: " << std::fixed << std::setprecision(6) 
              << pos.latitude_deg() << "°\n";
    std::cout << "  Lon: " << pos.longitude_deg() << "°\n";
    std::cout << "  Alt: " << std::setprecision(1) 
              << pos.altitude() * 3.28084 << " ft ("
              << pos.altitude() << " m)\n";
}

int main() {
    std::cout << std::fixed << std::setprecision(3);
    
    // ========================================================================
    // Example 1: Simple Flat-Earth Maneuver
    // ========================================================================
    std::cout << "Example 1: Simple Flat-Earth Maneuver\n";
    print_separator();
    
    // Starting position: San Francisco International (KSFO)
    auto ksfo = ECEFCoordinate::from_geodetic_deg(37.6213, -122.3790, 13.0);
    print_position("Starting at KSFO", ksfo);
    
    // Fly 100 nautical miles Northeast at 30,000 ft, then descend 10,000 ft
    std::cout << "\nManeuver: 100 NM NE at FL300, descend 10,000 ft\n";
    
    // In flat-earth ENU coordinates:
    // - Northeast = equal East and North components
    // - Start at 30,000 ft, descend to 20,000 ft
    double distance_m = 100.0 * NM_TO_M;
    double east_m = distance_m / std::sqrt(2.0);   // NE = 45°
    double north_m = distance_m / std::sqrt(2.0);
    double climb_to_cruise = 30000.0 * FT_TO_M - ksfo.altitude();
    double descent = -10000.0 * FT_TO_M;
    
    // Apply flat-earth displacement
    Vec3 enu_displacement{east_m, north_m, climb_to_cruise + descent};
    auto final_pos = ksfo.apply_enu_displacement(enu_displacement);
    
    print_position("\nFinal position", final_pos);
    
    // Verify the maneuver
    double actual_distance = ksfo.geodesic_distance_to(final_pos);
    double actual_bearing = ksfo.bearing_to(final_pos) * RAD_TO_DEG;
    
    std::cout << "\nVerification:\n";
    std::cout << "  Planned distance: " << 100.0 << " NM\n";
    std::cout << "  Actual distance: " << actual_distance / NM_TO_M << " NM\n";
    std::cout << "  Planned bearing: 45.0°\n";
    std::cout << "  Actual bearing: " << actual_bearing << "°\n";
    std::cout << "  Altitude change: " << (final_pos.altitude() - ksfo.altitude()) * 3.28084 
              << " ft\n";
    
    // ========================================================================
    // Example 2: Multi-Leg Flight Plan
    // ========================================================================
    std::cout << "\n\nExample 2: Multi-Leg Flight Plan\n";
    print_separator();
    
    // Starting position: Los Angeles (KLAX)
    auto klax = ECEFCoordinate::from_geodetic_deg(33.9425, -118.4081, 125.0);
    print_position("Starting at KLAX", klax);
    
    std::cout << "\nFlight plan (in flat-earth NED coordinates):\n";
    std::cout << "  Leg 1: Climb to FL350, fly 150 NM North\n";
    std::cout << "  Leg 2: Turn East, fly 100 NM\n";
    std::cout << "  Leg 3: Descend to FL250, fly 80 NM Northeast\n";
    std::cout << "  Leg 4: Final descent to 5000 ft, fly 50 NM North\n";
    
    // Define legs in NED (North-East-Down) - common in aerospace
    std::vector<Vec3> ned_legs;
    
    // Leg 1: Climb to FL350 (35,000 ft), fly 150 NM North
    double climb_1 = -(35000.0 * FT_TO_M - klax.altitude());  // Negative = up
    ned_legs.push_back(Vec3{150.0 * NM_TO_M, 0.0, climb_1});
    
    // Leg 2: Maintain altitude, fly 100 NM East
    ned_legs.push_back(Vec3{0.0, 100.0 * NM_TO_M, 0.0});
    
    // Leg 3: Descend to FL250, fly 80 NM Northeast
    double descent_1 = 10000.0 * FT_TO_M;  // Positive = down
    double ne_dist = 80.0 * NM_TO_M;
    ned_legs.push_back(Vec3{ne_dist / std::sqrt(2.0), ne_dist / std::sqrt(2.0), descent_1});
    
    // Leg 4: Descend to 5000 ft, fly 50 NM North
    double descent_2 = (25000.0 - 5000.0) * FT_TO_M;
    ned_legs.push_back(Vec3{50.0 * NM_TO_M, 0.0, descent_2});
    
    // Apply trajectory with frame updates every 50 km
    auto final_flight_pos = klax.apply_ned_trajectory(ned_legs, 50000.0);
    
    print_position("\nFinal position after flight plan", final_flight_pos);
    
    // Calculate total distance
    double total_planned = (150.0 + 100.0 + 80.0 + 50.0) * NM_TO_M;
    double total_actual = klax.geodesic_distance_to(final_flight_pos);
    
    std::cout << "\nFlight summary:\n";
    std::cout << "  Total planned distance: " << total_planned / NM_TO_M << " NM\n";
    std::cout << "  Actual ground distance: " << total_actual / NM_TO_M << " NM\n";
    std::cout << "  Starting altitude: " << klax.altitude() * 3.28084 << " ft\n";
    std::cout << "  Final altitude: " << final_flight_pos.altitude() * 3.28084 << " ft\n";
    
    // ========================================================================
    // Example 3: Dynamic Maneuver Simulation
    // ========================================================================
    std::cout << "\n\nExample 3: Dynamic Maneuver Simulation\n";
    print_separator();
    
    // Starting position: Denver (KDEN) - high altitude airport
    auto kden = ECEFCoordinate::from_geodetic_deg(39.8561, -104.6737, 1655.0);
    print_position("Starting at KDEN", kden);
    
    std::cout << "\nSimulating dynamic maneuver:\n";
    std::cout << "  - Takeoff and climb at 250 knots\n";
    std::cout << "  - 5-minute climb to 10,000 ft\n";
    std::cout << "  - Heading: 090° (East)\n";
    
    // Simulation parameters
    double airspeed_ms = 250.0 * KT_TO_MS;
    double climb_rate_fpm = 2000.0;  // feet per minute
    double climb_rate_ms = climb_rate_fpm * FT_TO_M / 60.0;
    double duration_s = 5.0 * 60.0;  // 5 minutes
    double dt = 1.0;  // 1 second time steps
    
    // Simulate in flat-earth NED frame
    ECEFCoordinate current_pos = kden;
    std::vector<Vec3> trajectory_points;
    
    std::cout << "\nSimulating " << duration_s << " seconds of flight...\n";
    
    for (double t = 0; t < duration_s; t += dt) {
        // Compute displacement in this time step (NED frame)
        double north_disp = 0.0;  // Heading 090° = pure East
        double east_disp = airspeed_ms * dt;
        double down_disp = -climb_rate_ms * dt;  // Negative = climbing
        
        Vec3 ned_step{north_disp, east_disp, down_disp};
        current_pos = current_pos.apply_ned_displacement(ned_step);
        
        // Store every 30 seconds
        if (static_cast<int>(t) % 30 == 0) {
            trajectory_points.push_back(current_pos.ecef());
        }
    }
    
    print_position("\nFinal position after maneuver", current_pos);
    
    // Calculate actual performance
    double ground_distance = kden.geodesic_distance_to(current_pos);
    double altitude_gain = current_pos.altitude() - kden.altitude();
    double maneuver_bearing = kden.bearing_to(current_pos) * RAD_TO_DEG;
    
    std::cout << "\nManeuver results:\n";
    std::cout << "  Ground distance: " << ground_distance / NM_TO_M << " NM\n";
    std::cout << "  Altitude gain: " << altitude_gain * 3.28084 << " ft\n";
    std::cout << "  Average ground speed: " 
              << (ground_distance / duration_s) / KT_TO_MS << " knots\n";
    std::cout << "  Actual bearing: " << maneuver_bearing << "°\n";
    std::cout << "  Trajectory points captured: " << trajectory_points.size() << "\n";
    
    // ========================================================================
    // Example 4: Accuracy Comparison - Flat Earth vs Geodesic
    // ========================================================================
    std::cout << "\n\nExample 4: Accuracy Analysis\n";
    print_separator();
    
    auto origin = ECEFCoordinate::from_geodetic_deg(40.0, -100.0, 10000.0 * FT_TO_M);
    
    std::cout << "Comparing flat-earth vs geodesic for various distances:\n\n";
    std::cout << std::setw(15) << "Distance (NM)" 
              << std::setw(20) << "Flat-Earth (NM)" 
              << std::setw(20) << "Geodesic (NM)"
              << std::setw(15) << "Error (m)\n";
    std::cout << std::string(70, '-') << "\n";
    
    std::vector<double> test_distances = {10, 50, 100, 200, 500, 1000};
    
    for (double dist_nm : test_distances) {
        // Apply flat-earth displacement (due East)
        Vec3 enu_disp{dist_nm * NM_TO_M, 0.0, 0.0};
        auto flat_pos = origin.apply_enu_displacement(enu_disp);
        
        // Measure actual geodesic distance
        double geodesic_dist = origin.geodesic_distance_to(flat_pos);
        double error = std::abs(geodesic_dist - dist_nm * NM_TO_M);
        
        std::cout << std::setw(15) << dist_nm
                  << std::setw(20) << dist_nm
                  << std::setw(20) << geodesic_dist / NM_TO_M
                  << std::setw(15) << error << "\n";
    }
    
    std::cout << "\nConclusion: Flat-earth approximation is excellent for < 100 NM\n";
    std::cout << "For longer distances, use frame updates (apply_enu_trajectory)\n";
    
    std::cout << "\n\nAll examples completed successfully!\n";
    
    return 0;
}
