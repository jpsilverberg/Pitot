// Debug test to see actual values
#include <dstl/Geodesy.h>
#include <gtest/gtest.h>
#include <iostream>
#include <iomanip>

using namespace dstl::geo;
using namespace dstl::linalg;

TEST(DebugTest, CheckFailingValues)
{
    std::cout << std::fixed << std::setprecision(10) << "\n";
    
    // 1. Antipodal distance
    {
        auto coord1 = ECEFCoordinate::from_geodetic_deg(0, 0, 0);
        auto coord2 = ECEFCoordinate::from_geodetic_deg(0, 179, 0);
        double geodesic_dist = coord1.geodesic_distance_to(coord2);
        double half_circumference = PI * WGS84::a;
        std::cout << "1. Antipodal (0,0) to (0,179):\n";
        std::cout << "   Actual: " << geodesic_dist << " m\n";
        std::cout << "   Expected (half circ): " << half_circumference << " m\n";
        std::cout << "   Diff: " << (geodesic_dist - half_circumference) << " m\n\n";
    }
    
    // 2. Bearing North
    {
        auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0);
        auto north = ECEFCoordinate::from_geodetic_deg(38.0, -122.0, 0);
        double bearing = origin.bearing_to(north);
        std::cout << "2. Bearing North:\n";
        std::cout << "   Actual: " << bearing << " rad (" << (bearing * RAD_TO_DEG) << "°)\n";
        std::cout << "   Expected: 0.0 rad (0°)\n";
        std::cout << "   Diff: " << bearing << " rad\n\n";
    }
    
    // 3. ApproxEqual
    {
        auto coord1 = ECEFCoordinate(1000000.0, 2000000.0, 3000000.0);
        auto coord3 = ECEFCoordinate(1000010.0, 2000000.0, 3000000.0);
        bool result = coord1.approx_equal(coord3, 10.0);
        double actual_dist = coord1.distance_to(coord3);
        std::cout << "3. ApproxEqual (10m tolerance):\n";
        std::cout << "   Actual distance: " << actual_dist << " m\n";
        std::cout << "   Tolerance: 10.0 m\n";
        std::cout << "   Result: " << (result ? "true" : "false") << "\n\n";
    }
    
    // 4. NED displacement altitude
    {
        auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 1000.0);
        Vec3 ned_disp{10000.0, 5000.0, 500.0};
        auto target = origin.apply_ned_displacement(ned_disp);
        double alt_change = target.altitude() - origin.altitude();
        std::cout << "4. NED Displacement (down=500m):\n";
        std::cout << "   Origin alt: " << origin.altitude() << " m\n";
        std::cout << "   Target alt: " << target.altitude() << " m\n";
        std::cout << "   Actual change: " << alt_change << " m\n";
        std::cout << "   Expected: -500.0 m\n";
        std::cout << "   Diff: " << (alt_change + 500.0) << " m\n\n";
    }
    
    // 5. ENU Trajectory altitude
    {
        auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 1000.0);
        std::vector<Vec3> legs = {
            {10000.0, 10000.0, 500.0},
            {5000.0, -5000.0, -200.0},
            {-3000.0, 8000.0, 100.0}
        };
        auto final_pos = origin.apply_enu_trajectory(legs);
        double total_alt_change = 500.0 - 200.0 + 100.0;
        double actual_alt_change = final_pos.altitude() - origin.altitude();
        std::cout << "5. ENU Trajectory:\n";
        std::cout << "   Origin alt: " << origin.altitude() << " m\n";
        std::cout << "   Final alt: " << final_pos.altitude() << " m\n";
        std::cout << "   Actual change: " << actual_alt_change << " m\n";
        std::cout << "   Expected: " << total_alt_change << " m\n";
        std::cout << "   Diff: " << (actual_alt_change - total_alt_change) << " m\n\n";
    }
    
    // 6. NED Trajectory altitude
    {
        auto origin = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 1000.0);
        std::vector<Vec3> legs = {
            {10000.0, 5000.0, -500.0},
            {-5000.0, 10000.0, 200.0}
        };
        auto final_pos = origin.apply_ned_trajectory(legs);
        double total_alt_change = 500.0 - 200.0;
        double actual_alt_change = final_pos.altitude() - origin.altitude();
        std::cout << "6. NED Trajectory:\n";
        std::cout << "   Origin alt: " << origin.altitude() << " m\n";
        std::cout << "   Final alt: " << final_pos.altitude() << " m\n";
        std::cout << "   Actual change: " << actual_alt_change << " m\n";
        std::cout << "   Expected: " << total_alt_change << " m\n";
        std::cout << "   Diff: " << (actual_alt_change - total_alt_change) << " m\n\n";
    }
}
