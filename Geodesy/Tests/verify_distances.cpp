// Quick program to verify actual computed distances
#include <dstl/Geodesy.h>
#include <iostream>
#include <iomanip>

using namespace dstl::geo;

int main()
{
    std::cout << std::fixed << std::setprecision(2);
    
    // Airport coordinates
    struct Airport {
        const char* code;
        double lat, lon, alt;
    };
    
    Airport airports[] = {
        {"SFO", 37.6213, -122.3790, 4.0},
        {"LAX", 33.9416, -118.4085, 38.0},
        {"JFK", 40.6413, -73.7781, 4.0},
        {"ORD", 41.9742, -87.9073, 205.0},
        {"LHR", 51.4700, -0.4543, 25.0},
        {"DXB", 25.2532, 55.3657, 19.0},
        {"SYD", -33.9461, 151.1772, 6.0},
        {"GRU", -23.4356, -46.4731, 749.0},
        {"NRT", 35.7647, 140.3864, 43.0}
    };
    
    // Test pairs
    struct Pair {
        const char* from;
        const char* to;
    } pairs[] = {
        {"SFO", "LAX"},
        {"JFK", "ORD"},
        {"LAX", "JFK"},
        {"JFK", "LHR"},
        {"SFO", "NRT"},
        {"LAX", "SYD"},
        {"DXB", "SFO"},
        {"GRU", "NRT"}
    };
    
    for (const auto& pair : pairs) {
        // Find airports
        Airport *from = nullptr, *to = nullptr;
        for (auto& ap : airports) {
            if (std::string(ap.code) == pair.from) from = &ap;
            if (std::string(ap.code) == pair.to) to = &ap;
        }
        
        if (!from || !to) continue;
        
        auto coord1 = ECEFCoordinate::from_geodetic_deg(from->lat, from->lon, from->alt);
        auto coord2 = ECEFCoordinate::from_geodetic_deg(to->lat, to->lon, to->alt);
        
        double distance = coord1.geodesic_distance_to(coord2);
        double bearing = coord1.bearing_to(coord2) * RAD_TO_DEG;
        
        std::cout << pair.from << " to " << pair.to << ":\n";
        std::cout << "  Distance: " << distance << " m (" << distance/1000.0 << " km)\n";
        std::cout << "  Bearing:  " << bearing << "°\n\n";
    }
    
    return 0;
}
