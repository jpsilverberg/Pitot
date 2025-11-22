// Temporary test to print actual computed distances
#include <dstl/Geodesy.h>
#include <gtest/gtest.h>
#include <iostream>
#include <iomanip>

using namespace dstl::geo;

TEST(PrintDistances, ComputeActualValues)
{
    std::cout << std::fixed << std::setprecision(2) << "\n";
    
    // SFO to LAX
    {
        auto sfo = ECEFCoordinate::from_geodetic_deg(37.6213, -122.3790, 4.0);
        auto lax = ECEFCoordinate::from_geodetic_deg(33.9416, -118.4085, 38.0);
        double dist = sfo.geodesic_distance_to(lax);
        double bear = sfo.bearing_to(lax) * RAD_TO_DEG;
        std::cout << "SFO-LAX: " << dist << " m (" << dist/1000.0 << " km), bearing " << bear << "°\n";
    }
    
    // JFK to ORD
    {
        auto jfk = ECEFCoordinate::from_geodetic_deg(40.6413, -73.7781, 4.0);
        auto ord = ECEFCoordinate::from_geodetic_deg(41.9742, -87.9073, 205.0);
        double dist = jfk.geodesic_distance_to(ord);
        double bear = jfk.bearing_to(ord) * RAD_TO_DEG;
        std::cout << "JFK-ORD: " << dist << " m (" << dist/1000.0 << " km), bearing " << bear << "°\n";
    }
    
    // LAX to JFK
    {
        auto lax = ECEFCoordinate::from_geodetic_deg(33.9416, -118.4085, 38.0);
        auto jfk = ECEFCoordinate::from_geodetic_deg(40.6413, -73.7781, 4.0);
        double dist = lax.geodesic_distance_to(jfk);
        double bear = lax.bearing_to(jfk) * RAD_TO_DEG;
        std::cout << "LAX-JFK: " << dist << " m (" << dist/1000.0 << " km), bearing " << bear << "°\n";
    }
    
    // JFK to LHR
    {
        auto jfk = ECEFCoordinate::from_geodetic_deg(40.6413, -73.7781, 4.0);
        auto lhr = ECEFCoordinate::from_geodetic_deg(51.4700, -0.4543, 25.0);
        double dist = jfk.geodesic_distance_to(lhr);
        double bear = jfk.bearing_to(lhr) * RAD_TO_DEG;
        std::cout << "JFK-LHR: " << dist << " m (" << dist/1000.0 << " km), bearing " << bear << "°\n";
    }
    
    // SFO to NRT
    {
        auto sfo = ECEFCoordinate::from_geodetic_deg(37.6213, -122.3790, 4.0);
        auto nrt = ECEFCoordinate::from_geodetic_deg(35.7647, 140.3864, 43.0);
        double dist = sfo.geodesic_distance_to(nrt);
        double bear = sfo.bearing_to(nrt) * RAD_TO_DEG;
        std::cout << "SFO-NRT: " << dist << " m (" << dist/1000.0 << " km), bearing " << bear << "°\n";
    }
    
    // LAX to SYD
    {
        auto lax = ECEFCoordinate::from_geodetic_deg(33.9416, -118.4085, 38.0);
        auto syd = ECEFCoordinate::from_geodetic_deg(-33.9461, 151.1772, 6.0);
        double dist = lax.geodesic_distance_to(syd);
        double bear = lax.bearing_to(syd) * RAD_TO_DEG;
        std::cout << "LAX-SYD: " << dist << " m (" << dist/1000.0 << " km), bearing " << bear << "°\n";
    }
    
    // DXB to SFO
    {
        auto dxb = ECEFCoordinate::from_geodetic_deg(25.2532, 55.3657, 19.0);
        auto sfo = ECEFCoordinate::from_geodetic_deg(37.6213, -122.3790, 4.0);
        double dist = dxb.geodesic_distance_to(sfo);
        double bear = dxb.bearing_to(sfo) * RAD_TO_DEG;
        std::cout << "DXB-SFO: " << dist << " m (" << dist/1000.0 << " km), bearing " << bear << "°\n";
    }
    
    // GRU to NRT
    {
        auto gru = ECEFCoordinate::from_geodetic_deg(-23.4356, -46.4731, 749.0);
        auto nrt = ECEFCoordinate::from_geodetic_deg(35.7647, 140.3864, 43.0);
        double dist = gru.geodesic_distance_to(nrt);
        double bear = gru.bearing_to(nrt) * RAD_TO_DEG;
        std::cout << "GRU-NRT: " << dist << " m (" << dist/1000.0 << " km), bearing " << bear << "°\n";
    }
    
    std::cout << "\n";
}
