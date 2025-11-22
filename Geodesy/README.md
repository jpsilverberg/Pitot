# Geodesy Library

A high-performance C++17 library for WGS84 geodetic coordinate transformations with minimal trigonometric operations.

## Overview

The Geodesy component provides efficient handling of Earth-centered, Earth-fixed (ECEF) coordinates with lazy conversion to geodetic (latitude/longitude/altitude) coordinates. The design prioritizes performance by:

- **Primary ECEF storage**: Fast vector operations using `dstl::linalg::Vec3`
- **Lazy geodetic computation**: Lat/lon/alt computed only when accessed
- **Cached trigonometric values**: sin/cos values cached to avoid recomputation
- **Zero-cost ECEF operations**: Distance, vectors, and offsets use pure Vec3 math

## Features

- **WGS84 ECEF coordinates**: Earth-Centered, Earth-Fixed coordinate system
- **Geodetic conversions**: Bidirectional conversion between ECEF and lat/lon/alt
- **Lazy evaluation**: Geodetic coordinates computed on-demand and cached
- **Trig caching**: sin/cos values cached for reuse
- **Local frames**: ENU (East-North-Up) tangent plane coordinate systems
- **Vector operations**: Fast distance, bearing, and offset calculations
- **Geodesic operations**: Great circle distance and bearing using Vincenty's formula
- **Navigation**: Move by distance and bearing along geodesics
- **Header-only**: No compilation required, just include and use

## Quick Start

```cpp
#include <dstl/Geodesy.h>

using namespace dstl::geo;
using namespace dstl::linalg;

// Create from geodetic coordinates (degrees)
auto sf = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 10.0);
auto la = ECEFCoordinate::from_geodetic_deg(34.0522, -118.2437, 50.0);

// Fast ECEF operations (no trig!)
double chord_dist = sf.distance_to(la);  // Straight-line distance
double geodesic_dist = sf.geodesic_distance_to(la);  // Surface distance ~559 km
Vec3 vector = sf.vector_to(la);

// Navigation operations
double bearing = sf.bearing_to(la);  // Bearing from North (radians)
auto target = sf.move_by_bearing(10000.0, bearing);  // Move 10km along bearing

// Lazy geodetic access (trig only on first access)
double lat = sf.latitude_deg();
double lon = sf.longitude_deg();
double alt = sf.altitude();

// Cached trig values
double sin_lat = sf.sin_lat();  // Cached, no recomputation
double cos_lat = sf.cos_lat();

// Local ENU frame for relative positioning
auto frame = sf.local_frame();
Vec3 enu_offset{100, 200, 50};  // 100m E, 200m N, 50m U
auto target = frame.enu_to_coord(enu_offset);
```

## Installation & Building

### As Part of DSTL

```bash
cd DSTL
mkdir build && cd build
cmake ..
make
```

### Standalone

```bash
cd Geodesy
mkdir build && cd build
cmake ..
make
```

### Running Tests

```bash
cd build
./Tests/GeodesyTests
```

## Usage Guide

### Creating Coordinates

#### From ECEF (Fast - No Trig)

```cpp
// Direct ECEF construction
ECEFCoordinate coord1(4000000, 3000000, 2000000);  // X, Y, Z in meters

// From Vec3
Vec3 ecef{4000000, 3000000, 2000000};
ECEFCoordinate coord2(ecef);
```

#### From Geodetic (Requires Trig)

```cpp
// From radians
double lat_rad = 0.7854;  // 45°
double lon_rad = -2.0944; // -120°
double alt_m = 100.0;
auto coord = ECEFCoordinate::from_geodetic(lat_rad, lon_rad, alt_m);

// From degrees (more convenient)
auto coord = ECEFCoordinate::from_geodetic_deg(45.0, -120.0, 100.0);
```

### Accessing Coordinates

#### ECEF Access (Fast - No Trig)

```cpp
const Vec3& ecef = coord.ecef();
double x = coord.x();
double y = coord.y();
double z = coord.z();
```

#### Geodetic Access (Lazy Trig)

```cpp
// First access computes and caches
double lat_rad = coord.latitude();
double lon_rad = coord.longitude();
double alt_m = coord.altitude();

// Subsequent accesses use cache (O(1))
double lat_rad2 = coord.latitude();  // Cache hit!

// Degrees
double lat_deg = coord.latitude_deg();
double lon_deg = coord.longitude_deg();
```

#### Cached Trigonometric Values

```cpp
// Access cached sin/cos (computed with geodetic)
double sin_lat = coord.sin_lat();
double cos_lat = coord.cos_lat();
double sin_lon = coord.sin_lon();
double cos_lon = coord.cos_lon();

// These are cached - no recomputation!
double sin_lat2 = coord.sin_lat();  // O(1)
```

### Vector Operations (Fast - No Trig)

```cpp
auto coord1 = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 10.0);
auto coord2 = ECEFCoordinate::from_geodetic_deg(34.0522, -118.2437, 50.0);

// Straight-line (chord) distance (meters)
double chord_dist = coord1.distance_to(coord2);

// ECEF vector from coord1 to coord2
Vec3 vec = coord1.vector_to(coord2);

// Add ECEF offset
Vec3 offset{1000, 2000, 500};
auto coord3 = coord1 + offset;

// In-place modification
coord1 += offset;
```

### Geodesic Operations

```cpp
auto sf = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 0);
auto la = ECEFCoordinate::from_geodetic_deg(34.0522, -118.2437, 0);

// Great circle distance using Vincenty's formula (accurate for all distances)
double geodesic_dist = sf.geodesic_distance_to(la);  // ~559 km

// Bearing from North (radians, 0 = North, π/2 = East)
double bearing = sf.bearing_to(la);  // ~136° (Southeast)

// Move along geodesic by distance and bearing
auto target = sf.move_by_bearing(10000.0, bearing);  // 10km toward LA

// For long distances, use accurate Vincenty direct formula
auto target_accurate = sf.move_by_bearing_accurate(100000.0, bearing);  // 100km
```

### Local ENU Frames

ENU (East-North-Up) frames provide a local tangent plane coordinate system for intuitive relative positioning.

```cpp
// Create reference point
auto origin = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 10.0);

// Create local frame (one-time trig cost)
auto frame = origin.local_frame();

// Define offset in intuitive ENU coordinates
Vec3 enu_offset{100, 200, 50};  // 100m East, 200m North, 50m Up

// Convert to ECEF and create new coordinate
auto target = frame.enu_to_coord(enu_offset);

// Convert another coordinate to ENU relative to origin
auto other = ECEFCoordinate::from_geodetic_deg(37.7850, -122.4094, 15.0);
Vec3 enu_relative = frame.coord_to_enu(other);

std::cout << "Other is " << enu_relative.x << "m East, "
          << enu_relative.y << "m North, "
          << enu_relative.z << "m Up from origin\n";
```

### Local NED Frames

NED (North-East-Down) frames are common in aerospace and marine applications.

```cpp
// Create reference point
auto origin = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 10.0);

// Create NED frame
auto ned_frame = origin.local_ned_frame();

// Define offset in NED coordinates
Vec3 ned_offset{100, 200, -50};  // 100m North, 200m East, 50m Up (negative down)

// Convert to coordinate
auto target = ned_frame.ned_to_coord(ned_offset);

// Convert coordinate to NED
Vec3 ned_coords = ned_frame.coord_to_ned(target);
```

### Flat-Earth Dynamics (Flight Simulation)

**Critical for flight dynamics, trajectory planning, and navigation systems!**

Compute maneuvers in a local flat-earth frame, then map back to WGS84:

```cpp
// Starting position: San Francisco
auto ksfo = ECEFCoordinate::from_geodetic_deg(37.6213, -122.3790, 13.0);

// Fly 100 nautical miles Northeast at 30,000 ft, descend 10,000 ft
constexpr double NM_TO_M = 1852.0;
constexpr double FT_TO_M = 0.3048;

double distance_m = 100.0 * NM_TO_M;
double east_m = distance_m / std::sqrt(2.0);   // NE = 45°
double north_m = distance_m / std::sqrt(2.0);
double altitude_change = (30000.0 - 10000.0) * FT_TO_M - ksfo.altitude();

// Apply flat-earth displacement in ENU
Vec3 enu_displacement{east_m, north_m, altitude_change};
auto final_pos = ksfo.apply_enu_displacement(enu_displacement);

// Or use NED (common in aerospace)
Vec3 ned_displacement{north_m, east_m, -altitude_change};  // Negative = up
auto final_pos_ned = ksfo.apply_ned_displacement(ned_displacement);
```

**Multi-leg flight plans with frame updates:**

```cpp
// For long trajectories, update frame periodically for accuracy
std::vector<Vec3> flight_legs = {
    {150000, 0, 10000},      // 150km N, climb 10km
    {0, 100000, 0},          // 100km E, maintain altitude
    {50000, 50000, -5000}    // 50km N, 50km E, descend 5km
};

// Automatically updates local frame every 50km
auto final_pos = origin.apply_enu_trajectory(flight_legs, 50000.0);
```

**Why this matters:**
- Flight dynamics are computed in flat-earth (simple physics)
- Trajectory integration uses local coordinates
- Final position correctly mapped to curved Earth
- Altitude changes are preserved accurately
- Excellent accuracy for < 100 NM segments

### Geometric Utilities

**Surface normal, midpoints, and interpolation:**

```cpp
// Get surface normal (for gravity direction, ray tracing)
auto coord = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 0);
Vec3 normal = coord.surface_normal();  // Unit vector pointing up

// Compute midpoint between two locations
auto sf = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 0);
auto la = ECEFCoordinate::from_geodetic_deg(34.0522, -118.2437, 0);
auto mid = midpoint(sf, la);  // Geodesic midpoint

// Interpolate along geodesic (0 = start, 1 = end)
auto quarter = interpolate(sf, la, 0.25);  // 25% of the way
auto three_quarter = interpolate(sf, la, 0.75);  // 75% of the way

// Compute centroid of multiple points
std::vector<ECEFCoordinate> waypoints = {wp1, wp2, wp3, wp4};
auto center = centroid(waypoints);

// Commutative addition
Vec3 offset{1000, 2000, 500};
auto result1 = coord + offset;
auto result2 = offset + coord;  // Same result
```

### Performance Optimization Patterns

#### Pattern 1: Batch ECEF Operations

```cpp
// Create trajectory in ECEF space (no trig per point!)
auto start = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 10.0);
std::vector<ECEFCoordinate> trajectory;

for (int i = 0; i < 1000; ++i) {
    trajectory.push_back(start + Vec3{i * 10.0, 0, 0});  // Fast!
}

// Compute distances (no trig!)
for (size_t i = 1; i < trajectory.size(); ++i) {
    double dist = trajectory[i-1].distance_to(trajectory[i]);
}
```

#### Pattern 2: Reuse Local Frames

```cpp
auto origin = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 10.0);
auto frame = origin.local_frame();  // One-time trig cost

// Process many points relative to origin
std::vector<Vec3> enu_points;
for (const auto& coord : many_coordinates) {
    enu_points.push_back(frame.coord_to_enu(coord));  // Fast!
}
```

#### Pattern 3: Pre-cache Geodetic

```cpp
// If you know you'll need geodetic coordinates later
auto coord = ECEFCoordinate(4000000, 3000000, 2000000);
coord.precompute_geodetic();  // Compute now

// Later accesses are instant
double lat = coord.latitude();  // Cache hit
double lon = coord.longitude(); // Cache hit
```

## API Reference

### ECEFCoordinate

#### Constructors

```cpp
ECEFCoordinate(double x, double y, double z);
ECEFCoordinate(const Vec3& ecef);
ECEFCoordinate();  // Origin

static ECEFCoordinate from_geodetic(double lat_rad, double lon_rad, double alt_m);
static ECEFCoordinate from_geodetic_deg(double lat_deg, double lon_deg, double alt_m);
```

#### ECEF Access

```cpp
const Vec3& ecef() const;
double x() const;
double y() const;
double z() const;
```

#### Geodetic Access

```cpp
double latitude() const;        // Radians
double longitude() const;       // Radians
double altitude() const;        // Meters
double latitude_deg() const;    // Degrees
double longitude_deg() const;   // Degrees
```

#### Cached Trigonometry

```cpp
double sin_lat() const;
double cos_lat() const;
double sin_lon() const;
double cos_lon() const;
```

#### Vector Operations

```cpp
Vec3 vector_to(const ECEFCoordinate& other) const;
double distance_to(const ECEFCoordinate& other) const;  // Chord distance

ECEFCoordinate operator+(const Vec3& offset) const;
ECEFCoordinate operator-(const Vec3& offset) const;
ECEFCoordinate& operator+=(const Vec3& offset);
ECEFCoordinate& operator-=(const Vec3& offset);
```

#### Geodesic Operations

```cpp
double geodesic_distance_to(const ECEFCoordinate& other, int max_iterations = 20) const;
double bearing_to(const ECEFCoordinate& other) const;  // Radians from North
ECEFCoordinate move_by_bearing(double distance_m, double bearing_rad) const;
ECEFCoordinate move_by_bearing_accurate(double distance_m, double bearing_rad, int max_iterations = 20) const;
```

#### Local Frames

```cpp
ENUFrame local_frame() const;
NEDFrame local_ned_frame() const;
```

#### Flat-Earth Dynamics

```cpp
// Single displacement
ECEFCoordinate apply_enu_displacement(const Vec3& enu_displacement) const;
ECEFCoordinate apply_ned_displacement(const Vec3& ned_displacement) const;

// Multi-leg trajectory with frame updates
ECEFCoordinate apply_enu_trajectory(const std::vector<Vec3>& enu_displacements, 
                                    double update_interval_m = 50000.0) const;
ECEFCoordinate apply_ned_trajectory(const std::vector<Vec3>& ned_displacements,
                                    double update_interval_m = 50000.0) const;
```

#### Geometric Utilities

```cpp
Vec3 surface_normal() const;                              // Unit normal to ellipsoid
ECEFCoordinate midpoint(const ECEFCoordinate& other) const;  // Geodesic midpoint
```

#### Free Functions

```cpp
ECEFCoordinate operator+(const Vec3& offset, const ECEFCoordinate& coord);  // Commutative
ECEFCoordinate midpoint(const ECEFCoordinate& a, const ECEFCoordinate& b);
ECEFCoordinate interpolate(const ECEFCoordinate& a, const ECEFCoordinate& b, double t);
ECEFCoordinate centroid(const std::vector<ECEFCoordinate>& coords);
```

#### Utilities

```cpp
bool has_geodetic_cache() const;
void precompute_geodetic() const;
std::string to_string(int precision = 6, bool use_degrees = true) const;
std::string to_string_ecef(int precision = 3) const;
```

### ENUFrame

```cpp
struct ENUFrame {
    Vec3 east;   // Unit vector pointing East
    Vec3 north;  // Unit vector pointing North
    Vec3 up;     // Unit vector pointing Up
    Vec3 origin_ecef;  // Origin in ECEF
    
    Vec3 to_enu(const Vec3& ecef_vec) const;
    Vec3 to_ecef(const Vec3& enu_vec) const;
    Vec3 coord_to_enu(const ECEFCoordinate& coord) const;
    ECEFCoordinate enu_to_coord(const Vec3& enu_offset) const;
};
```

### NEDFrame

```cpp
struct NEDFrame {
    Vec3 north;  // Unit vector pointing North
    Vec3 east;   // Unit vector pointing East
    Vec3 down;   // Unit vector pointing Down
    Vec3 origin_ecef;  // Origin in ECEF
    
    Vec3 to_ned(const Vec3& ecef_vec) const;
    Vec3 to_ecef(const Vec3& ned_vec) const;
    Vec3 coord_to_ned(const ECEFCoordinate& coord) const;
    ECEFCoordinate ned_to_coord(const Vec3& ned_offset) const;
};
```

### WGS84 Constants

```cpp
struct WGS84 {
    static constexpr double a = 6378137.0;           // Semi-major axis (m)
    static constexpr double f = 1.0 / 298.257223563; // Flattening
    static constexpr double b = 6356752.314245;      // Semi-minor axis (m)
    static constexpr double e2 = 0.00669437999014;   // First eccentricity²
    static constexpr double ep2 = 0.00673949674228;  // Second eccentricity²
};
```

## Performance Characteristics

### Operation Costs

| Operation | Trig Calls | Complexity | Notes |
|-----------|-----------|------------|-------|
| ECEF construction | 0 | O(1) | Instant |
| Geodetic construction | 6 | O(1) | One-time cost, cached |
| ECEF access | 0 | O(1) | Direct memory access |
| Geodetic access (first) | 8-10 | O(1) | Iterative, then cached |
| Geodetic access (cached) | 0 | O(1) | Cache hit |
| Trig access (cached) | 0 | O(1) | Computed with geodetic |
| Chord distance | 0 | O(1) | Vec3 subtraction + norm |
| Geodesic distance | ~20-30 | O(1) | Vincenty, 1-3 iterations |
| Bearing calculation | 4 | O(1) | Uses local frame |
| Move by bearing | 6-8 | O(1) | Approximate, fast |
| Move by bearing (accurate) | ~20-30 | O(1) | Vincenty direct |
| Vector operations | 0 | O(1) | Pure Vec3 math |
| Local frame creation | 4 | O(1) | One-time, reuse frame |
| ENU/NED conversion | 0 | O(1) | Dot products only |

### Memory Usage

- **ECEFCoordinate**: 80 bytes
  - Vec3 (24 bytes)
  - 7 × std::optional<double> (56 bytes)
- **ENUFrame**: 104 bytes
  - 3 × Vec3 (72 bytes)
  - ECEFCoordinate (80 bytes, but shared)

### Optimization Tips

1. **Work in ECEF when possible**: Avoid geodetic conversions in tight loops
2. **Reuse local frames**: Create once, use many times
3. **Pre-cache if needed**: Call `precompute_geodetic()` before hot loops
4. **Batch operations**: Process multiple coordinates in ECEF space
5. **Cache trig values**: Use `sin_lat()`, `cos_lat()` instead of recomputing

## Algorithm Details

### Geodetic to ECEF

Uses standard WGS84 transformation:

```
N = a / sqrt(1 - e² sin²(lat))
X = (N + h) cos(lat) cos(lon)
Y = (N + h) cos(lat) sin(lon)
Z = (N(1 - e²) + h) sin(lat)
```

### ECEF to Geodetic

Uses Bowring's iterative method:
- Converges in 2-3 iterations for Earth-surface points
- Handles polar singularities
- Accuracy: < 1mm for altitudes < 10,000 km

## Thread Safety

- **Read-only operations**: Thread-safe
- **Lazy computation**: Thread-safe (uses mutable cache)
- **Modifications**: Not thread-safe (use external synchronization)

## Dependencies

- **dstl::LinAlg**: Vec3 vector operations
- **C++17**: std::optional for caching
- **Standard library**: <cmath> for trigonometric functions

## Examples

### Example 1: GPS Trajectory Processing

```cpp
#include <dstl/Geodesy.h>
#include <vector>

using namespace dstl::geo;

// Process GPS trajectory
std::vector<ECEFCoordinate> load_gps_trajectory() {
    std::vector<ECEFCoordinate> trajectory;
    
    // Load from GPS data (lat/lon/alt)
    trajectory.push_back(ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 10));
    trajectory.push_back(ECEFCoordinate::from_geodetic_deg(37.7750, -122.4193, 12));
    // ... more points
    
    return trajectory;
}

void analyze_trajectory(const std::vector<ECEFCoordinate>& traj) {
    // Compute total distance (fast - no trig!)
    double total_distance = 0;
    for (size_t i = 1; i < traj.size(); ++i) {
        total_distance += traj[i-1].distance_to(traj[i]);
    }
    
    std::cout << "Total distance: " << total_distance << " meters\n";
}
```

### Example 2: Local Area Mapping

```cpp
#include <dstl/Geodesy.h>

using namespace dstl::geo;
using namespace dstl::linalg;

void map_local_area() {
    // Reference point
    auto origin = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 0);
    auto frame = origin.local_frame();
    
    // Define points in intuitive ENU coordinates
    std::vector<Vec3> landmarks = {
        {100, 0, 0},      // 100m East
        {0, 200, 0},      // 200m North
        {-50, 150, 10},   // 50m West, 150m North, 10m Up
    };
    
    // Convert to ECEF coordinates
    std::vector<ECEFCoordinate> ecef_landmarks;
    for (const auto& enu : landmarks) {
        ecef_landmarks.push_back(frame.enu_to_coord(enu));
    }
    
    // Now can compute distances, etc. in ECEF
    for (const auto& landmark : ecef_landmarks) {
        double dist = origin.distance_to(landmark);
        std::cout << "Distance: " << dist << " m\n";
    }
}
```

### Example 3: Satellite Ground Track

```cpp
#include <dstl/Geodesy.h>

using namespace dstl::geo;
using namespace dstl::linalg;

void compute_ground_track() {
    // Satellite position in ECEF
    ECEFCoordinate sat_pos(7000000, 0, 0);  // ~600 km altitude
    
    // Ground station
    auto station = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 0);
    
    // Compute slant range (fast!)
    double range = station.distance_to(sat_pos);
    Vec3 los_vector = station.vector_to(sat_pos);
    
    std::cout << "Slant range: " << range / 1000.0 << " km\n";
    
    // Get satellite geodetic position (lazy trig)
    std::cout << "Satellite lat: " << sat_pos.latitude_deg() << "°\n";
    std::cout << "Satellite lon: " << sat_pos.longitude_deg() << "°\n";
    std::cout << "Satellite alt: " << sat_pos.altitude() / 1000.0 << " km\n";
}
```

## See Also

- [LinAlg Component](../LinAlg/README.md) - Vector and matrix operations
- [DSTL Documentation Index](../DOCUMENTATION_INDEX.md) - Complete library overview

## References

- WGS84 Specification: NIMA TR8350.2
- Bowring, B.R. (1976). "Transformation from spatial to geographical coordinates"
- Zhu, J. (1994). "Conversion of Earth-centered Earth-fixed coordinates to geodetic coordinates"
