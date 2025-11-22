# Geodesy Quick Reference

## Construction

```cpp
// From ECEF (fast, no trig)
ECEFCoordinate coord(x, y, z);
ECEFCoordinate coord(vec3);

// From geodetic (radians)
auto coord = ECEFCoordinate::from_geodetic(lat_rad, lon_rad, alt_m);

// From geodetic (degrees) - most common
auto coord = ECEFCoordinate::from_geodetic_deg(lat_deg, lon_deg, alt_m);
```

## Access

```cpp
// ECEF (fast, no trig)
const Vec3& ecef = coord.ecef();
double x = coord.x();
double y = coord.y();
double z = coord.z();

// Geodetic (lazy, cached)
double lat_rad = coord.latitude();
double lon_rad = coord.longitude();
double alt_m = coord.altitude();
double lat_deg = coord.latitude_deg();
double lon_deg = coord.longitude_deg();

// Cached trig (free after geodetic access)
double sl = coord.sin_lat();
double cl = coord.cos_lat();
double sn = coord.sin_lon();
double cn = coord.cos_lon();
```

## Distance

```cpp
// Chord distance (fast, no trig)
double chord = coord1.distance_to(coord2);

// Geodesic distance (accurate, ~500ns)
double geodesic = coord1.geodesic_distance_to(coord2);
```

## Navigation

```cpp
// Bearing from North (radians, 0-2π)
double bearing = coord1.bearing_to(coord2);

// Move by distance and bearing
auto target = origin.move_by_bearing(10000.0, bearing);  // Fast, < 100km
auto target = origin.move_by_bearing_accurate(100000.0, bearing);  // Accurate, all distances
```

## Vector Operations

```cpp
// ECEF vector (fast, no trig)
Vec3 vec = coord1.vector_to(coord2);

// Add/subtract ECEF offsets
auto coord3 = coord1 + Vec3{dx, dy, dz};
auto coord4 = coord1 - Vec3{dx, dy, dz};
coord1 += Vec3{dx, dy, dz};
coord1 -= Vec3{dx, dy, dz};
```

## Local Frames

### ENU (East-North-Up)
```cpp
auto frame = origin.local_frame();

// Convert to/from ENU
Vec3 enu = frame.coord_to_enu(target);
auto target = frame.enu_to_coord(Vec3{east, north, up});

// Convert vectors
Vec3 enu_vec = frame.to_enu(ecef_vec);
Vec3 ecef_vec = frame.to_ecef(enu_vec);
```

### NED (North-East-Down)
```cpp
auto frame = origin.local_ned_frame();

// Convert to/from NED
Vec3 ned = frame.coord_to_ned(target);
auto target = frame.ned_to_coord(Vec3{north, east, down});

// Convert vectors
Vec3 ned_vec = frame.to_ned(ecef_vec);
Vec3 ecef_vec = frame.to_ecef(ned_vec);
```

## Comparison

```cpp
// Equality (1 micron tolerance)
bool same = coord1 == coord2;
bool different = coord1 != coord2;

// Custom tolerance
bool close = coord1.approx_equal(coord2, 1.0);  // Within 1m
```

## Utilities

```cpp
// String conversion
std::string str = coord.to_string();  // Geodetic
std::string str = coord.to_string_ecef();  // ECEF

// Stream output
std::cout << coord << "\n";

// Cache management
bool cached = coord.has_geodetic_cache();
coord.precompute_geodetic();  // Force computation

// Container support
std::unordered_map<ECEFCoordinate, std::string> map;
std::unordered_set<ECEFCoordinate> set;
```

## Constants

```cpp
dstl::geo::PI                    // π
dstl::geo::DEG_TO_RAD           // π/180
dstl::geo::RAD_TO_DEG           // 180/π
dstl::geo::WGS84::a             // 6378137.0 m (semi-major axis)
dstl::geo::WGS84::b             // 6356752.3 m (semi-minor axis)
dstl::geo::WGS84::f             // 1/298.257223563 (flattening)
dstl::geo::WGS84::e2            // First eccentricity²
```

## Performance Guide

| Operation | Time | Trig | When to Use |
|-----------|------|------|-------------|
| `distance_to()` | ~11 ns | 0 | Short ranges, performance critical |
| `geodesic_distance_to()` | ~500 ns | 20-30 | Accurate surface distance |
| `bearing_to()` | ~100 ns | 4 | Navigation, direction finding |
| `move_by_bearing()` | ~150 ns | 6-8 | Fast movement, < 100 km |
| `move_by_bearing_accurate()` | ~500 ns | 20-30 | Accurate movement, all distances |
| `local_frame()` | ~100 ns | 4 | Create once, reuse many times |

## Common Patterns

### Pattern 1: Batch ECEF Operations
```cpp
auto start = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 10.0);
for (int i = 0; i < 1000; ++i) {
    auto point = start + Vec3{i * 10.0, 0, 0};  // Fast!
    double dist = start.distance_to(point);     // Fast!
}
```

### Pattern 2: Reuse Local Frames
```cpp
auto origin = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 10.0);
auto frame = origin.local_frame();  // One-time cost

for (const auto& coord : many_coords) {
    Vec3 enu = frame.coord_to_enu(coord);  // Fast!
}
```

### Pattern 3: Navigation Loop
```cpp
auto current = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0);
auto waypoint = ECEFCoordinate::from_geodetic_deg(38.0, -121.0, 0);

while (current.geodesic_distance_to(waypoint) > 10.0) {
    double bearing = current.bearing_to(waypoint);
    current = current.move_by_bearing(1.0, bearing);  // Move 1m
}
```

### Pattern 4: Pre-cache for Hot Loops
```cpp
auto coord = ECEFCoordinate(x, y, z);
coord.precompute_geodetic();  // Compute now

for (int i = 0; i < 1000000; ++i) {
    double lat = coord.latitude();  // Cache hit, instant!
    double sl = coord.sin_lat();    // Cache hit, instant!
}
```

## Accuracy Guide

| Operation | Accuracy | Range |
|-----------|----------|-------|
| ECEF ↔ Geodetic | < 1 mm | All altitudes < 10,000 km |
| Chord distance | Exact | All distances |
| Geodesic distance | < 1 mm | All Earth-surface distances |
| Bearing | < 0.01° | All distances |
| Movement (fast) | ~1 m | < 100 km |
| Movement (accurate) | < 1 mm | All Earth-surface distances |

## Edge Cases

### Poles
```cpp
auto north_pole = ECEFCoordinate::from_geodetic_deg(90.0, 0.0, 0.0);
auto south_pole = ECEFCoordinate::from_geodetic_deg(-90.0, 0.0, 0.0);
// Handled correctly
```

### Date Line
```cpp
auto west = ECEFCoordinate::from_geodetic_deg(0.0, 179.9, 0.0);
auto east = ECEFCoordinate::from_geodetic_deg(0.0, -179.9, 0.0);
double dist = west.geodesic_distance_to(east);  // Correct
```

### Antipodal Points
```cpp
auto coord1 = ECEFCoordinate::from_geodetic_deg(0.0, 0.0, 0.0);
auto coord2 = ECEFCoordinate::from_geodetic_deg(0.0, 180.0, 0.0);
double dist = coord1.geodesic_distance_to(coord2);  // ~20,000 km (approximate)
```

## Troubleshooting

### "Distance seems wrong"
- Use `geodesic_distance_to()` for surface distance, not `distance_to()`
- `distance_to()` returns chord (straight-line) distance

### "Bearing is negative"
- `bearing_to()` returns [0, 2π), never negative
- 0 = North, π/2 = East, π = South, 3π/2 = West

### "Movement doesn't preserve altitude"
- `move_by_bearing()` and `move_by_bearing_accurate()` maintain origin altitude
- Compute altitude changes separately if needed

### "Performance is slow"
- Pre-cache geodetic: `coord.precompute_geodetic()`
- Reuse local frames instead of recreating
- Use `distance_to()` instead of `geodesic_distance_to()` when appropriate

## Examples

### GPS Waypoint Navigation
```cpp
auto current = ECEFCoordinate::from_geodetic_deg(gps_lat, gps_lon, gps_alt);
auto waypoint = ECEFCoordinate::from_geodetic_deg(wp_lat, wp_lon, wp_alt);

double distance = current.geodesic_distance_to(waypoint);
double bearing = current.bearing_to(waypoint);

std::cout << "Distance: " << distance << " m\n";
std::cout << "Bearing: " << bearing * RAD_TO_DEG << "°\n";
```

### Local Mapping
```cpp
auto origin = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 0);
auto frame = origin.local_frame();

// Define landmarks in ENU
std::vector<Vec3> landmarks = {
    {100, 0, 0},    // 100m East
    {0, 200, 0},    // 200m North
    {-50, 150, 10}  // 50m West, 150m North, 10m Up
};

for (const auto& enu : landmarks) {
    auto coord = frame.enu_to_coord(enu);
    std::cout << coord << "\n";
}
```

### Dead Reckoning
```cpp
auto position = ECEFCoordinate::from_geodetic_deg(37.0, -122.0, 0);
double heading = 45.0 * DEG_TO_RAD;  // Northeast
double speed = 10.0;  // m/s
double dt = 1.0;      // 1 second

// Update position
position = position.move_by_bearing(speed * dt, heading);
```

---

**Quick Tip**: For best performance, work in ECEF space and only convert to geodetic when displaying to users!
