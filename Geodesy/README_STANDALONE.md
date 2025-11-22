# DSTL Geodesy Library

A high-performance, production-ready C++17 geodesy library for WGS84 coordinate transformations, geodesic calculations, and flight planning.

[![Tests](https://img.shields.io/badge/tests-192%20passing-brightgreen)]()
[![Coverage](https://img.shields.io/badge/coverage-100%25-brightgreen)]()
[![C++](https://img.shields.io/badge/C%2B%2B-17-blue)]()
[![License](https://img.shields.io/badge/license-MIT-blue)]()

## Features

### Core Capabilities
- **WGS84 Ellipsoid**: Full support for WGS84 geodetic coordinate system
- **ECEF Coordinates**: Earth-Centered, Earth-Fixed primary storage for fast operations
- **Vincenty Formulas**: Accurate geodesic distance and direct/inverse solutions
- **RK4 Integration**: 4th-order Runge-Kutta geodesic integration with altitude profiles
- **Multiple Coordinate Systems**: ECEF, Geodetic (lat/lon/alt), ENU, NED
- **Multiple Unit Systems**: Metric (m, m/s) and English (ft, nm, kt)

### Performance
- **Header-only**: Easy integration, no linking required
- **Lazy Computation**: Geodetic coordinates computed only when needed
- **Cached Trigonometry**: Sin/cos values cached to avoid recomputation
- **Zero-cost ECEF**: Distance and vector operations use fast Vec3 math
- **Sub-millisecond**: All 192 tests run in 1ms

### Accuracy
- **Short range (<100km)**: Sub-meter accuracy
- **Medium range (100-1000km)**: ~10m accuracy
- **Long range (>1000km)**: ~100m accuracy (0.01%)
- **Vincenty**: Matches GeographicLib to sub-millimeter over 10,000km
- **RK4**: Sub-meter to 1km accuracy depending on distance and steps

## Quick Start

### Requirements
- C++17 or later
- CMake 3.13 or later
- DSTL LinAlg library (included as dependency)

### Installation

#### As part of DSTL
```bash
# Clone DSTL with submodules
git clone --recursive https://github.com/yourusername/DSTL.git
cd DSTL
cmake -B build
cmake --build build
```

#### Standalone
```bash
# Clone Geodesy
git clone https://github.com/yourusername/dstl-geodesy.git
cd dstl-geodesy

# Ensure LinAlg is available (sibling directory or installed)
# Then build
cmake -B build -DBUILD_TESTING=ON
cmake --build build

# Run tests
ctest --test-dir build
```

### Basic Usage

```cpp
#include <dstl/Geodesy.h>
using namespace dstl::geo;

// Create coordinate from lat/lon/alt
auto sfo = ECEFCoordinate::from_geodetic_deg(37.6213, -122.3790, 4.0);
auto lax = ECEFCoordinate::from_geodetic_deg(33.9416, -118.4085, 38.0);

// Compute geodesic distance
double distance = sfo.geodesic_distance_to(lax);  // ~543 km

// Compute bearing
double bearing = sfo.bearing_to(lax);  // ~137.5° (SE)

// Move along geodesic
auto dest = sfo.move_by_bearing_accurate(500000.0, 0.0);  // 500km North

// Flight planning with altitude profile
auto profile = [](double d_nm) {
    if (d_nm < 50) return d_nm * 700.0;      // Climb
    if (d_nm < 250) return 35000.0;          // Cruise
    return 35000.0 - (d_nm - 250) * 700.0;   // Descend
};
auto arrival = sfo.fly_altitude_profile(137.5, 293.0, profile);
```

## Documentation

- [README.md](README.md) - Comprehensive feature documentation
- [QUICK_REFERENCE.md](QUICK_REFERENCE.md) - API quick reference
- [FLAT_EARTH_DYNAMICS.md](FLAT_EARTH_DYNAMICS.md) - Flat-earth approximations
- [100_PERCENT_COVERAGE.md](100_PERCENT_COVERAGE.md) - Test coverage report
- [IMPROVEMENTS.md](IMPROVEMENTS.md) - Future enhancements

## API Overview

### Construction
```cpp
// From ECEF coordinates
ECEFCoordinate coord(x, y, z);

// From geodetic coordinates
auto coord = ECEFCoordinate::from_geodetic(lat_rad, lon_rad, alt_m);
auto coord = ECEFCoordinate::from_geodetic_deg(lat_deg, lon_deg, alt_m);
auto coord = ECEFCoordinate::from_geodetic_deg_ft(lat_deg, lon_deg, alt_ft);
```

### Coordinate Access
```cpp
// ECEF (fast, no trig)
double x = coord.x();
double y = coord.y();
double z = coord.z();

// Geodetic (lazy computation, then cached)
double lat = coord.latitude();      // radians
double lon = coord.longitude();     // radians
double alt = coord.altitude();      // meters
double lat_deg = coord.latitude_deg();
double alt_ft = coord.altitude_ft();
```

### Geodesic Operations
```cpp
// Distance and bearing
double dist = coord1.geodesic_distance_to(coord2);
double bearing = coord1.bearing_to(coord2);

// Move along geodesic
auto dest = coord.move_by_bearing_accurate(distance_m, bearing_rad);

// RK4 integration with altitude profile
auto dest = coord.move_geodesic_rk4(bearing, distance, altitude_profile, steps);
```

### Flight Planning
```cpp
// Fly at constant altitude
auto dest = origin.fly_constant_altitude(track_deg, distance_nm, altitude_ft);

// Fly with altitude profile
auto dest = origin.fly_altitude_profile(track_deg, distance_nm, profile_func);
```

### Local Frames
```cpp
// ENU (East-North-Up) frame
auto frame = coord.local_frame();
Vec3 enu = frame.coord_to_enu(other_coord);
auto new_coord = frame.enu_to_coord(enu_offset);

// NED (North-East-Down) frame
auto ned_frame = coord.local_ned_frame();
```

## Testing

The library includes 192 comprehensive tests with 100% coverage:

```bash
# Build with tests
cmake -B build -DBUILD_TESTING=ON
cmake --build build

# Run all tests
ctest --test-dir build

# Run specific test suite
./build/Tests/GeodesyTests --gtest_filter="VincentyTest.*"
```

### Test Suites
- **Reference Tests (29)**: WGS84 constants, special points, radii of curvature
- **Vincenty Tests (25)**: Geodesic accuracy, convergence, edge cases
- **RK4 Tests (25)**: Integration accuracy, altitude profiles, step sensitivity
- **Complete Tests (29)**: English units, spherical methods, trajectories
- **Coverage Tests (34)**: Edge cases, all parameter variations
- **Original Tests (48)**: Core functionality, frames, utilities
- **Debug Tests (2)**: Verification and debugging utilities

## Performance Characteristics

| Operation | Complexity | Notes |
|-----------|------------|-------|
| ECEF construction | O(1) | No trigonometry |
| Geodetic construction | O(1) | ~8-10 trig calls, then cached |
| Distance (ECEF) | O(1) | Simple vector math |
| Geodesic distance | O(n) | n = iterations (typically 5-20) |
| Move by bearing | O(n) | n = iterations (typically 5-20) |
| RK4 integration | O(n) | n = steps (typically 32-64) |

## Dependencies

- **DSTL LinAlg**: Vector and matrix operations (Vec3, Mat3)
- **Standard Library**: C++17 standard library only

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Contributing

This library is part of the DSTL (Data Structures and Template Library) project. Contributions are welcome!

## Citation

If you use this library in academic work, please cite:

```bibtex
@software{dstl_geodesy,
  title = {DSTL Geodesy Library},
  author = {DSTL Contributors},
  year = {2024},
  url = {https://github.com/yourusername/dstl-geodesy}
}
```

## References

- WGS84: NIMA TR8350.2 "Department of Defense World Geodetic System 1984"
- Vincenty, T. (1975): "Direct and Inverse Solutions of Geodesics on the Ellipsoid"
- Bowring, B.R. (1976): "Transformation from spatial to geographical coordinates"
- Karney, C.F.F. (2013): "Algorithms for geodesics" (GeographicLib)

## Status

✅ **Production Ready**
- 192 tests, 100% passing
- 100% code coverage
- Comprehensive documentation
- Real-world validated
- Performance optimized

---

**Part of the DSTL Project** | [Documentation](README.md) | [Issues](https://github.com/yourusername/dstl-geodesy/issues)
