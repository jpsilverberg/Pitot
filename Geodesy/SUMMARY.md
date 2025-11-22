# Geodesy Module Summary

## Overview

The Geodesy module is a new DSTL component providing high-performance WGS84 coordinate transformations with a focus on minimizing expensive trigonometric operations.

## Key Design Decisions

### 1. ECEF as Primary Storage
- **Decision**: Store coordinates in ECEF (X, Y, Z) format
- **Benefit**: Enables zero-cost vector operations using `dstl::linalg::Vec3`
- **Trade-off**: Geodetic access requires computation (mitigated by caching)

### 2. Cached Trigonometric Values
- **Decision**: Cache sin(lat), cos(lat), sin(lon), cos(lon) alongside geodetic coordinates
- **Benefit**: Eliminates repeated expensive trig calculations
- **Memory Cost**: 32 bytes (4 doubles) - negligible for the performance gain
- **Use Case**: Essential for local frame creation and repeated geodetic operations

### 3. Lazy Geodetic Computation
- **Decision**: Compute geodetic coordinates only when accessed
- **Benefit**: Batch ECEF operations avoid all trigonometry
- **Implementation**: Uses `std::optional` for cache with mutable members

## Performance Characteristics

| Operation | Trig Calls | Complexity | Notes |
|-----------|-----------|------------|-------|
| ECEF construction | 0 | O(1) | Instant |
| Geodetic construction | 6 | O(1) | One-time, cached |
| Distance calculation | 0 | O(1) | Vec3 math only |
| Geodetic access (first) | 8-10 | O(1) | Iterative, then cached |
| Geodetic access (cached) | 0 | O(1) | Cache hit |
| Trig access (cached) | 0 | O(1) | Computed with geodetic |
| Local frame creation | 4 | O(1) | One-time, reuse frame |

## Memory Layout

```
ECEFCoordinate: 80 bytes total
├─ Vec3 m_ecef: 24 bytes (3 × double)
├─ optional<double> m_lat: 8 bytes
├─ optional<double> m_lon: 8 bytes
├─ optional<double> m_alt: 8 bytes
├─ optional<double> m_sin_lat: 8 bytes
├─ optional<double> m_cos_lat: 8 bytes
├─ optional<double> m_sin_lon: 8 bytes
└─ optional<double> m_cos_lon: 8 bytes
```

## Usage Patterns

### Pattern 1: Fast ECEF Operations
```cpp
// Create from geodetic (one-time trig cost)
auto start = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 10.0);

// Fast operations (no trig!)
for (int i = 0; i < 1000; ++i) {
    auto point = start + Vec3{i * 10.0, 0, 0};
    double dist = start.distance_to(point);  // Pure Vec3 math
}
```

### Pattern 2: Reuse Local Frames
```cpp
auto origin = ECEFCoordinate::from_geodetic_deg(37.7749, -122.4194, 10.0);
auto frame = origin.local_frame();  // One-time trig cost

// Process many points (no trig per point!)
for (const auto& coord : many_coordinates) {
    Vec3 enu = frame.coord_to_enu(coord);  // Fast!
}
```

### Pattern 3: Pre-cache When Needed
```cpp
auto coord = ECEFCoordinate(4000000, 3000000, 2000000);
coord.precompute_geodetic();  // Compute now

// Later accesses are instant
double lat = coord.latitude();  // Cache hit
double sin_lat = coord.sin_lat();  // Cache hit
```

## Why Cache Trigonometry?

### Scenario: Creating Local Frames

Without caching:
```cpp
auto frame = coord.local_frame();
// Computes: sin(lat), cos(lat), sin(lon), cos(lon) = 4 trig calls
```

With caching:
```cpp
auto frame = coord.local_frame();
// Uses cached values from geodetic computation = 0 additional trig calls
```

### Scenario: Repeated Geodetic Access

Without caching:
```cpp
for (int i = 0; i < 100; ++i) {
    double lat = coord.latitude();  // Recomputes every time
    double sin_lat = std::sin(lat);  // Additional trig call
}
// Total: 100 × (8-10 + 1) = 900-1100 trig calls
```

With caching:
```cpp
for (int i = 0; i < 100; ++i) {
    double lat = coord.latitude();  // Cache hit after first
    double sin_lat = coord.sin_lat();  // Cache hit
}
// Total: 8-10 trig calls (first iteration only)
```

## Comparison with Alternative Designs

### Alternative 1: Geodetic Primary Storage
```cpp
struct GeodeticCoordinate {
    double lat, lon, alt;  // 24 bytes
};
```
- **Pros**: Smaller (24 vs 80 bytes), human-readable
- **Cons**: Every distance calculation requires trig, slow vector operations
- **Verdict**: Poor for computational geometry

### Alternative 2: No Caching
```cpp
struct ECEFCoordinate {
    Vec3 ecef;  // 24 bytes
    // Compute geodetic on every access
};
```
- **Pros**: Smaller (24 vs 80 bytes)
- **Cons**: Repeated trig calculations, no trig value reuse
- **Verdict**: Slower for typical use cases

### Alternative 3: Cache Geodetic Only
```cpp
struct ECEFCoordinate {
    Vec3 ecef;  // 24 bytes
    optional<double> lat, lon, alt;  // 24 bytes
    // Total: 48 bytes
};
```
- **Pros**: Smaller (48 vs 80 bytes)
- **Cons**: Must recompute trig for local frames
- **Verdict**: Misses key optimization opportunity

### Our Design: Cache Everything
```cpp
struct ECEFCoordinate {
    Vec3 ecef;  // 24 bytes
    optional<double> lat, lon, alt;  // 24 bytes
    optional<double> sin_lat, cos_lat, sin_lon, cos_lon;  // 32 bytes
    // Total: 80 bytes
};
```
- **Pros**: Zero trig for cached access, optimal for local frames
- **Cons**: Larger memory footprint
- **Verdict**: Best for performance-critical applications

## Benchmarks (Estimated)

Based on typical CPU performance (Intel/AMD x86-64):

| Operation | Cycles | Time @ 3GHz |
|-----------|--------|-------------|
| Vec3 addition | ~3 | 1 ns |
| Vec3 dot product | ~10 | 3 ns |
| Vec3 norm (sqrt) | ~30 | 10 ns |
| sin/cos | ~100-200 | 30-60 ns |
| atan2 | ~150-300 | 50-100 ns |

**Distance calculation (ECEF)**:
- Vec3 subtraction + norm = ~33 cycles = ~11 ns

**Distance calculation (geodetic, no cache)**:
- Compute lat/lon (8-10 trig) + Vincenty formula (4 trig) = ~1500 cycles = ~500 ns

**Speedup**: 45x faster with ECEF!

## Future Enhancements

See [IMPROVEMENTS.md](IMPROVEMENTS.md) for detailed improvement plans:

1. **Great circle distance** - Vincenty/Karney algorithms
2. **Velocity/acceleration** - Support for dynamic systems
3. **UTM coordinates** - Universal Transverse Mercator
4. **Multiple ellipsoids** - GRS80, NAD83, etc.
5. **Serialization** - JSON, binary formats

## Conclusion

The Geodesy module's design prioritizes performance through:
1. ECEF primary storage for fast vector operations
2. Lazy geodetic computation to avoid unnecessary work
3. Aggressive caching of both geodetic and trigonometric values

The cached trigonometry design is particularly noteworthy - it adds only 32 bytes but eliminates repeated expensive calculations, making it ideal for applications that need local frames or repeated geodetic access. This is a pattern worth considering for other components dealing with expensive computations.
