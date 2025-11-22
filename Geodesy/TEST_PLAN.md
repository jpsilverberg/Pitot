# Geodesy Comprehensive Test Plan - 100% Coverage

## Objective
Achieve 100% test coverage of all Geodesy.h functions with mathematically known reference points and edge cases.

## Current Coverage Analysis

### ✅ Already Tested (from test_geodesy.cpp)
- Basic construction (ECEF, geodetic, degrees)
- Round-trip conversions
- Caching mechanisms
- Vector operations (distance_to, vector_to, offsets)
- ENU/NED frame creation and conversions
- Special cases (equator, poles)
- Geodesic distance (short, long, antipodal)
- Bearing calculations (N, E, S, W)
- Basic movement (move_by_bearing)
- Flat-earth displacements (ENU/NED)
- Trajectories
- Geometric utilities (surface_normal, midpoint, interpolate, centroid)
- String conversions
- Equality and hashing
- Container usage (unordered_map)

### ❌ Missing Coverage - Functions Not Yet Tested

#### 1. **Vincenty Direct Formula**
- `move_by_bearing_accurate()` - only basic test exists
- Need: Long-distance accuracy tests with known geodetic solutions
- Need: Convergence behavior tests
- Need: Edge cases (near-poles, antipodal)

#### 2. **RK4 Geodesic Integration**
- `move_geodesic_rk4()` - NOT TESTED
- `move_geodesic_rk4_constant_altitude()` - NOT TESTED
- `fly_constant_altitude()` - NOT TESTED
- `fly_altitude_profile()` - NOT TESTED
- Need: Accuracy comparison vs Vincenty
- Need: Altitude profile integration tests
- Need: Step count sensitivity analysis

#### 3. **Spherical Displacement Methods**
- `apply_enu_displacement_spherical()` - NOT TESTED
- `apply_enu_displacement_spherical_ft()` - NOT TESTED
- Need: Accuracy vs flat-earth at various distances
- Need: Iterative refinement convergence tests
- Need: Altitude handling verification

#### 4. **English Units Conversions**
- `from_geodetic_ft()` - NOT TESTED
- `altitude_ft()` - NOT TESTED
- `apply_enu_displacement_ft()` - NOT TESTED
- `apply_ned_displacement_ft()` - NOT TESTED
- Need: Conversion accuracy tests
- Need: Round-trip ft ↔ m tests

#### 5. **Trajectory Frame Updates**
- `apply_enu_trajectory()` with `update_interval_m` - NOT TESTED
- `apply_ned_trajectory()` with `update_interval_m` - NOT TESTED
- Need: Frame update frequency tests
- Need: Accuracy vs single-frame application

#### 6. **Radii of Curvature**
- `meridian_radius()` - NOT TESTED
- `prime_vertical_radius()` - NOT TESTED
- Need: Known values at equator, poles, 45°
- Need: Comparison with WGS84 formulas

#### 7. **Edge Cases & Stress Tests**
- Longitude wrapping (±180°)
- Very high altitudes (LEO orbit ~400km)
- Very low altitudes (below sea level)
- Numerical stability near poles (89.9°, 89.99°, 89.999°)
- Zero-distance operations
- Extremely long distances (>10,000 km)
- Rapid altitude changes in trajectories

#### 8. **Performance Characteristics**
- Cache hit rates
- ECEF operations remain trig-free
- Lazy computation verification

## Test Data Sources

### 1. **GeographicLib Test Data**
Use Charles Karney's GeographicLib test suite for Vincenty validation:
- https://geographiclib.sourceforge.io/html/geodesic.html
- Known accurate solutions for geodesic problems

### 2. **NIMA/NGA Test Points**
Standard geodetic test points from National Geospatial-Intelligence Agency:
- Equator: (0°, 0°, 0m)
- North Pole: (90°, 0°, 0m)
- South Pole: (-90°, 0°, 0m)
- Prime Meridian: (51.4778°, 0°, 0m) - Greenwich
- International Date Line: (0°, 180°, 0m)

### 3. **Known Airport Pairs**
Real-world navigation test cases with published distances:
- JFK (New York) to LHR (London): ~5,585 km
- LAX (Los Angeles) to SYD (Sydney): ~12,051 km
- DXB (Dubai) to SFO (San Francisco): ~13,006 km
- GRU (São Paulo) to NRT (Tokyo): ~18,374 km

### 4. **WGS84 Reference Values**
- Equatorial radius: 6,378,137.0 m
- Polar radius: 6,356,752.314245 m
- Meridian radius at equator: 6,335,439.327 m
- Prime vertical radius at equator: 6,378,137.0 m
- Meridian radius at poles: 6,399,593.626 m

### 5. **Vincenty Test Cases**
From Vincenty's original 1975 paper:
- Antipodal points
- Near-antipodal points
- Equatorial paths
- Meridional paths

## Test Organization

### Phase 1: Reference Point Validation (test_geodesy_reference.cpp)
- WGS84 constants verification
- Known ECEF ↔ Geodetic conversions
- Radii of curvature at special latitudes
- Airport coordinate validation

### Phase 2: Vincenty Accuracy Tests (test_geodesy_vincenty.cpp)
- Short distance (<100 km)
- Medium distance (100-1000 km)
- Long distance (1000-10000 km)
- Very long distance (>10000 km)
- Antipodal and near-antipodal
- Convergence behavior
- Comparison with GeographicLib

### Phase 3: RK4 Integration Tests (test_geodesy_rk4.cpp)
- Constant altitude paths
- Linear altitude profiles
- Nonlinear altitude profiles (climb/descent)
- Step count sensitivity
- Accuracy vs Vincenty direct
- Long-range integration (>1000 km)

### Phase 4: Spherical Methods Tests (test_geodesy_spherical.cpp)
- Iterative refinement convergence
- Accuracy vs flat-earth
- Accuracy vs ellipsoidal methods
- Altitude handling
- Distance range validation

### Phase 5: Units & Conversions Tests (test_geodesy_units.cpp)
- Feet ↔ meters conversions
- Nautical miles ↔ meters
- Knots ↔ m/s
- Round-trip accuracy
- English units API coverage

### Phase 6: Edge Cases & Stress Tests (test_geodesy_edge_cases.cpp)
- Polar regions (>89°)
- Longitude wrapping
- High altitudes (LEO)
- Below sea level
- Zero distances
- Numerical stability
- Extreme trajectories

### Phase 7: Performance Tests (test_geodesy_performance.cpp)
- Cache behavior verification
- ECEF operation trig-free guarantee
- Lazy computation validation
- Benchmark comparisons

## Success Criteria

1. ✅ 100% function coverage
2. ✅ All tests pass with strict tolerances
3. ✅ Known reference points validated to <1m accuracy
4. ✅ Vincenty matches GeographicLib to <1mm over 10,000 km
5. ✅ RK4 matches Vincenty to <1m over 10,000 km
6. ✅ All edge cases handled gracefully
7. ✅ No warnings or errors in test compilation
8. ✅ Tests run in <5 seconds total

## Implementation Order

1. **First**: Create reference data file with known coordinates
2. **Second**: Implement Phase 1 (reference validation)
3. **Third**: Implement Phase 2 (Vincenty accuracy)
4. **Fourth**: Implement Phase 3 (RK4 integration)
5. **Fifth**: Implement Phases 4-6 (remaining coverage)
6. **Sixth**: Implement Phase 7 (performance validation)
7. **Finally**: Generate coverage report and fill gaps

## Tools & Infrastructure

- **Test Framework**: Google Test (already in use)
- **Coverage Tool**: gcov/lcov for coverage analysis
- **Reference Data**: CSV files with known solutions
- **Validation**: Compare against GeographicLib where possible
- **CI Integration**: All tests must pass before merge

## Next Steps

1. Create `test_reference_data.h` with known coordinates
2. Implement `test_geodesy_reference.cpp`
3. Implement `test_geodesy_vincenty.cpp`
4. Implement `test_geodesy_rk4.cpp`
5. Continue with remaining test files
6. Generate coverage report
7. Fill any remaining gaps
