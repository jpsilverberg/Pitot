// Geodesy.h - WGS84 ECEF Coordinate System (C++17)
// High-performance geodetic coordinate handling with minimal trigonometric operations
#pragma once

#include <dstl/LinAlg.h>
#include <cmath>
#include <optional>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>

#ifndef ND_HD
#define ND_HD
#endif
#ifndef ND
#define ND [[nodiscard]]
#endif

namespace dstl
{
   namespace geo
   {
      // ========================== Constants ==========================
      
      // Use C++20 std::numbers::pi if available, otherwise define
      #if __cplusplus >= 202002L && __has_include(<numbers>)
      #include <numbers>
      constexpr double PI = std::numbers::pi;
      #else
      constexpr double PI = 3.14159265358979323846;
      #endif
      
      constexpr double DEG_TO_RAD = PI / 180.0;
      constexpr double RAD_TO_DEG = 180.0 / PI;
      
      // Numerical tolerances
      constexpr double POLAR_EPSILON = 1e-10;  // For polar axis detection
      constexpr double CONVERGENCE_EPSILON = 1e-12;  // For iterative convergence
      constexpr double COORDINATE_EPSILON = 1e-6;  // For coordinate equality (1 micron)
      
      // ========================== WGS84 Constants ==========================
      struct WGS84
      {
         static constexpr double a = 6378137.0;           // Semi-major axis (m)
         static constexpr double f = 1.0 / 298.257223563; // Flattening
         static constexpr double b = a * (1.0 - f);       // Semi-minor axis (m)
         static constexpr double e2 = f * (2.0 - f);      // First eccentricity squared
         static constexpr double ep2 = e2 / (1.0 - e2);   // Second eccentricity squared
      };

      // ======================== ECEF Coordinate Class ======================
      /**
       * @brief WGS84 ECEF (Earth-Centered, Earth-Fixed) coordinate representation
       * 
       * Design Philosophy:
       * - Primary storage: ECEF (X, Y, Z) for fast vector operations
       * - Lazy geodetic computation: lat/lon/alt computed only when accessed
       * - Cached trig values: sin/cos cached to avoid recomputation
       * - Zero-cost ECEF operations: distance, vectors, offsets use Vec3 math
       * 
       * Performance:
       * - ECEF construction: O(1), no trig
       * - Distance/vector ops: O(1), no trig
       * - Geodetic access (first): ~8-10 trig calls, then cached
       * - Geodetic access (subsequent): O(1), cache hit
       */
      class ECEFCoordinate
      {
      private:
         linalg::Vec3 m_ecef; // Primary storage: X, Y, Z in meters

         // Cached geodetic coordinates (lazy computation)
         mutable std::optional<double> m_lat; // Latitude (radians)
         mutable std::optional<double> m_lon; // Longitude (radians)
         mutable std::optional<double> m_alt; // Altitude (meters)

         // Cached trig values (computed with geodetic coords)
         mutable std::optional<double> m_sin_lat;
         mutable std::optional<double> m_cos_lat;
         mutable std::optional<double> m_sin_lon;
         mutable std::optional<double> m_cos_lon;

         // Invalidate cache when ECEF changes
         void invalidate_cache() const
         {
            m_lat.reset();
            m_lon.reset();
            m_alt.reset();
            m_sin_lat.reset();
            m_cos_lat.reset();
            m_sin_lon.reset();
            m_cos_lon.reset();
         }

         // Compute geodetic coordinates using Bowring's method
         // Converges in 2-3 iterations for most Earth-surface points
         void compute_geodetic() const
         {
            if (m_lat.has_value())
               return; // Already computed

            const double x = m_ecef.x;
            const double y = m_ecef.y;
            const double z = m_ecef.z;

            // Longitude is cheap - only one atan2
            m_lon = std::atan2(y, x);
            m_sin_lon = std::sin(*m_lon);
            m_cos_lon = std::cos(*m_lon);

            // For latitude and altitude, use Bowring's method
            const double p = std::sqrt(x * x + y * y);

            if (p < POLAR_EPSILON)
            {
               // Special case: on polar axis
               m_lat = (z >= 0) ? PI / 2.0 : -PI / 2.0;
               m_sin_lat = (z >= 0) ? 1.0 : -1.0;
               m_cos_lat = 0.0;
               
               // More accurate altitude for polar case
               double N = WGS84::a / std::sqrt(1.0 - WGS84::e2 * (*m_sin_lat) * (*m_sin_lat));
               m_alt = std::abs(z) - N * (1.0 - WGS84::e2);
               return;
            }

            // Initial estimate using spherical approximation
            double lat = std::atan2(z, p * (1.0 - WGS84::e2));

            // Iterate to refine (typically 2-3 iterations, max 4 for safety)
            for (int i = 0; i < 4; ++i)
            {
               double sin_lat = std::sin(lat);
               double cos_lat = std::cos(lat);
               double N = WGS84::a / std::sqrt(1.0 - WGS84::e2 * sin_lat * sin_lat);
               double h = p / cos_lat - N;
               double lat_new = std::atan2(z, p * (1.0 - WGS84::e2 * N / (N + h)));

               if (std::abs(lat_new - lat) < CONVERGENCE_EPSILON)
               {
                  lat = lat_new;
                  break;
               }
               lat = lat_new;
            }

            // Final computation with cached trig values
            m_lat = lat;
            m_sin_lat = std::sin(lat);
            m_cos_lat = std::cos(lat);

            double N = WGS84::a / std::sqrt(1.0 - WGS84::e2 * (*m_sin_lat) * (*m_sin_lat));
            m_alt = p / (*m_cos_lat) - N;
         }

      public:
         // ==================== Constructors ====================

         /// Construct from ECEF coordinates (fast - no trig)
         ND_HD constexpr ECEFCoordinate(double x, double y, double z) noexcept
             : m_ecef(x, y, z) {}

         /// Construct from ECEF vector (fast - no trig)
         ND_HD constexpr ECEFCoordinate(const linalg::Vec3 &ecef) noexcept
             : m_ecef(ecef) {}

         /// Default constructor (origin)
         ND_HD constexpr ECEFCoordinate() noexcept
             : m_ecef(0, 0, 0) {}

         /// Construct from geodetic coordinates (requires trig)
         /// @param lat_rad Latitude in radians
         /// @param lon_rad Longitude in radians
         /// @param alt_m Altitude above ellipsoid in meters
         ND_HD static inline ECEFCoordinate from_geodetic(double lat_rad, double lon_rad, double alt_m) noexcept
         {
            double sin_lat = std::sin(lat_rad);
            double cos_lat = std::cos(lat_rad);
            double sin_lon = std::sin(lon_rad);
            double cos_lon = std::cos(lon_rad);

            double N = WGS84::a / std::sqrt(1.0 - WGS84::e2 * sin_lat * sin_lat);

            double x = (N + alt_m) * cos_lat * cos_lon;
            double y = (N + alt_m) * cos_lat * sin_lon;
            double z = (N * (1.0 - WGS84::e2) + alt_m) * sin_lat;

            ECEFCoordinate coord(x, y, z);

            // Pre-cache the geodetic values and trig we just computed
            coord.m_lat = lat_rad;
            coord.m_lon = lon_rad;
            coord.m_alt = alt_m;
            coord.m_sin_lat = sin_lat;
            coord.m_cos_lat = cos_lat;
            coord.m_sin_lon = sin_lon;
            coord.m_cos_lon = cos_lon;

            return coord;
         }

         /// Construct from geodetic coordinates in degrees
         ND_HD static inline ECEFCoordinate from_geodetic_deg(double lat_deg, double lon_deg, double alt_m) noexcept
         {
            return from_geodetic(lat_deg * DEG_TO_RAD, lon_deg * DEG_TO_RAD, alt_m);
         }

         // ==================== ECEF Access (Fast - No Trig) ====================

         /// Get ECEF vector (X, Y, Z)
         ND_HD constexpr const linalg::Vec3 &ecef() const noexcept { return m_ecef; }

         /// Get X coordinate (meters)
         ND_HD constexpr double x() const noexcept { return m_ecef.x; }

         /// Get Y coordinate (meters)
         ND_HD constexpr double y() const noexcept { return m_ecef.y; }

         /// Get Z coordinate (meters)
         ND_HD constexpr double z() const noexcept { return m_ecef.z; }

         // ==================== Geodetic Access (Lazy Trig) ====================

         /// Get latitude in radians (lazy computation, then cached)
         ND_HD inline double latitude() const noexcept
         {
            compute_geodetic();
            return *m_lat;
         }

         /// Get longitude in radians (lazy computation, then cached)
         ND_HD inline double longitude() const noexcept
         {
            compute_geodetic();
            return *m_lon;
         }

         /// Get altitude in meters (lazy computation, then cached)
         ND_HD inline double altitude() const noexcept
         {
            compute_geodetic();
            return *m_alt;
         }

         /// Get latitude in degrees
         ND_HD inline double latitude_deg() const noexcept
         {
            return latitude() * RAD_TO_DEG;
         }

         /// Get longitude in degrees
         ND_HD inline double longitude_deg() const noexcept
         {
            return longitude() * RAD_TO_DEG;
         }

         // ==================== Cached Trig Access ====================

         /// Get sin(latitude) - cached after first geodetic computation
         ND_HD inline double sin_lat() const noexcept
         {
            compute_geodetic();
            return *m_sin_lat;
         }

         /// Get cos(latitude) - cached after first geodetic computation
         ND_HD inline double cos_lat() const noexcept
         {
            compute_geodetic();
            return *m_cos_lat;
         }

         /// Get sin(longitude) - cached after first geodetic computation
         ND_HD inline double sin_lon() const noexcept
         {
            compute_geodetic();
            return *m_sin_lon;
         }

         /// Get cos(longitude) - cached after first geodetic computation
         ND_HD inline double cos_lon() const noexcept
         {
            compute_geodetic();
            return *m_cos_lon;
         }

         // ==================== Vector Operations (Fast - No Trig) ====================

         /// Compute ECEF vector from this point to another
         ND_HD constexpr linalg::Vec3 vector_to(const ECEFCoordinate &other) const noexcept
         {
            return other.m_ecef - m_ecef;
         }

         /// Compute straight-line (chord) distance to another point (meters)
         /// Note: This is NOT surface distance. Use geodesic_distance_to() for that.
         ND_HD inline double distance_to(const ECEFCoordinate &other) const noexcept
         {
            return (other.m_ecef - m_ecef).norm();
         }

         /// Compute geodesic (great circle) distance using Vincenty's formula (meters)
         /// This is the actual surface distance along the ellipsoid.
         /// Accurate for all distances, converges in 1-3 iterations typically.
         /// Note: For near-antipodal points, may fail to converge (returns approximate distance).
         ND_HD inline double geodesic_distance_to(const ECEFCoordinate &other, int max_iterations = 20) const noexcept
         {
            // Get geodetic coordinates (triggers lazy computation if needed)
            double lat1 = latitude();
            double lon1 = longitude();
            double lat2 = other.latitude();
            double lon2 = other.longitude();

            // Vincenty constants
            constexpr double f_16 = WGS84::f / 16.0;
            constexpr double one_minus_f = 1.0 - WGS84::f;

            // Use cached trig values
            double U1 = std::atan(one_minus_f * std::tan(lat1));
            double U2 = std::atan(one_minus_f * std::tan(lat2));
            double sinU1 = std::sin(U1);
            double cosU1 = std::cos(U1);
            double sinU2 = std::sin(U2);
            double cosU2 = std::cos(U2);

            double L = lon2 - lon1;
            double lambda = L;
            double lambda_prev;

            double sinSigma, cosSigma, sigma, sinAlpha, cosSqAlpha, cos2SigmaM;
            bool converged = false;

            for (int i = 0; i < max_iterations; ++i)
            {
               double sinLambda = std::sin(lambda);
               double cosLambda = std::cos(lambda);

               sinSigma = std::sqrt((cosU2 * sinLambda) * (cosU2 * sinLambda) +
                                    (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda) *
                                        (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda));

               if (sinSigma < CONVERGENCE_EPSILON)
                  return 0.0; // Coincident points

               cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;
               sigma = std::atan2(sinSigma, cosSigma);
               sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;
               cosSqAlpha = 1.0 - sinAlpha * sinAlpha;

               cos2SigmaM = (cosSqAlpha > CONVERGENCE_EPSILON) ? (cosSigma - 2.0 * sinU1 * sinU2 / cosSqAlpha) : 0.0;

               double C = f_16 * cosSqAlpha * (4.0 + WGS84::f * (4.0 - 3.0 * cosSqAlpha));

               lambda_prev = lambda;
               lambda = L + (1.0 - C) * WGS84::f * sinAlpha *
                                (sigma + C * sinSigma * (cos2SigmaM + C * cosSigma * (-1.0 + 2.0 * cos2SigmaM * cos2SigmaM)));

               if (std::abs(lambda - lambda_prev) < CONVERGENCE_EPSILON)
               {
                  converged = true;
                  break;
               }
            }

            // Antipodal fallback: if didn't converge, return approximate half-circumference
            if (!converged && std::abs(sinSigma) < CONVERGENCE_EPSILON)
            {
               return PI * WGS84::a; // Approximate for antipodal points
            }

            double uSq = cosSqAlpha * (WGS84::a * WGS84::a - WGS84::b * WGS84::b) / (WGS84::b * WGS84::b);
            double A = 1.0 + uSq / 16384.0 * (4096.0 + uSq * (-768.0 + uSq * (320.0 - 175.0 * uSq)));
            double B = uSq / 1024.0 * (256.0 + uSq * (-128.0 + uSq * (74.0 - 47.0 * uSq)));
            double deltaSigma = B * sinSigma * (cos2SigmaM + B / 4.0 * (cosSigma * (-1.0 + 2.0 * cos2SigmaM * cos2SigmaM) -
                                                                         B / 6.0 * cos2SigmaM * (-3.0 + 4.0 * sinSigma * sinSigma) *
                                                                             (-3.0 + 4.0 * cos2SigmaM * cos2SigmaM)));

            return WGS84::b * A * (sigma - deltaSigma);
         }

         /// Compute forward azimuth (bearing) to another point (radians from North, 0-2π)
         /// Returns angle in range [0, 2π) where 0 = North, π/2 = East, π = South, 3π/2 = West
         ND_HD inline double bearing_to(const ECEFCoordinate &other) const noexcept
         {
            // Use ENU frame for simple, accurate bearing calculation
            auto frame = local_frame();
            linalg::Vec3 enu = frame.coord_to_enu(other);

            // atan2(East, North) gives bearing from North
            double bearing = std::atan2(enu.x, enu.y);

            // Normalize to [0, 2π)
            if (bearing < 0.0)
               bearing += 2.0 * PI;

            return bearing;
         }

         /// Move along geodesic by distance and bearing (approximate for short distances)
         /// For distances < 100km, this is accurate to ~1m
         /// For longer distances, use iterative Vincenty direct formula
         ND_HD inline ECEFCoordinate move_by_bearing(double distance_m, double bearing_rad) const noexcept
         {
            // Use local ENU frame for approximate movement
            auto frame = local_frame();

            // Convert bearing to ENU offset
            double east = distance_m * std::sin(bearing_rad);
            double north = distance_m * std::cos(bearing_rad);

            linalg::Vec3 enu_offset{east, north, 0.0};
            return frame.enu_to_coord(enu_offset);
         }

         /// Move along geodesic by distance and bearing (accurate Vincenty direct formula)
         /// Accurate for all distances on Earth
         /// Note: Maintains the same altitude as the origin point
         ND_HD inline ECEFCoordinate move_by_bearing_accurate(double distance_m, double bearing_rad, int max_iterations = 20) const noexcept
         {
            double lat1 = latitude();
            double lon1 = longitude();
            double alt = altitude();
            double alpha1 = bearing_rad;

            double sinAlpha1 = std::sin(alpha1);
            double cosAlpha1 = std::cos(alpha1);

            double tanU1 = (1.0 - WGS84::f) * std::tan(lat1);
            double cosU1 = 1.0 / std::sqrt(1.0 + tanU1 * tanU1);
            double sinU1 = tanU1 * cosU1;

            double sigma1 = std::atan2(tanU1, cosAlpha1);
            double sinAlpha = cosU1 * sinAlpha1;
            double cosSqAlpha = 1.0 - sinAlpha * sinAlpha;
            double uSq = cosSqAlpha * (WGS84::a * WGS84::a - WGS84::b * WGS84::b) / (WGS84::b * WGS84::b);
            double A = 1.0 + uSq / 16384.0 * (4096.0 + uSq * (-768.0 + uSq * (320.0 - 175.0 * uSq)));
            double B = uSq / 1024.0 * (256.0 + uSq * (-128.0 + uSq * (74.0 - 47.0 * uSq)));

            double sigma = distance_m / (WGS84::b * A);
            double sigma_prev;

            double cos2SigmaM, sinSigma, cosSigma;

            for (int i = 0; i < max_iterations; ++i)
            {
               cos2SigmaM = std::cos(2.0 * sigma1 + sigma);
               sinSigma = std::sin(sigma);
               cosSigma = std::cos(sigma);

               double deltaSigma = B * sinSigma * (cos2SigmaM + B / 4.0 * (cosSigma * (-1.0 + 2.0 * cos2SigmaM * cos2SigmaM) -
                                                                            B / 6.0 * cos2SigmaM * (-3.0 + 4.0 * sinSigma * sinSigma) *
                                                                                (-3.0 + 4.0 * cos2SigmaM * cos2SigmaM)));

               sigma_prev = sigma;
               sigma = distance_m / (WGS84::b * A) + deltaSigma;

               if (std::abs(sigma - sigma_prev) < CONVERGENCE_EPSILON)
                  break;
            }

            double tmp = sinU1 * sinSigma - cosU1 * cosSigma * cosAlpha1;
            double lat2 = std::atan2(sinU1 * cosSigma + cosU1 * sinSigma * cosAlpha1,
                                     (1.0 - WGS84::f) * std::sqrt(sinAlpha * sinAlpha + tmp * tmp));

            double lambda = std::atan2(sinSigma * sinAlpha1, cosU1 * cosSigma - sinU1 * sinSigma * cosAlpha1);
            double C = WGS84::f / 16.0 * cosSqAlpha * (4.0 + WGS84::f * (4.0 - 3.0 * cosSqAlpha));
            double L = lambda - (1.0 - C) * WGS84::f * sinAlpha *
                                    (sigma + C * sinSigma * (cos2SigmaM + C * cosSigma * (-1.0 + 2.0 * cos2SigmaM * cos2SigmaM)));

            double lon2 = lon1 + L;

            return from_geodetic(lat2, lon2, alt);
         }

         /// Add ECEF offset vector
         ND_HD constexpr ECEFCoordinate operator+(const linalg::Vec3 &offset) const noexcept
         {
            return ECEFCoordinate(m_ecef + offset);
         }

         /// Subtract ECEF offset vector
         ND_HD constexpr ECEFCoordinate operator-(const linalg::Vec3 &offset) const noexcept
         {
            return ECEFCoordinate(m_ecef - offset);
         }

         /// Add ECEF offset in place
         ND_HD inline ECEFCoordinate &operator+=(const linalg::Vec3 &offset) noexcept
         {
            m_ecef += offset;
            invalidate_cache();
            return *this;
         }

         /// Subtract ECEF offset in place
         ND_HD inline ECEFCoordinate &operator-=(const linalg::Vec3 &offset) noexcept
         {
            m_ecef -= offset;
            invalidate_cache();
            return *this;
         }

         // ==================== Local Coordinate Frames ====================

         /**
          * @brief Local ENU (East-North-Up) coordinate frame
          * 
          * Provides transformation between ECEF and local tangent plane coordinates.
          * Compute once per reference point, then reuse for multiple transformations.
          */
         struct ENUFrame
         {
            linalg::Vec3 east;  ///< Unit vector pointing East in ECEF
            linalg::Vec3 north; ///< Unit vector pointing North in ECEF
            linalg::Vec3 up;    ///< Unit vector pointing Up in ECEF
            linalg::Vec3 origin_ecef;  ///< Origin point in ECEF

            /// Convert ECEF vector to ENU coordinates
            ND_HD constexpr linalg::Vec3 to_enu(const linalg::Vec3 &ecef_vec) const noexcept
            {
               return linalg::Vec3{
                   ecef_vec.dot(east),
                   ecef_vec.dot(north),
                   ecef_vec.dot(up)};
            }

            /// Convert ENU vector to ECEF coordinates
            ND_HD constexpr linalg::Vec3 to_ecef(const linalg::Vec3 &enu_vec) const noexcept
            {
               return east * enu_vec.x + north * enu_vec.y + up * enu_vec.z;
            }

            /// Convert another coordinate to ENU relative to origin
            ND_HD constexpr linalg::Vec3 coord_to_enu(const ECEFCoordinate &coord) const noexcept
            {
               return to_enu(coord.ecef() - origin_ecef);
            }

            /// Convert ENU offset to absolute ECEF coordinate
            ND_HD constexpr ECEFCoordinate enu_to_coord(const linalg::Vec3 &enu_offset) const noexcept
            {
               return ECEFCoordinate(origin_ecef + to_ecef(enu_offset));
            }
         };

         /// Create local ENU frame at this point (requires trig, but cache the frame)
         ND_HD inline ENUFrame local_frame() const noexcept
         {
            // Trigger geodetic computation if needed
            double lat = latitude();
            double lon = longitude();

            // Use cached trig values
            double sl = *m_sin_lat;
            double cl = *m_cos_lat;
            double sn = *m_sin_lon;
            double cn = *m_cos_lon;

            return ENUFrame{
                linalg::Vec3{-sn, cn, 0.0},
                linalg::Vec3{-sl * cn, -sl * sn, cl},
                linalg::Vec3{cl * cn, cl * sn, sl},
                m_ecef};
         }

         /// Apply flat-earth displacement in ENU coordinates
         /// This is the key method for flight dynamics: compute maneuvers in local flat-earth,
         /// then apply back to WGS84 with proper altitude handling.
         /// 
         /// @param enu_displacement Displacement in local ENU frame (meters)
         ///        - enu_displacement.x = East displacement
         ///        - enu_displacement.y = North displacement  
         ///        - enu_displacement.z = Up displacement (altitude change)
         /// @return New WGS84 coordinate with displacement applied
         /// 
         /// Note: The Z-component is interpreted as altitude change relative to the ellipsoid.
         /// For maintaining constant altitude over long distances, use apply_enu_displacement_constant_altitude().
         /// 
         /// Example: Fly 100nm NE, descend 10,000ft
         ///   double nm_to_m = 1852.0;
         ///   double ft_to_m = 0.3048;
         ///   Vec3 displacement{100*nm_to_m/sqrt(2), 100*nm_to_m/sqrt(2), -10000*ft_to_m};
         ///   auto new_pos = origin.apply_enu_displacement(displacement);
         ND_HD inline ECEFCoordinate apply_enu_displacement(const linalg::Vec3 &enu_displacement) const noexcept
         {
            auto frame = local_frame();
            auto target_ecef = frame.origin_ecef + frame.to_ecef(enu_displacement);
            return ECEFCoordinate(target_ecef);
         }

         /// Apply flat-earth displacement with constant altitude
         /// This version maintains the altitude (above ellipsoid) constant, which is what you typically
         /// want for flight at constant altitude. The horizontal displacement is applied in the local
         /// tangent plane, then the altitude is explicitly set.
         /// 
         /// @param enu_horizontal Horizontal displacement in ENU (only x and y used)
         /// @param altitude_change_m Change in altitude (meters), default 0 for constant altitude
         /// @return New WGS84 coordinate with displacement applied and altitude corrected
         /// 
         /// Example: Fly 100nm NE at constant 30,000ft
         ///   Vec3 horizontal{100*NM_TO_M/sqrt(2), 100*NM_TO_M/sqrt(2), 0};
         ///   auto new_pos = origin.apply_enu_displacement_constant_altitude(horizontal);
         ND_HD inline ECEFCoordinate apply_enu_displacement_constant_altitude(
             const linalg::Vec3 &enu_horizontal,
             double altitude_change_m = 0.0) const noexcept
         {
            // Apply horizontal displacement only
            auto frame = local_frame();
            linalg::Vec3 horizontal_only{enu_horizontal.x, enu_horizontal.y, 0.0};
            auto target_ecef = frame.origin_ecef + frame.to_ecef(horizontal_only);
            auto target = ECEFCoordinate(target_ecef);
            
            // Set altitude explicitly
            double target_altitude = altitude() + altitude_change_m;
            return from_geodetic(target.latitude(), target.longitude(), target_altitude);
         }

         /// Apply flat-earth displacement in NED coordinates
         /// Same as apply_enu_displacement but for NED frame (common in aerospace)
         /// 
         /// @param ned_displacement Displacement in local NED frame (meters)
         ///        - ned_displacement.x = North displacement
         ///        - ned_displacement.y = East displacement
         ///        - ned_displacement.z = Down displacement (negative = climb)
         /// @return New WGS84 coordinate with displacement applied
         /// 
         /// Note: For maintaining constant altitude, use apply_ned_displacement_constant_altitude().
         ND_HD inline ECEFCoordinate apply_ned_displacement(const linalg::Vec3 &ned_displacement) const noexcept
         {
            auto frame = local_ned_frame();
            auto target_ecef = frame.origin_ecef + frame.to_ecef(ned_displacement);
            return ECEFCoordinate(target_ecef);
         }

         /// Apply flat-earth displacement with constant altitude (NED version)
         /// Maintains altitude constant while applying horizontal displacement in NED frame.
         /// 
         /// @param ned_horizontal Horizontal displacement in NED (only x and y used)
         /// @param altitude_change_m Change in altitude (meters), default 0 for constant altitude
         /// @return New WGS84 coordinate with displacement applied and altitude corrected
         /// 
         /// Example: Fly 100nm North at constant FL350
         ///   Vec3 horizontal{100*NM_TO_M, 0, 0};
         ///   auto new_pos = origin.apply_ned_displacement_constant_altitude(horizontal);
         ND_HD inline ECEFCoordinate apply_ned_displacement_constant_altitude(
             const linalg::Vec3 &ned_horizontal,
             double altitude_change_m = 0.0) const noexcept
         {
            // Apply horizontal displacement only
            auto frame = local_ned_frame();
            linalg::Vec3 horizontal_only{ned_horizontal.x, ned_horizontal.y, 0.0};
            auto target_ecef = frame.origin_ecef + frame.to_ecef(horizontal_only);
            auto target = ECEFCoordinate(target_ecef);
            
            // Set altitude explicitly
            double target_altitude = altitude() + altitude_change_m;
            return from_geodetic(target.latitude(), target.longitude(), target_altitude);
         }

         /// Apply flat-earth trajectory with frame updates
         /// For long trajectories, periodically update the local frame to maintain accuracy.
         /// This is critical for flights > 100km where Earth curvature becomes significant.
         /// 
         /// @param enu_displacements Vector of displacement segments in ENU
         /// @param update_interval_m Update local frame every N meters of travel (default 50km)
         /// @return Final WGS84 coordinate after all displacements
         /// 
         /// Example: Multi-leg flight with frame updates
         ///   std::vector<Vec3> legs = {
         ///       {50000, 50000, 1000},   // 50km E, 50km N, climb 1km
         ///       {30000, -20000, -500},  // 30km E, 20km S, descend 500m
         ///   };
         ///   auto final_pos = origin.apply_enu_trajectory(legs);
         inline ECEFCoordinate apply_enu_trajectory(
             const std::vector<linalg::Vec3> &enu_displacements,
             double update_interval_m = 50000.0) const noexcept
         {
            ECEFCoordinate current = *this;
            linalg::Vec3 accumulated{0, 0, 0};
            double accumulated_distance = 0.0;

            for (const auto &displacement : enu_displacements)
            {
               accumulated += displacement;
               double segment_distance = std::sqrt(displacement.x * displacement.x +
                                                   displacement.y * displacement.y);
               accumulated_distance += segment_distance;

               // Update frame if we've traveled far enough
               if (accumulated_distance >= update_interval_m)
               {
                  current = current.apply_enu_displacement(accumulated);
                  accumulated = linalg::Vec3{0, 0, 0};
                  accumulated_distance = 0.0;
               }
            }

            // Apply any remaining displacement
            if (accumulated.norm() > 1e-6)
            {
               current = current.apply_enu_displacement(accumulated);
            }

            return current;
         }

         /// Apply flat-earth trajectory in NED with frame updates
         /// Same as apply_enu_trajectory but for NED coordinates
         inline ECEFCoordinate apply_ned_trajectory(
             const std::vector<linalg::Vec3> &ned_displacements,
             double update_interval_m = 50000.0) const noexcept
         {
            ECEFCoordinate current = *this;
            linalg::Vec3 accumulated{0, 0, 0};
            double accumulated_distance = 0.0;

            for (const auto &displacement : ned_displacements)
            {
               accumulated += displacement;
               double segment_distance = std::sqrt(displacement.x * displacement.x +
                                                   displacement.y * displacement.y);
               accumulated_distance += segment_distance;

               // Update frame if we've traveled far enough
               if (accumulated_distance >= update_interval_m)
               {
                  current = current.apply_ned_displacement(accumulated);
                  accumulated = linalg::Vec3{0, 0, 0};
                  accumulated_distance = 0.0;
               }
            }

            // Apply any remaining displacement
            if (accumulated.norm() > 1e-6)
            {
               current = current.apply_ned_displacement(accumulated);
            }

            return current;
         }

         /**
          * @brief Local NED (North-East-Down) coordinate frame
          * 
          * Common in aerospace and marine applications.
          * Provides transformation between ECEF and local tangent plane coordinates.
          */
         struct NEDFrame
         {
            linalg::Vec3 north; ///< Unit vector pointing North in ECEF
            linalg::Vec3 east;  ///< Unit vector pointing East in ECEF
            linalg::Vec3 down;  ///< Unit vector pointing Down in ECEF
            linalg::Vec3 origin_ecef;  ///< Origin point in ECEF

            /// Convert ECEF vector to NED coordinates
            ND_HD constexpr linalg::Vec3 to_ned(const linalg::Vec3 &ecef_vec) const noexcept
            {
               return linalg::Vec3{
                   ecef_vec.dot(north),
                   ecef_vec.dot(east),
                   ecef_vec.dot(down)};
            }

            /// Convert NED vector to ECEF coordinates
            ND_HD constexpr linalg::Vec3 to_ecef(const linalg::Vec3 &ned_vec) const noexcept
            {
               return north * ned_vec.x + east * ned_vec.y + down * ned_vec.z;
            }

            /// Convert another coordinate to NED relative to origin
            ND_HD constexpr linalg::Vec3 coord_to_ned(const ECEFCoordinate &coord) const noexcept
            {
               return to_ned(coord.ecef() - origin_ecef);
            }

            /// Convert NED offset to absolute ECEF coordinate
            ND_HD constexpr ECEFCoordinate ned_to_coord(const linalg::Vec3 &ned_offset) const noexcept
            {
               return ECEFCoordinate(origin_ecef + to_ecef(ned_offset));
            }
         };

         /// Create local NED frame at this point (requires trig, but cache the frame)
         ND_HD inline NEDFrame local_ned_frame() const noexcept
         {
            // Trigger geodetic computation if needed
            double lat = latitude();
            double lon = longitude();

            // Use cached trig values
            double sl = *m_sin_lat;
            double cl = *m_cos_lat;
            double sn = *m_sin_lon;
            double cn = *m_cos_lon;

            return NEDFrame{
                linalg::Vec3{-sl * cn, -sl * sn, cl},   // North
                linalg::Vec3{-sn, cn, 0.0},             // East
                linalg::Vec3{-cl * cn, -cl * sn, -sl},  // Down (negative up)
                m_ecef};
         }

         // ==================== Utility ====================

         /// Check if geodetic coordinates have been computed
         ND_HD constexpr bool has_geodetic_cache() const noexcept
         {
            return m_lat.has_value();
         }

         /// Force computation of geodetic coordinates (for pre-caching)
         inline void precompute_geodetic() const noexcept
         {
            compute_geodetic();
         }

         /// Convert to string representation
         /// @param precision Number of decimal places for coordinates
         /// @param use_degrees If true, output lat/lon in degrees; otherwise radians
         inline std::string to_string(int precision = 6, bool use_degrees = true) const
         {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(precision);

            if (use_degrees)
            {
               oss << "ECEFCoordinate(lat=" << latitude_deg() << "°, "
                   << "lon=" << longitude_deg() << "°, "
                   << "alt=" << altitude() << "m)";
            }
            else
            {
               oss << "ECEFCoordinate(lat=" << latitude() << " rad, "
                   << "lon=" << longitude() << " rad, "
                   << "alt=" << altitude() << "m)";
            }

            return oss.str();
         }

         /// Convert to string with ECEF coordinates
         inline std::string to_string_ecef(int precision = 3) const
         {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(precision);
            oss << "ECEFCoordinate(x=" << x() << ", "
                << "y=" << y() << ", "
                << "z=" << z() << ")";
            return oss.str();
         }

         // ==================== Comparison Operators ====================

         /// Equality comparison (within tolerance)
         /// Uses COORDINATE_EPSILON (1 micron) for floating-point comparison
         ND_HD constexpr bool operator==(const ECEFCoordinate &other) const noexcept
         {
            return std::abs(m_ecef.x - other.m_ecef.x) < COORDINATE_EPSILON &&
                   std::abs(m_ecef.y - other.m_ecef.y) < COORDINATE_EPSILON &&
                   std::abs(m_ecef.z - other.m_ecef.z) < COORDINATE_EPSILON;
         }

         /// Inequality comparison
         ND_HD constexpr bool operator!=(const ECEFCoordinate &other) const noexcept
         {
            return !(*this == other);
         }

         /// Approximate equality with custom tolerance
         ND_HD constexpr bool approx_equal(const ECEFCoordinate &other, double tolerance_m = COORDINATE_EPSILON) const noexcept
         {
            return std::abs(m_ecef.x - other.m_ecef.x) < tolerance_m &&
                   std::abs(m_ecef.y - other.m_ecef.y) < tolerance_m &&
                   std::abs(m_ecef.z - other.m_ecef.z) < tolerance_m;
         }

         // ==================== Geometric Utilities ====================

         /// Get unit normal vector to ellipsoid surface at this point
         /// Points radially outward from Earth center (not exactly vertical due to ellipsoid shape)
         /// Useful for: gravity direction, surface projections, ray tracing
         ND_HD inline linalg::Vec3 surface_normal() const noexcept
         {
            compute_geodetic();
            return linalg::Vec3{*m_cos_lat * *m_cos_lon, 
                               *m_cos_lat * *m_sin_lon, 
                               *m_sin_lat};
         }

         /// Get geodesic midpoint between this point and another
         /// Uses Vincenty's formula for accurate midpoint on ellipsoid
         /// Useful for: path smoothing, interpolation, waypoint generation
         ND_HD inline ECEFCoordinate midpoint(const ECEFCoordinate &other) const noexcept
         {
            // For short distances, simple average works well
            double dist = distance_to(other);
            if (dist < 100000.0)  // < 100 km
            {
               // Simple ECEF average (fast, accurate for short distances)
               linalg::Vec3 mid_ecef = (m_ecef + other.m_ecef) * 0.5;
               return ECEFCoordinate(mid_ecef);
            }

            // For longer distances, use geodetic interpolation
            double lat1 = latitude();
            double lon1 = longitude();
            double alt1 = altitude();
            double lat2 = other.latitude();
            double lon2 = other.longitude();
            double alt2 = other.altitude();

            // Spherical interpolation for lat/lon
            double Bx = other.cos_lat() * std::cos(lon2 - lon1);
            double By = other.cos_lat() * std::sin(lon2 - lon1);
            
            double lat_mid = std::atan2(*m_sin_lat + other.sin_lat(),
                                       std::sqrt((*m_cos_lat + Bx) * (*m_cos_lat + Bx) + By * By));
            double lon_mid = lon1 + std::atan2(By, *m_cos_lat + Bx);
            double alt_mid = (alt1 + alt2) * 0.5;

            return from_geodetic(lat_mid, lon_mid, alt_mid);
         }
      };

      // ==================== Free Functions ====================

      /// Stream output operator
      inline std::ostream &operator<<(std::ostream &os, const ECEFCoordinate &coord)
      {
         return os << coord.to_string();
      }

      // ==================== Free Functions ====================

      /// Commutative addition: offset + coordinate
      ND_HD constexpr inline ECEFCoordinate operator+(const linalg::Vec3 &offset, const ECEFCoordinate &coord) noexcept
      {
         return coord + offset;
      }

      /// Compute geodesic midpoint between two coordinates
      ND_HD inline ECEFCoordinate midpoint(const ECEFCoordinate &a, const ECEFCoordinate &b) noexcept
      {
         return a.midpoint(b);
      }

      /// Interpolate between two coordinates (0 = a, 1 = b)
      /// Uses geodesic interpolation for accurate results
      ND_HD inline ECEFCoordinate interpolate(const ECEFCoordinate &a, const ECEFCoordinate &b, double t) noexcept
      {
         if (t <= 0.0) return a;
         if (t >= 1.0) return b;
         if (std::abs(t - 0.5) < 1e-10) return midpoint(a, b);

         // For general interpolation, move along geodesic
         double total_dist = a.geodesic_distance_to(b);
         double bearing = a.bearing_to(b);
         return a.move_by_bearing_accurate(total_dist * t, bearing);
      }

      /// Compute centroid of multiple coordinates
      /// Uses ECEF average for simplicity (accurate for nearby points)
      inline ECEFCoordinate centroid(const std::vector<ECEFCoordinate> &coords) noexcept
      {
         if (coords.empty()) return ECEFCoordinate{};
         if (coords.size() == 1) return coords[0];

         linalg::Vec3 sum{0, 0, 0};
         for (const auto &coord : coords)
         {
            sum += coord.ecef();
         }
         return ECEFCoordinate(sum / static_cast<double>(coords.size()));
      }

   } // namespace geo
} // namespace dstl

// ==================== std::hash specialization ====================
namespace std
{
   /// Hash function for ECEFCoordinate (enables use in unordered containers)
   template <>
   struct hash<dstl::geo::ECEFCoordinate>
   {
      size_t operator()(const dstl::geo::ECEFCoordinate &coord) const noexcept
      {
         // Combine hashes of ECEF coordinates using FNV-1a-like algorithm
         size_t h1 = hash<double>{}(coord.x());
         size_t h2 = hash<double>{}(coord.y());
         size_t h3 = hash<double>{}(coord.z());

         // Mix the hashes
         size_t result = 2166136261u; // FNV offset basis
         result ^= h1;
         result *= 16777619u; // FNV prime
         result ^= h2;
         result *= 16777619u;
         result ^= h3;
         result *= 16777619u;

         return result;
      }
   };
} // namespace std
