// vec.h  (C++17)
// Zero-overhead vectors with CRTP base. No unions; no UB; standard layout.
#pragma once
#include <assert.h>
#include <cstddef>
#include <type_traits>
#include <algorithm>
#include <ostream>
#include <cmath>

#include <dstl/dprint.h>
#include <dstl/Numbers.h>

#ifndef ND_HD
#define ND_HD
#endif
#ifndef ND
#define ND [[nodiscard]]
#endif

namespace dstl
{
   namespace linalg
   {

      // =============================== VecBase =====================================
      template <class Derived, std::size_t N, class T>
      struct VecBase
      {
         using value_type = T;
         static constexpr std::size_t kN = N;

         // storage access that derived must implement
         ND_HD constexpr T *data() noexcept { return static_cast<Derived *>(this)->data_impl(); }
         ND_HD constexpr const T *data() const noexcept { return static_cast<const Derived *>(this)->data_impl(); }

         // element access / iterators
         ND_HD constexpr T &operator[](std::size_t i) noexcept
         {
            assert(i < N);
            return data()[i];
         }
         ND_HD constexpr const T &operator[](std::size_t i) const noexcept
         {
            assert(i < N);
            return data()[i];
         }
         ND_HD constexpr T *begin() noexcept { return data(); }
         ND_HD constexpr T *end() noexcept { return data() + N; }
         ND_HD constexpr const T *begin() const noexcept { return data(); }
         ND_HD constexpr const T *end() const noexcept { return data() + N; }

         // unary
         ND_HD constexpr Derived operator+() const noexcept { return derived(); }
         ND_HD constexpr Derived operator-() const noexcept
         {
            Derived r{};
            for (std::size_t i = 0; i < N; ++i)
               r[i] = -(*this)[i];
            return r;
         }

         // elementwise compound
         ND_HD constexpr Derived &operator+=(const Derived &b) noexcept
         {
            for (std::size_t i = 0; i < N; ++i)
               (*this)[i] += b[i];
            return derived();
         }
         ND_HD constexpr Derived &operator-=(const Derived &b) noexcept
         {
            for (std::size_t i = 0; i < N; ++i)
               (*this)[i] -= b[i];
            return derived();
         }
         ND_HD constexpr Derived &operator*=(const Derived &b) noexcept
         {
            for (std::size_t i = 0; i < N; ++i)
               (*this)[i] *= b[i];
            return derived();
         }
         ND_HD constexpr Derived &operator/=(const Derived &b) noexcept
         {
            for (std::size_t i = 0; i < N; ++i)
               (*this)[i] /= b[i];
            return derived();
         }

         // scalar compound
         template <class S, std::enable_if_t<std::is_arithmetic<S>::value, int> = 0>
         ND_HD constexpr Derived &operator*=(S s) noexcept
         {
            for (std::size_t i = 0; i < N; ++i)
               (*this)[i] = T((*this)[i] * s);
            return derived();
         }
         template <class S, std::enable_if_t<std::is_arithmetic<S>::value, int> = 0>
         ND_HD constexpr Derived &operator/=(S s) noexcept
         {
            for (std::size_t i = 0; i < N; ++i)
               (*this)[i] = T((*this)[i] / s);
            return derived();
         }

         // dot / norms
         ND_HD ND constexpr T dot(const Derived &b) const noexcept
         {
            Accum a{};
            // T a = T(0);
            for (std::size_t i = 0; i < N; ++i)
               a += Accum((*this)[i] * b[i]);
            return a.Double();
         }
         ND_HD ND constexpr T sqnorm() const noexcept { return dot(derived()); }
         ND_HD ND constexpr T norm() const noexcept
         {
            using std::sqrt;
            return T(sqrt(double(sqnorm())));
         }

         // normalized (runtime sqrt; if you want constexpr, add a constexpr sqrt)
         ND_HD ND constexpr Derived normalized(T eps = T(0)) const noexcept
         {
            T n = norm();
            return (n > eps) ? (derived() / n) : Derived{};
         }

         // utilities
         ND_HD static constexpr Derived zero() noexcept { return Derived{}; }
         ND_HD static constexpr Derived unit(std::size_t i) noexcept
         {
            Derived r{};
            if (i < N)
               r[i] = T(1);
            return r;
         }

         void fill(T s) noexcept
         {
            for (std::size_t i = 0; i < N; ++i)
               (*this)[i] = s;
         }

      private:
         ND_HD constexpr Derived &derived() noexcept { return *static_cast<Derived *>(this); }
         ND_HD constexpr const Derived &derived() const noexcept { return *static_cast<const Derived *>(this); }
      };

      // ========================== generic N storage =================================
      template <std::size_t N, class T = double>
      struct Vec : VecBase<Vec<N, T>, N, T>
      {
         using Base = VecBase<Vec<N, T>, N, T>;
         T a[N];

         // ctors
         ND_HD constexpr Vec() noexcept
         {
            for (std::size_t i = 0; i < N; ++i)
               a[i] = T(0);
         }
         ND_HD explicit constexpr Vec(T s) noexcept
         {
            for (std::size_t i = 0; i < N; ++i)
               a[i] = s;
         }
         template <class... U,
                   std::enable_if_t<(sizeof...(U) == N) && std::conjunction<std::is_convertible<U, T>...>::value, int> = 0>
         ND_HD constexpr Vec(U... xs) noexcept : a{T(xs)...} {}

         // CRTP storage hooks
         ND_HD constexpr T *data_impl() noexcept { return a; }
         ND_HD constexpr const T *data_impl() const noexcept { return a; }
      };

      // CTAD for Vec<N,T> from N args
      template <class... U>
      Vec(U...) -> Vec<sizeof...(U), std::common_type_t<U...>>;

      // ====================== 2/3/4 storage with named fields ======================
      struct Vec2 : VecBase<Vec2, 2, double>
      {
         using Base = VecBase<Vec2, 2, double>;
         double x, y; // Guaranteed to be laid out contiguously
         ND_HD constexpr Vec2() noexcept : x(0), y(0) {}
         ND_HD explicit constexpr Vec2(double s) noexcept : x(s), y(s) {}
         ND_HD constexpr Vec2(double x_, double y_) noexcept : x(x_), y(y_) {}
         // storage hook: strictly standard; no union/type-punning
         ND_HD constexpr double *data_impl() noexcept { return &x; }
         ND_HD constexpr const double *data_impl() const noexcept { return &x; }
      };

      struct Vec3 : VecBase<Vec3, 3, double>
      {
         using Base = VecBase<Vec3, 3, double>;
         double x, y, z; // Guaranteed to be laid out contiguously
         ND_HD constexpr Vec3() noexcept : x(0), y(0), z(0) {}
         ND_HD explicit constexpr Vec3(double s) noexcept : x(s), y(s), z(s) {}
         ND_HD constexpr Vec3(double x_, double y_, double z_) noexcept : x(x_), y(y_), z(z_) {}
         ND_HD constexpr double *data_impl() noexcept { return &x; }
         ND_HD constexpr const double *data_impl() const noexcept { return &x; }
      };

      struct Vec4 : VecBase<Vec4, 4, double>
      {
         using Base = VecBase<Vec4, 4, double>;
         double x, y, z, w; // Guaranteed to be laid out contiguously
         ND_HD constexpr Vec4() noexcept : x(0), y(0), z(0), w(0) {}
         ND_HD explicit constexpr Vec4(double s) noexcept : x(s), y(s), z(s), w(s) {}
         ND_HD constexpr Vec4(double x_, double y_, double z_, double w_) noexcept : x(x_), y(y_), z(z_), w(w_) {}
         ND_HD constexpr double *data_impl() noexcept { return &x; }
         ND_HD constexpr const double *data_impl() const noexcept { return &x; }
      };

      // ============================= free operators ================================
      template <class D, std::size_t N, class T>
      ND_HD ND constexpr D operator+(D a, const VecBase<D, N, T> &b) noexcept { return a += static_cast<const D &>(b), a; }
      template <class D, std::size_t N, class T>
      ND_HD ND constexpr D operator-(D a, const VecBase<D, N, T> &b) noexcept { return a -= static_cast<const D &>(b), a; }

      // v * s
      template <class D, std::size_t N, class T, class S,
                std::enable_if_t<std::is_arithmetic<S>::value, int> = 0>
      ND_HD ND constexpr D operator*(const VecBase<D, N, T> &a, S s) noexcept
      {
         D r = static_cast<const D &>(a);
         r *= s;
         return r;
      }

      // s * v
      template <class D, std::size_t N, class T, class S,
                std::enable_if_t<std::is_arithmetic<S>::value, int> = 0>
      ND_HD ND constexpr D operator*(S s, const VecBase<D, N, T> &a) noexcept
      {
         D r = static_cast<const D &>(a);
         r *= s;
         return r;
      }

      // v / s
      template <class D, std::size_t N, class T, class S,
                std::enable_if_t<std::is_arithmetic<S>::value, int> = 0>
      ND_HD ND constexpr D operator/(const VecBase<D, N, T> &a, S s) noexcept
      {
         D r = static_cast<const D &>(a);
         r /= s;
         return r;
      }

      template <class D, std::size_t N, class T, class S, std::enable_if_t<std::is_arithmetic<S>::value, int> = 0>
      ND_HD ND constexpr D operator*(D a, S s) noexcept { return a *= s, a; }
      template <class D, std::size_t N, class T, class S, std::enable_if_t<std::is_arithmetic<S>::value, int> = 0>
      ND_HD ND constexpr D operator*(S s, D a) noexcept { return a *= s, a; }
      template <class D, std::size_t N, class T, class S, std::enable_if_t<std::is_arithmetic<S>::value, int> = 0>
      ND_HD ND constexpr D operator/(D a, S s) noexcept { return a /= s, a; }

      // ============================ helpers & ops ==================================
      template <class D, std::size_t N, class T>
      ND_HD ND constexpr D vmin(D a, const D &b) noexcept
      {
         for (std::size_t i = 0; i < N; ++i)
            a[i] = (a[i] < b[i] ? a[i] : b[i]);
         return a;
      }
      template <class D, std::size_t N, class T>
      ND_HD ND constexpr D vmax(D a, const D &b) noexcept
      {
         for (std::size_t i = 0; i < N; ++i)
            a[i] = (a[i] > b[i] ? a[i] : b[i]);
         return a;
      }
      template <class D, std::size_t N, class T>
      ND_HD ND constexpr D clamp(D x, const D &lo, const D &hi) noexcept { return vmax<D, N, T>(lo, vmin<D, N, T>(x, hi)); }
      template <class D, std::size_t N, class T>
      ND_HD ND constexpr D lerp(const D &a, const D &b, T t) noexcept { return a + (b - a) * t; }

      // 2D
      ND_HD ND constexpr double det2(const Vec2 &a, const Vec2 &b) noexcept { return a.x * b.y - a.y * b.x; }
      ND_HD ND inline double norm2(const Vec2 &v) noexcept { return v.norm(); }
      ND_HD ND constexpr double crossz(const Vec2 &a, const Vec2 &b) noexcept { return det2(a, b); }

      // 3D
      ND_HD ND constexpr double dot(const Vec3 &a, const Vec3 &b) noexcept { return a.dot(b); }
      ND_HD ND constexpr Vec3 cross(const Vec3 &a, const Vec3 &b) noexcept
      {
         return Vec3{
             a.y * b.z - a.z * b.y,
             a.z * b.x - a.x * b.z,
             a.x * b.y - a.y * b.x};
      }

      // 4D (face area-normal from columns A,B,C)
      ND_HD ND constexpr double dot(const Vec4 &a, const Vec4 &b) noexcept { return a.dot(b); }

      ND_HD ND constexpr Vec4 area_normal4(const Vec4 &A, const Vec4 &B, const Vec4 &C) noexcept
      {
         const double n0 = A.y * (B.z * C.w - B.w * C.z) - A.z * (B.y * C.w - B.w * C.y) + A.w * (B.y * C.z - B.z * C.y);
         const double n1 = -(A.x * (B.z * C.w - B.w * C.z) - A.z * (B.x * C.w - B.w * C.x) + A.w * (B.x * C.z - B.z * C.x));
         const double n2 = A.x * (B.y * C.w - B.w * C.y) - A.y * (B.x * C.w - B.w * C.x) + A.w * (B.x * C.y - B.y * C.x);
         const double n3 = -(A.x * (B.y * C.z - B.z * C.y) - A.y * (B.x * C.z - B.z * C.x) + A.z * (B.x * C.y - B.y * C.x));
         return Vec4{n0, n1, n2, n3};
      }

      // streamers
      template <class D, std::size_t N, class T>
      inline std::ostream &operator<<(std::ostream &os, const VecBase<D, N, T> &v)
      {
         os << "Vec (" << N << ") type " << typeid(T).name() << ":\n";

         os << '{';
         for (std::size_t i = 0; i < N; ++i)
         {
            if (i)
               os << ", ";
            os << static_cast<const D &>(v)[i];
         }
         return os << "}T";
      }

      // Free norm that forwards to the member (enables norm(v) style)
      template <class D, std::size_t N, class T>
      ND_HD ND inline T norm(const VecBase<D, N, T> &a) noexcept
      {
         return static_cast<const D &>(a).norm();
      }

      // sanity
      static_assert(std::is_standard_layout<Vec2>::value && std::is_trivially_copyable<Vec2>::value, "Vec2 must be trivial");
      static_assert(std::is_standard_layout<Vec3>::value && std::is_trivially_copyable<Vec3>::value, "Vec3 must be trivial");
      static_assert(std::is_standard_layout<Vec4>::value && std::is_trivially_copyable<Vec4>::value, "Vec4 must be trivial");

   }
} // namespace dstl::math
