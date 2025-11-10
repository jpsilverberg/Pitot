#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <ostream>
#include <type_traits>

// --------------------------------------------------------------------------------------
// ND/HD macros (kept no-ops unless you've defined them elsewhere)
// --------------------------------------------------------------------------------------
#ifndef ND_HD
#define ND_HD
#endif

#if __cplusplus >= 201703L
#ifndef ND
#define ND [[nodiscard]]
#endif
#else
#ifndef ND
#define ND
#endif
#endif

namespace dstl
{
   namespace math
   {

      // ======================================================================================
      // Vec<N,T> : compact value type with raw array storage (trivially copyable).
      // - No exceptions, no RTTI, no virtuals
      // - Elementwise ops, scalar ops, dot/sqnorm/norm, min/max/clamp, lerp
      // - Iterators for std:: algorithms; stream operator for debugging
      // ======================================================================================
      template <std::size_t N, class T = double>
      struct Vec
      {
         T v[N];

         // Ctors
         ND_HD constexpr Vec() noexcept
         {
            for (std::size_t i = 0; i < N; ++i)
               v[i] = T(0);
         }

         ND_HD explicit constexpr Vec(T s) noexcept
         {
            for (std::size_t i = 0; i < N; ++i)
               v[i] = s;
         }

         constexpr Vec(const Vec &) noexcept = default;
         constexpr Vec &operator=(const Vec &) noexcept = default;

// Gate the initializer_list ctor to C++20+ only
#if __cplusplus >= 202002L
         ND_HD constexpr Vec(std::initializer_list<T> il) noexcept
         {
            std::size_t i = 0;
            for (auto x : il)
            {
               if (i < N)
                  v[i++] = x;
            }
            for (; i < N; ++i)
               v[i] = T(0);
         }
#else
         // Variadic ctor: exactly N args, common-type converted to T
         template <class... U,
                   typename = typename std::enable_if<(sizeof...(U) == N) &&
                                                      std::conjunction<std::is_convertible<U, T>...>::value>::type>
         ND_HD constexpr Vec(U... xs) noexcept : v{T(xs)...} {}
#endif

         // Iterators
         ND_HD T *begin() noexcept { return v; }
         ND_HD T *end() noexcept { return v + N; }
         ND_HD const T *begin() const noexcept { return v; }
         ND_HD const T *end() const noexcept { return v + N; }

         // Element access
         ND_HD constexpr T &operator[](std::size_t i) noexcept { return v[i]; }
         ND_HD constexpr const T &operator[](std::size_t i) const noexcept { return v[i]; }

         // Unary
         ND_HD constexpr Vec operator+() const noexcept { return *this; }
         ND_HD constexpr Vec operator-() const noexcept
         {
            Vec r{};
            for (std::size_t i = 0; i < N; ++i)
               r.v[i] = -v[i];
            return r;
         }

         // Compound (elementwise)
         ND_HD constexpr Vec &operator+=(const Vec &b) noexcept
         {
            for (std::size_t i = 0; i < N; ++i)
               v[i] += b.v[i];
            return *this;
         }
         ND_HD constexpr Vec &operator-=(const Vec &b) noexcept
         {
            for (std::size_t i = 0; i < N; ++i)
               v[i] -= b.v[i];
            return *this;
         }
         ND_HD constexpr Vec &operator*=(const Vec &b) noexcept
         {
            for (std::size_t i = 0; i < N; ++i)
               v[i] *= b.v[i];
            return *this;
         }
         ND_HD constexpr Vec &operator/=(const Vec &b) noexcept
         {
            for (std::size_t i = 0; i < N; ++i)
               v[i] /= b.v[i];
            return *this;
         }

         // Compound (scalar)
         template <class S, typename std::enable_if<std::is_arithmetic<S>::value, int>::type = 0>
         ND_HD constexpr Vec &operator*=(S s) noexcept
         {
            for (std::size_t i = 0; i < N; ++i)
               v[i] *= T(s);
            return *this;
         }
         template <class S, typename std::enable_if<std::is_arithmetic<S>::value, int>::type = 0>
         ND_HD constexpr Vec &operator/=(S s) noexcept
         {
            for (std::size_t i = 0; i < N; ++i)
               v[i] /= T(s);
            return *this;
         }

         // Dot/norms
         ND_HD ND constexpr T dot(const Vec &b) const noexcept
         {
            T a = T(0);
            for (std::size_t i = 0; i < N; ++i)
               a += v[i] * b.v[i];
            return a;
         }
         ND_HD ND constexpr T sqnorm() const noexcept { return this->dot(*this); }

         // norm() stays non-constexpr unless we add constexpr sqrt (see below)
         ND_HD ND constexpr T norm() const noexcept
         {
            using std::sqrt;
            return T(sqrt(double(sqnorm())));
         }

         // normalized() can be constexpr if we use constexpr sqrt; otherwise leave runtime
         ND_HD ND constexpr Vec normalized(T eps = T(0)) const noexcept
         {
            T n = norm();
            return (n > eps) ? (*this) / n : Vec{};
         }

         // Utilities
         ND_HD static constexpr Vec zero() noexcept { return Vec{}; }
         ND_HD static constexpr Vec unit(std::size_t i) noexcept
         {
            Vec r{};
            if (i < N)
               r.v[i] = T(1);
            return r;
         }
      };

      // Free operators (elementwise)
      template <std::size_t N, class T>
      ND_HD ND constexpr Vec<N, T> operator+(Vec<N, T> a, const Vec<N, T> &b) noexcept { return a += b, a; }
      template <std::size_t N, class T>
      ND_HD ND constexpr Vec<N, T> operator-(Vec<N, T> a, const Vec<N, T> &b) noexcept { return a -= b, a; }
      template <std::size_t N, class T>
      ND_HD ND constexpr Vec<N, T> operator*(Vec<N, T> a, const Vec<N, T> &b) noexcept { return a *= b, a; }
      template <std::size_t N, class T>
      ND_HD ND constexpr Vec<N, T> operator/(Vec<N, T> a, const Vec<N, T> &b) noexcept { return a /= b, a; }

      // Scalar on either side
      template <std::size_t N, class T, class S, typename std::enable_if<std::is_arithmetic<S>::value, int>::type = 0>
      ND_HD ND constexpr Vec<N, T> operator*(Vec<N, T> a, S s) noexcept { return a *= s, a; }
      template <std::size_t N, class T, class S, typename std::enable_if<std::is_arithmetic<S>::value, int>::type = 0>
      ND_HD ND constexpr Vec<N, T> operator*(S s, Vec<N, T> a) noexcept { return a *= s, a; }
      template <std::size_t N, class T, class S, typename std::enable_if<std::is_arithmetic<S>::value, int>::type = 0>
      ND_HD ND constexpr Vec<N, T> operator/(Vec<N, T> a, S s) noexcept { return a /= s, a; }

      // vmin/vmax/clamp/lerp
      template <std::size_t N, class T>
      ND_HD ND constexpr Vec<N, T> vmin(Vec<N, T> a, const Vec<N, T> &b) noexcept
      {
         for (std::size_t i = 0; i < N; ++i)
            a[i] = std::min(a[i], b[i]);
         return a;
      }
      template <std::size_t N, class T>
      ND_HD ND constexpr Vec<N, T> vmax(Vec<N, T> a, const Vec<N, T> &b) noexcept
      {
         for (std::size_t i = 0; i < N; ++i)
            a[i] = std::max(a[i], b[i]);
         return a;
      }
      template <std::size_t N, class T>
      ND_HD ND constexpr Vec<N, T> clamp(Vec<N, T> x, const Vec<N, T> &lo, const Vec<N, T> &hi) noexcept
      {
         return vmax(lo, vmin(x, hi)); // or use std::clamp elementwise
      }
      template <std::size_t N, class T>
      ND_HD ND constexpr Vec<N, T> lerp(const Vec<N, T> &a, const Vec<N, T> &b, T t) noexcept
      {
         return a + (b - a) * t;
      }

      // Stream
      template <std::size_t N, class T>
      inline std::ostream &operator<<(std::ostream &os, const Vec<N, T> &x)
      {
         os << '(';
         for (std::size_t i = 0; i < N; ++i)
         {
            if (i)
               os << ',';
            os << x.v[i];
         }
         return os << ')';
      }

      // --------------------------------------------------------------------------------------
      // Ergonomic Vec2/Vec3/Vec4 wrappers (same layout; adds .x()/.y()/.z()/.w() accessors)
      // --------------------------------------------------------------------------------------
      struct Vec2 : Vec<2, double>
      {
         using Base = Vec<2, double>;
         using Base::Vec;
         ND_HD double &x() noexcept { return v[0]; }
         ND_HD const double &x() const noexcept { return v[0]; }
         ND_HD double &y() noexcept { return v[1]; }
         ND_HD const double &y() const noexcept { return v[1]; }
      };
      struct Vec3 : Vec<3, double>
      {
         using Base = Vec<3, double>;
         using Base::Vec;
         ND_HD double &x() noexcept { return v[0]; }
         ND_HD const double &x() const noexcept { return v[0]; }
         ND_HD double &y() noexcept { return v[1]; }
         ND_HD const double &y() const noexcept { return v[1]; }
         ND_HD double &z() noexcept { return v[2]; }
         ND_HD const double &z() const noexcept { return v[2]; }
      };
      struct Vec4 : Vec<4, double>
      {
         using Base = Vec<4, double>;
         using Base::Vec;
         ND_HD double &x() noexcept { return v[0]; }
         ND_HD const double &x() const noexcept { return v[0]; }
         ND_HD double &y() noexcept { return v[1]; }
         ND_HD const double &y() const noexcept { return v[1]; }
         ND_HD double &z() noexcept { return v[2]; }
         ND_HD const double &z() const noexcept { return v[2]; }
         ND_HD double &w() noexcept { return v[3]; }
         ND_HD const double &w() const noexcept { return v[3]; }
      };

      // --------------------------------------------------------------------------------------
      // 2D helpers
      // --------------------------------------------------------------------------------------
      ND_HD ND constexpr double dot(const Vec2 &a, const Vec2 &b) noexcept { return a.dot(b); }
      ND_HD ND constexpr double det2(const Vec2 &a, const Vec2 &b) noexcept { return a.v[0] * b.v[1] - a.v[1] * b.v[0]; }
      // z-component of 3D cross where z=det2 (useful for orientation tests)
      ND_HD ND constexpr double crossz(const Vec2 &a, const Vec2 &b) noexcept { return det2(a, b); }

      // --------------------------------------------------------------------------------------
      // 3D helpers
      // --------------------------------------------------------------------------------------
      ND_HD ND constexpr double dot(const Vec3 &a, const Vec3 &b) noexcept { return a.dot(b); }
      ND_HD ND constexpr Vec3 cross(const Vec3 &a, const Vec3 &b) noexcept
      {
         return Vec3{
             a.v[1] * b.v[2] - a.v[2] * b.v[1],
             a.v[2] * b.v[0] - a.v[0] * b.v[2],
             a.v[0] * b.v[1] - a.v[1] * b.v[0]};
      }

      // --------------------------------------------------------------------------------------
      // 4D: area-normal of a 3D face spanned by columns A,B,C (your formula, layout-consistent)
      // --------------------------------------------------------------------------------------
      ND_HD ND constexpr double dot(const Vec4 &a, const Vec4 &b) noexcept { return a.dot(b); }
      ND_HD ND constexpr Vec4 area_normal4(const Vec4 &A, const Vec4 &B, const Vec4 &C) noexcept
      {
         const double n0 = A.v[1] * (B.v[2] * C.v[3] - B.v[3] * C.v[2]) - A.v[2] * (B.v[1] * C.v[3] - B.v[3] * C.v[1]) + A.v[3] * (B.v[1] * C.v[2] - B.v[2] * C.v[1]);
         const double n1 = -(A.v[0] * (B.v[2] * C.v[3] - B.v[3] * C.v[2]) - A.v[2] * (B.v[0] * C.v[3] - B.v[3] * C.v[0]) + A.v[3] * (B.v[0] * C.v[2] - B.v[2] * C.v[0]));
         const double n2 = A.v[0] * (B.v[1] * C.v[3] - B.v[3] * C.v[1]) - A.v[1] * (B.v[0] * C.v[3] - B.v[3] * C.v[0]) + A.v[3] * (B.v[0] * C.v[1] - B.v[1] * C.v[0]);
         const double n3 = -(A.v[0] * (B.v[1] * C.v[2] - B.v[2] * C.v[1]) - A.v[1] * (B.v[0] * C.v[2] - B.v[2] * C.v[0]) + A.v[2] * (B.v[0] * C.v[1] - B.v[1] * C.v[0]));
         return Vec4{n0, n1, n2, n3};
      }

      // --------------------------------------------------------------------------------------
      // Traits & sanity checks
      // --------------------------------------------------------------------------------------
      static_assert(std::is_trivially_copyable<Vec<2, double>>::value, "Vec must be trivial");
      static_assert(std::is_trivially_copyable<Vec<3, double>>::value, "Vec must be trivial");
      static_assert(std::is_trivially_copyable<Vec<4, double>>::value, "Vec must be trivial");

      // Deduce N and T from argument list: Vec v{1,2,3} -> Vec<3,int>
      template <class... U>
      Vec(U...) -> Vec<sizeof...(U), typename std::common_type<U...>::type>;

      // template <class... U>
      // Vec(U...) -> Vec<sizeof...(U), std::common_type_t<U..., double>>;

   }
} // namespace dstl::math
