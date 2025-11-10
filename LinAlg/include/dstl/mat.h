// include/dstl/math/mat.h  (C++17)
// Matrix with CRTP base, row-major storage, zero-overhead interop with Vec.
#pragma once
#include <array>
#include <cstddef>
#include <type_traits>
#include <algorithm>
#include <initializer_list>
#include <iomanip>
#include <ostream>
#include <typeinfo>
#include <cmath>

#include <dstl/dprint.h>
#include <dstl/Numbers.h>

#ifndef ND_HD
#define ND_HD
#endif
#ifndef ND
#define ND [[nodiscard]]
#endif

#include "vec.h" // uses dstl::math VecBase/Vec/Vec2/Vec3/Vec4

namespace dstl
{
   namespace linalg
   {

      // Forward declaration so MatBase can name Mat in member functions like transpose()
      template <std::size_t M, std::size_t N, class T>
      struct Mat;

      // =============================== MatBase =====================================
      template <class Derived, std::size_t M, std::size_t N, class T>
      struct MatBase
      {
         using value_type = T;
         static constexpr std::size_t kRows = M;
         static constexpr std::size_t kCols = N;

         // storage hooks the derived must provide
         ND_HD constexpr T *data() noexcept { return static_cast<Derived *>(this)->data_impl(); }
         ND_HD constexpr const T *data() const noexcept { return static_cast<const Derived *>(this)->data_impl(); }

         // element access (row-major)
         ND_HD constexpr T &operator()(std::size_t r, std::size_t c) noexcept { return data()[r * N + c]; }
         ND_HD constexpr const T &operator()(std::size_t r, std::size_t c) const noexcept { return data()[r * N + c]; }

         // element access / iterators
         ND_HD constexpr T &operator[](std::size_t i) noexcept
         {
            assert(i < M * N);
            return data()[i];
         }
         ND_HD constexpr const T &operator[](std::size_t i) const noexcept
         {
            assert(i < M * N);
            return data()[i];
         }

         // iterators
         ND_HD constexpr T *begin() noexcept { return data(); }
         ND_HD constexpr T *end() noexcept { return data() + M * N; }
         ND_HD constexpr const T *begin() const noexcept { return data(); }
         ND_HD constexpr const T *end() const noexcept { return data() + M * N; }

         // unary
         ND_HD constexpr Derived operator+() const noexcept { return derived(); }
         ND_HD constexpr Derived operator-() const noexcept
         {
            Derived r{};
            for (std::size_t i = 0; i < M * N; ++i)
               r.data()[i] = -data()[i];
            return r;
         }

         // elementwise compound
         ND_HD constexpr Derived &operator+=(const Derived &b) noexcept
         {
            for (std::size_t i = 0; i < M * N; ++i)
               data()[i] += b.data()[i];
            return derived();
         }
         ND_HD constexpr Derived &operator-=(const Derived &b) noexcept
         {
            for (std::size_t i = 0; i < M * N; ++i)
               data()[i] -= b.data()[i];
            return derived();
         }
         ND_HD constexpr Derived &operator*=(const Derived &b) noexcept
         { // Hadamard
            for (std::size_t i = 0; i < M * N; ++i)
               data()[i] *= b.data()[i];
            return derived();
         }
         ND_HD constexpr Derived &operator/=(const Derived &b) noexcept
         {
            for (std::size_t i = 0; i < M * N; ++i)
               data()[i] /= b.data()[i];
            return derived();
         }

         // scalar compound
         template <class S, std::enable_if_t<std::is_arithmetic<S>::value, int> = 0>
         ND_HD constexpr Derived &operator*=(S s) noexcept
         {
            for (auto &x : *this)
               x = T(x * s);
            return derived();
         }
         template <class S, std::enable_if_t<std::is_arithmetic<S>::value, int> = 0>
         ND_HD constexpr Derived &operator/=(S s) noexcept
         {
            for (auto &x : *this)
               x = T(x / s);
            return derived();
         }

         // row/col accessors using your Vec types
         ND_HD ND constexpr Vec<N, T> row(std::size_t r) const noexcept
         {
            Vec<N, T> v{};
            for (std::size_t j = 0; j < N; ++j)
               v[j] = (*this)(r, j);
            return v;
         }
         ND_HD ND constexpr Vec<M, T> col(std::size_t c) const noexcept
         {
            Vec<M, T> v{};
            for (std::size_t i = 0; i < M; ++i)
               v[i] = (*this)(i, c);
            return v;
         }
         ND_HD constexpr void set_row(std::size_t r, const Vec<N, T> &v) noexcept
         {
            for (std::size_t j = 0; j < N; ++j)
               (*this)(r, j) = v[j];
         }
         ND_HD constexpr void set_col(std::size_t c, const Vec<M, T> &v) noexcept
         {
            for (std::size_t i = 0; i < M; ++i)
               (*this)(i, c) = v[i];
         }

         // trace (square only)
         ND_HD ND constexpr T trace() const noexcept
         {
            static_assert(M == N, "trace requires square matrix");
            T s = T(0);
            for (std::size_t i = 0; i < N; ++i)
               s += (*this)(i, i);
            return s;
         }

         // transpose
         ND_HD ND constexpr auto transpose() const noexcept
         {
            Mat<N, M, T> r{};
            for (std::size_t i = 0; i < M; ++i)
               for (std::size_t j = 0; j < N; ++j)
                  r(j, i) = (*this)(i, j);
            return r;
         }

         // identity & zeros
         ND_HD static constexpr Derived zeros() noexcept { return Derived{}; }
         ND_HD static constexpr Derived eye() noexcept
         {
            static_assert(M == N, "eye() requires square matrix");
            Derived I{};
            for (std::size_t i = 0; i < N; ++i)
               I(i, i) = T(1);
            return I;
         }

         // builders from rows/cols
         ND_HD static constexpr Derived from_rows(const Vec<N, T> (&rows)[M]) noexcept
         {
            Derived A{};
            for (std::size_t i = 0; i < M; ++i)
               for (std::size_t j = 0; j < N; ++j)
                  A(i, j) = rows[i][j];
            return A;
         }
         ND_HD static constexpr Derived from_cols(const Vec<M, T> (&cols)[N]) noexcept
         {
            Derived A{};
            for (std::size_t j = 0; j < N; ++j)
               for (std::size_t i = 0; i < M; ++i)
                  A(i, j) = cols[j][i];
            return A;
         }

         void fill(T s) noexcept
         {
            for (std::size_t i = 0; i < M * N; ++i)
               data()[i] = s;
         }


      private:
         ND_HD constexpr Derived &derived() noexcept { return *static_cast<Derived *>(this); }
         ND_HD constexpr const Derived &derived() const noexcept { return *static_cast<const Derived *>(this); }
      };

      // ============================= Mat storage ===================================
      template <std::size_t M, std::size_t N, class T = double>
      struct Mat : MatBase<Mat<M, N, T>, M, N, T>
      {
         using Base = MatBase<Mat<M, N, T>, M, N, T>;
         T a[M * N];

         ND_HD constexpr Mat() noexcept
         {
            for (std::size_t i = 0; i < M * N; ++i)
               a[i] = T(0);
         }
         ND_HD explicit constexpr Mat(T s) noexcept
         {
            for (std::size_t i = 0; i < M * N; ++i)
               a[i] = s;
         }

         // C++17 variadic ctor: exactly M*N args
         template <class... U,
                   std::enable_if_t<(sizeof...(U) == M * N) && std::conjunction<std::is_convertible<U, T>...>::value, int> = 0>
         ND_HD constexpr Mat(U... xs) noexcept : a{T(xs)...} {}

         // (Optional) initializer_list ctor (not constexpr pre-C++20)
#if __cplusplus >= 202002L
         ND_HD constexpr Mat(std::initializer_list<T> il) noexcept : a{}
         {
            std::size_t i = 0;
            for (auto x : il)
            {
               if (i < M * N)
                  a[i++] = x;
            }
         }
#endif

         // storage hooks
         ND_HD constexpr T *data_impl() noexcept { return a; }
         ND_HD constexpr const T *data_impl() const noexcept { return a; }
      };

      // ============================= Free operators ================================
      // elementwise
      template <class D, std::size_t M, std::size_t N, class T>
      ND_HD ND constexpr D operator+(D a, const MatBase<D, M, N, T> &b) noexcept
      {
         return a += static_cast<const D &>(b), a;
      }
      template <class D, std::size_t M, std::size_t N, class T>
      ND_HD ND constexpr D operator-(D a, const MatBase<D, M, N, T> &b) noexcept
      {
         return a -= static_cast<const D &>(b), a;
      }
      template <class D, std::size_t M, std::size_t N, class T>
      ND_HD ND constexpr D hadamard(D a, const MatBase<D, M, N, T> &b) noexcept
      {
         return a *= static_cast<const D &>(b), a;
      }

      // scalar
      template <class D, std::size_t M, std::size_t N, class T, class S,
                std::enable_if_t<std::is_arithmetic<S>::value, int> = 0>
      ND_HD ND constexpr D operator*(const MatBase<D, M, N, T> &A, S s) noexcept
      {
         D r = static_cast<const D &>(A);
         r *= s;
         return r;
      }
      template <class D, std::size_t M, std::size_t N, class T, class S,
                std::enable_if_t<std::is_arithmetic<S>::value, int> = 0>
      ND_HD ND constexpr D operator*(S s, const MatBase<D, M, N, T> &A) noexcept
      {
         D r = static_cast<const D &>(A);
         r *= s;
         return r;
      }
      template <class D, std::size_t M, std::size_t N, class T, class S,
                std::enable_if_t<std::is_arithmetic<S>::value, int> = 0>
      ND_HD ND constexpr D operator/(const MatBase<D, M, N, T> &A, S s) noexcept
      {
         D r = static_cast<const D &>(A);
         r /= s;
         return r;
      }

      // mat × vec
      template <std::size_t M, std::size_t N, class T>
      ND_HD ND constexpr Vec<M, T> operator*(const MatBase<Mat<M, N, T>, M, N, T> &A, const Vec<N, T> &x) noexcept
      {
         Vec<M, T> y{};
         for (std::size_t i = 0; i < M; ++i)
         {
            T s = T(0);
            for (std::size_t j = 0; j < N; ++j)
               s += static_cast<const Mat<M, N, T> &>(A)(i, j) * x[j];
            y[i] = s;
         }
         return y;
      }

      // vec × mat
      template <std::size_t M, std::size_t N, class T>
      ND_HD ND constexpr Vec<N, T> operator*(const Vec<M, T> &x, const MatBase<Mat<M, N, T>, M, N, T> &A) noexcept
      {
         Vec<N, T> y{};
         for (std::size_t j = 0; j < N; ++j)
         {
            T s = T(0);
            for (std::size_t i = 0; i < M; ++i)
               s += x[i] * static_cast<const Mat<M, N, T> &>(A)(i, j);
            y[j] = s;
         }
         return y;
      }

      // mat × mat
      template <std::size_t M, std::size_t K, std::size_t N, class T>
      ND_HD ND constexpr Mat<M, N, T> operator*(const MatBase<Mat<M, K, T>, M, K, T> &A,
                                                const MatBase<Mat<K, N, T>, K, N, T> &B) noexcept
      {
         Mat<M, N, T> C{};
         const auto &AA = static_cast<const Mat<M, K, T> &>(A);
         const auto &BB = static_cast<const Mat<K, N, T> &>(B);
         for (std::size_t i = 0; i < M; ++i)
         {
            for (std::size_t j = 0; j < N; ++j)
            {
               T s = T(0);
               for (std::size_t k = 0; k < K; ++k)
                  s += AA(i, k) * BB(k, j);
               C(i, j) = s;
            }
         }
         return C;
      }

      // ========================== Determinant / Inverse ============================

      // 2×2 det/inv (fast paths)
      template <class T>
      ND_HD ND constexpr T det(const MatBase<Mat<2, 2, T>, 2, 2, T> &A_) noexcept
      {
         const auto &A = static_cast<const Mat<2, 2, T> &>(A_);
         return A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
      }
      template <class T>
      ND_HD ND inline Mat<2, 2, T> inv(const MatBase<Mat<2, 2, T>, 2, 2, T> &A_) noexcept
      {
         const auto &A = static_cast<const Mat<2, 2, T> &>(A_);
         T d = det(A);
         // (Optional: assert |d|>eps)
         return (T(1) / d) * Mat<2, 2, T>{A(1, 1), -A(0, 1),
                                          -A(1, 0), A(0, 0)};
      }

      // 3×3 det/inv (cofactor form)
      template <class T>
      ND_HD ND constexpr T det(const MatBase<Mat<3, 3, T>, 3, 3, T> &A_) noexcept
      {
         const auto &A = static_cast<const Mat<3, 3, T> &>(A_);
         return A(0, 0) * (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)) -
                A(0, 1) * (A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0)) +
                A(0, 2) * (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0));
      }
      template <class T>
      ND_HD ND inline Mat<3, 3, T> inv(const MatBase<Mat<3, 3, T>, 3, 3, T> &A_) noexcept
      {
         const auto &A = static_cast<const Mat<3, 3, T> &>(A_);
         const T d = det(A);
         // (Optional: assert |d|>eps)
         Mat<3, 3, T> C{};
         C(0, 0) = (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1));
         C(0, 1) = -(A(0, 1) * A(2, 2) - A(0, 2) * A(2, 1));
         C(0, 2) = (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1));
         C(1, 0) = -(A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0));
         C(1, 1) = (A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0));
         C(1, 2) = -(A(0, 0) * A(1, 2) - A(0, 2) * A(1, 0));
         C(2, 0) = (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0));
         C(2, 1) = -(A(0, 0) * A(2, 1) - A(0, 1) * A(2, 0));
         C(2, 2) = (A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0));
         return (T(1) / d) * C.transpose();
      }

      // Generic square det/inv via Gauss–Jordan with partial pivoting (small N)
      template <std::size_t N, class T>
      ND_HD ND inline T det(const MatBase<Mat<N, N, T>, N, N, T> &A_) noexcept
      {
         Mat<N, N, T> U = static_cast<const Mat<N, N, T> &>(A_);
         T d = T(1);
         for (std::size_t k = 0; k < N; ++k)
         {
            // pivot
            std::size_t p = k;
            T amax = std::abs(U(k, k));
            for (std::size_t i = k + 1; i < N; ++i)
            {
               T v = std::abs(U(i, k));
               if (v > amax)
               {
                  amax = v;
                  p = i;
               }
            }
            if (amax == T(0))
               return T(0);
            if (p != k)
            { // swap rows
               for (std::size_t j = k; j < N; ++j)
                  std::swap(U(k, j), U(p, j));
               d = -d;
            }
            T piv = U(k, k);
            d *= piv;
            // eliminate below
            for (std::size_t i = k + 1; i < N; ++i)
            {
               T f = U(i, k) / piv;
               for (std::size_t j = k + 1; j < N; ++j)
                  U(i, j) -= f * U(k, j);
            }
         }
         return d;
      }

      template <std::size_t N, class T>
      ND_HD ND inline Mat<N, N, T> inv(const MatBase<Mat<N, N, T>, N, N, T> &A_) noexcept
      {
         const Mat<N, N, T> &A = static_cast<const Mat<N, N, T> &>(A_);
         Mat<N, N, T> U = A;
         Mat<N, N, T> I = Mat<N, N, T>::eye();

         // Gauss–Jordan
         for (std::size_t k = 0; k < N; ++k)
         {
            // pivot
            std::size_t p = k;
            T amax = std::abs(U(k, k));
            for (std::size_t i = k + 1; i < N; ++i)
            {
               T v = std::abs(U(i, k));
               if (v > amax)
               {
                  amax = v;
                  p = i;
               }
            }
            // Singular?
            // (Optional) if (amax <= std::numeric_limits<T>::epsilon()) { /* handle */ }
            if (p != k)
            {
               for (std::size_t j = 0; j < N; ++j)
               {
                  std::swap(U(k, j), U(p, j));
                  std::swap(I(k, j), I(p, j));
               }
            }
            const T piv = U(k, k);
            // normalize row k
            const T invp = T(1) / piv;
            for (std::size_t j = 0; j < N; ++j)
            {
               U(k, j) *= invp;
               I(k, j) *= invp;
            }
            // eliminate others
            for (std::size_t i = 0; i < N; ++i)
               if (i != k)
               {
                  const T f = U(i, k);
                  if (f != T(0))
                  {
                     for (std::size_t j = 0; j < N; ++j)
                     {
                        U(i, j) -= f * U(k, j);
                        I(i, j) -= f * I(k, j);
                     }
                  }
               }
         }
         return I;
      }

      // ================================ LU (Doolittle, partial pivoting) ================================
      // Left-looking Doolittle with partial pivoting, Accum used on the two reductions.
      template <std::size_t N, class T>
      struct LU_Left
      {
         Mat<N, N, T> lu;      // Packed L (strictly below diag) and U (on/above)
         Vec<N, uint32_t> piv; // Row permutation P (applied as P*A = L*U)
         int sign = 1;         // 0 => singular

         // Stream
         std::string ToString() const
         {
            std::string o_Str = dstl::print("LU Decomposition >>\n");
            std::stringstream l_ss;
            l_ss << lu;
            o_Str += l_ss.str();
            o_Str += dstl::print("\nPivots %s: ", sign > 0 ? "+" : "-");
            l_ss.str("");
            l_ss << piv;
            o_Str += l_ss.str();
            return o_Str;
         }

         template <class Acc = Accum>
         ND_HD ND static inline LU_Left factor(const MatBase<Mat<N, N, T>, N, N, T> &A_) noexcept
         {
            LU_Left d{};
            d.lu = static_cast<const Mat<N, N, T> &>(A_);
            for (std::size_t i = 0; i < N; ++i)
               d.piv[i] = i;
            d.sign = 1;

            for (std::size_t k = 0; k < N; ++k)
            {
               // --- Pivot search on column k (in current panel) ---
               std::size_t p = k;
               T amax = std::abs(d.lu(k, k));
               for (std::size_t i = k + 1; i < N; ++i)
               {
                  T v = std::abs(d.lu(i, k));
                  if (v > amax)
                  {
                     amax = v;
                     p = i;
                  }
               }
               if (amax == T(0))
               {
                  d.sign = 0;
                  return d;
               } // singular

               if (p != k)
               {
                  for (std::size_t j = 0; j < N; ++j)
                     std::swap(d.lu(k, j), d.lu(p, j));
                  std::swap(d.piv[k], d.piv[p]);
                  d.sign = -d.sign;
               }

               // --- Compute U row k: U(k,j) for j=k..N-1 (left-looking dot) ---
               for (std::size_t j = k; j < N; ++j)
               {
                  Acc s = Acc(d.lu(k, j));
                  for (std::size_t q = 0; q < k; ++q)
                     s -= Acc(d.lu(k, q) * d.lu(q, j)); // L(k,q)*U(q,j)
                  d.lu(k, j) = s.Double();
               }

               const T ukk = d.lu(k, k);
               if (ukk == T(0))
               {
                  d.sign = 0;
                  return d;
               }

               // --- Compute L column k: L(i,k) for i=k+1..N-1 (left-looking dot) ---
               for (std::size_t i = k + 1; i < N; ++i)
               {
                  Acc s = Acc(d.lu(i, k));
                  for (std::size_t q = 0; q < k; ++q)
                     s -= Acc(d.lu(i, q) * d.lu(q, k)); // L(i,q)*U(q,k)
                  d.lu(i, k) = s.Double() / ukk;
               }
               // No trailing rank-1 update needed here—the left-looking dots already account for it.
            }
            return d;
         }

         // Reuse your Accum-enabled solve from before (works unchanged).
         template <class DV, class Acc = Accum>
         ND_HD ND inline Vec<N, T> solve(const VecBase<DV, N, T> &b) const noexcept
         {
            Vec<N, T> x{};
            for (std::size_t i = 0; i < N; ++i)
               x[i] = static_cast<const DV &>(b)[piv[i]];
            if (sign == 0)
               return Vec<N, T>{};

            // Forward: L y = P b
            for (std::size_t i = 0; i < N; ++i)
            {
               Acc s = Acc(x[i]);
               for (std::size_t j = 0; j < i; ++j)
                  s -= Acc(lu(i, j) * x[j]);
               x[i] = s.Double();
            }
            // Backward: U x = y
            for (std::size_t i = N; i-- > 0;)
            {
               Acc s = Acc(x[i]);
               for (std::size_t j = i + 1; j < N; ++j)
                  s -= Acc(lu(i, j) * x[j]);
               x[i] = s.Double() / lu(i, i);
            }
            return x;
         }

         ND_HD ND inline T det() const noexcept
         {
            if (sign == 0)
               return T(0);
            T d = (sign > 0 ? T(1) : T(-1));
            for (std::size_t i = 0; i < N; ++i)
               d *= lu(i, i);
            return d;
         }

         ND_HD ND inline Mat<N, N, T> inverse() const noexcept
         {
            Mat<N, N, T> invA{};
            if (sign == 0)
               return invA;
            for (std::size_t c = 0; c < N; ++c)
            {
               Vec<N, T> e{};
               e[c] = T(1);
               auto xc = solve<Vec<N, T>, Accum>(e);
               for (std::size_t i = 0; i < N; ++i)
                  invA(i, c) = xc[i];
            }
            return invA;
         }
      };

      template <std::size_t N, class T>
      struct LU
      {
         Mat<N, N, T> lu;      // Packed: i>j => L_ij, i<=j => U_ij; L has unit diagonal (not stored)
         Vec<N, uint32_t> piv; // Row permutation P (applied as P*A = L*U)
         int sign = 1;         // +1/-1 for even/odd permutations, 0 => singular

         // Stream
         std::string ToString() const
         {
            std::string o_Str = dstl::print("LU Decomposition >>\n");
            std::stringstream l_ss;
            l_ss << lu;
            o_Str += l_ss.str();
            o_Str += dstl::print("\nPivots %s: ", sign > 0 ? "+" : "-");
            l_ss.str("");
            l_ss << piv;
            o_Str += l_ss.str();
            return o_Str;
         }

         ND_HD ND static inline LU factor(const MatBase<Mat<N, N, T>, N, N, T> &A_) noexcept
         {
            LU d{};
            d.lu = static_cast<const Mat<N, N, T> &>(A_);
            for (std::size_t i = 0; i < N; ++i)
               d.piv[i] = i;
            d.sign = 1;

            for (std::size_t k = 0; k < N; ++k)
            {
               // Pivot selection (max |entry|)
               std::size_t p = k;
               T amax = std::abs(d.lu(k, k));
               for (std::size_t i = k + 1; i < N; ++i)
               {
                  T v = std::abs(d.lu(i, k));
                  if (v > amax)
                  {
                     amax = v;
                     p = i;
                  }
               }
               // Singular?
               if (amax == T(0))
               {
                  d.sign = 0;
                  return d;
               }

               // Row swap if needed
               if (p != k)
               {
                  for (std::size_t j = 0; j < N; ++j)
                     std::swap(d.lu(k, j), d.lu(p, j));
                  std::swap(d.piv[k], d.piv[p]);
                  d.sign = -d.sign;
               }

               const T piv = d.lu(k, k);

               // Compute multipliers for L (below diagonal)
               for (std::size_t i = k + 1; i < N; ++i)
               {
                  d.lu(i, k) /= piv;
               }
               // Schur complement update: A(i,j) -= L(i,k) * U(k,j)
               for (std::size_t i = k + 1; i < N; ++i)
               {
                  const T lik = d.lu(i, k);
                  if (lik != T(0))
                  {
                     for (std::size_t j = k + 1; j < N; ++j)
                     {
                        d.lu(i, j) -= lik * d.lu(k, j);
                     }
                  }
               }
            }
            return d;
         }

         // Solve A x = b using PA = LU  =>  L U x = P b
         template <class DV>
         ND_HD ND inline Vec<N, T> solve(const VecBase<DV, N, T> &b) const noexcept
         {
            Vec<N, T> x{}; // will hold Pb -> y -> x

            // Permute: x = P b
            for (std::size_t i = 0; i < N; ++i)
               x[i] = static_cast<const DV &>(b)[piv[i]];

            if (sign == 0)
            {
               return Vec<N, T>{}; // singular: return zero (caller can check det()==0)
            }

            // Forward solve: L y = P b  (unit diagonal L)
            for (std::size_t i = 0; i < N; ++i)
            {
               Accum s = Accum(x[i]);
               for (std::size_t j = 0; j < i; ++j)
                  s -= Accum(lu(i, j) * x[j]);
               x[i] = s.Double(); // y_i
            }

            // Backward solve: U x = y
            for (std::size_t i = N; i-- > 0;)
            {
               Accum s = Accum(x[i]);
               for (std::size_t j = i + 1; j < N; ++j)
                  s -= Accum(lu(i, j) * x[j]);
               // s /= lu(i, i);
               x[i] = s.Double() / lu(i, i);
            }
            return x;
         }

         // Determinant: det(A) = sign * prod(diag(U)); sign==0 => det=0
         ND_HD ND inline T det() const noexcept
         {
            if (sign == 0)
               return T(0);
            T d = (sign > 0 ? T(1) : T(-1));
            for (std::size_t i = 0; i < N; ++i)
               d *= lu(i, i);
            return d;
         }

         // Inverse via column solves
         ND_HD ND inline Mat<N, N, T> inverse() const noexcept
         {
            Mat<N, N, T> invA{};
            if (sign == 0)
            {
               // Singular: return zero (caller can check by det()==0)
               return invA;
            }
            for (std::size_t col = 0; col < N; ++col)
            {
               Vec<N, T> e{};
               e[col] = T(1);
               auto x = solve(e);
               for (std::size_t i = 0; i < N; ++i)
                  invA(i, col) = x[i];
            }
            return invA;
         }
      };

      // Free convenience wrappers
      template <std::size_t N, class T>
      ND_HD ND inline LU<N, T> lu_factor(const MatBase<Mat<N, N, T>, N, N, T> &A) noexcept
      {
         return LU<N, T>::factor(A);
      }

      template <std::size_t N, class T>
      ND_HD ND inline Vec<N, T> lu_solve(const MatBase<Mat<N, N, T>, N, N, T> &A, const Vec<N, T> &b) noexcept
      {
         auto lu = lu_factor<A.kRows, T>(A); // not usable: A is a type; write explicitly below
         // Simpler: factor then solve
         auto LUF = LU<N, T>::factor(A);
         return LUF.solve(b);
      }

      template <std::size_t N, class T>
      ND_HD ND inline T lu_det(const MatBase<Mat<N, N, T>, N, N, T> &A) noexcept
      {
         auto LUF = LU<N, T>::factor(A);
         return LUF.det();
      }

      template <std::size_t N, class T>
      ND_HD ND inline Mat<N, N, T> lu_inverse(const MatBase<Mat<N, N, T>, N, N, T> &A) noexcept
      {
         auto LUF = LU<N, T>::factor(A);
         return LUF.inverse();
      }

      // ================================ Utilities ==================================

      // Outer product: u (M) × v^T (N) → M×N
      template <std::size_t M, std::size_t N, class T>
      ND_HD ND constexpr Mat<M, N, T> outer(const Vec<M, T> &u, const Vec<N, T> &v) noexcept
      {
         Mat<M, N, T> A{};
         for (std::size_t i = 0; i < M; ++i)
            for (std::size_t j = 0; j < N; ++j)
               A(i, j) = u[i] * v[j];
         return A;
      }

      // Stream
      template <class D, std::size_t M, std::size_t N, class T>
      inline std::ostream &operator<<(std::ostream &os, const MatBase<D, M, N, T> &A)
      {
         os << "Matrix (" << M << ", " << N << ") type " << typeid(T).name() << ":\n";

         // First pass: compute column widths
         std::size_t col_widths[N] = {};
         for (std::size_t j = 0; j < N; ++j)
         {
            for (std::size_t i = 0; i < M; ++i)
            {
               std::ostringstream temp;
               temp << static_cast<const D &>(A)(i, j);
               col_widths[j] = std::max(col_widths[j], temp.str().length());
            }
         }

         // Second pass: print with alignment
         for (std::size_t i = 0; i < M; ++i)
         {
            os << (i ? "\n[ " : "[ ");
            for (std::size_t j = 0; j < N; ++j)
            {
               if (j)
                  os << ", ";
               os << std::setw(col_widths[j]) << static_cast<const D &>(A)(i, j);
            }
            os << " ]";
         }
         return os;
      }
   }
} // namespace dstl::math
