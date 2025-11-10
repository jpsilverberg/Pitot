#pragma once
#include <array>
#include <cmath>
#include <cassert>
#include <utility>
#include <initializer_list>

// detection: el.A0?
template <class T>
struct has_A0
{
private:
   template <class U>
   static auto check(int) -> decltype((void)std::declval<U &>().A0, std::true_type{});
   template <class>
   static std::false_type check(...);

public:
   static constexpr bool value = decltype(check<T>(0))::value;
};

template <class T>
struct has_B0
{
private:
   template <class U>
   static auto check(int) -> decltype((void)std::declval<U &>().B0, std::true_type{});
   template <class>
   static std::false_type check(...);

public:
   static constexpr bool value = decltype(check<T>(0))::value;
};

// get_A0(el): returns el.A0 if it exists, else 0.0
template <class El>
inline typename std::enable_if<has_A0<El>::value, double>::type
get_A0(El &el) { return el.A0; }

template <class El>
inline typename std::enable_if<!has_A0<El>::value, double>::type
get_A0(El &) { return 0.0; }

// get_B0(el): returns el.B0 if it exists, else 0.0
template <class El>
inline typename std::enable_if<has_B0<El>::value, double>::type
get_B0(El &el) { return el.B0; }

template <class El>
inline typename std::enable_if<!has_B0<El>::value, double>::type
get_B0(El &) { return 0.0; }

namespace quad
{

   // -----------------------------
   // Small helpers / math
   // -----------------------------
   struct Mat2
   {
      double a11, a12, a21, a22; // [ [a11 a12], [a21 a22] ]
      static Mat2 rotAlign(const std::array<double, 2> &p0,
                           const std::array<double, 2> &p1)
      {
         double dx = p1[0] - p0[0], dy = p1[1] - p0[1];
         double n = std::sqrt(dx * dx + dy * dy);
         double c = 1.0, s = 0.0;
         if (n > 1e-16)
         {
            c = dx / n;
            s = dy / n;
         }
         return {c, -s, s, c};
      }
      std::array<double, 2> Tmul(const std::array<double, 2> &v) const
      {
         return {a11 * v[0] + a21 * v[1], a12 * v[0] + a22 * v[1]}; // R^T v
      }
      std::array<double, 2> mul(const std::array<double, 2> &v) const
      {
         return {a11 * v[0] + a12 * v[1], a21 * v[0] + a22 * v[1]}; // R v
      }
   };

   // unit-normal pair in AB frame (x-axis along AB)
   inline std::pair<Mat2, double> normVec2D(const std::array<double, 2> &a,
                                            const std::array<double, 2> &b)
   {
      Mat2 R = Mat2::rotAlign(a, b);
      double dx = b[0] - a[0], dy = b[1] - a[1];
      return {R, std::sqrt(dx * dx + dy * dy)};
   }

   // -----------------------------
   // Gauss–Legendre nodes/weights
   // (constexpr tables for n = 2..8)
   // -----------------------------
   template <int N>
   struct GL;

   template <>
   struct GL<2>
   {
      static constexpr std::array<double, 2> x{{-0.57735026918962576451, 0.57735026918962576451}};
      static constexpr std::array<double, 2> w{{1.0, 1.0}};
   };

   template <>
   struct GL<3>
   {
      static constexpr std::array<double, 3> x{{-0.77459666924148337704, 0.0, 0.77459666924148337704}};
      static constexpr std::array<double, 3> w{{0.55555555555555555556, 0.88888888888888888889, 0.55555555555555555556}};
   };

   template <>
   struct GL<4>
   {
      static constexpr std::array<double, 4> x{{-0.86113631159405257522, -0.33998104358485626480,
                                                0.33998104358485626480, 0.86113631159405257522}};
      static constexpr std::array<double, 4> w{{0.34785484513745385737, 0.65214515486254614263,
                                                0.65214515486254614263, 0.34785484513745385737}};
   };

   template <>
   struct GL<5>
   {
      static constexpr std::array<double, 5> x{{-0.90617984593866399280, -0.53846931010568309104, 0.0,
                                                0.53846931010568309104, 0.90617984593866399280}};
      static constexpr std::array<double, 5> w{{0.23692688505618908751, 0.47862867049936646804, 0.56888888888888888889,
                                                0.47862867049936646804, 0.23692688505618908751}};
   };

   template <>
   struct GL<6>
   {
      static constexpr std::array<double, 6> x{{-0.93246951420315202781, -0.66120938646626451366, -0.23861918608319690863,
                                                0.23861918608319690863, 0.66120938646626451366, 0.93246951420315202781}};
      static constexpr std::array<double, 6> w{{0.17132449237917034504, 0.36076157304813860757, 0.46791393457269104739,
                                                0.46791393457269104739, 0.36076157304813860757, 0.17132449237917034504}};
   };

   template <>
   struct GL<7>
   {
      static constexpr std::array<double, 7> x{{-0.94910791234275852453, -0.74153118559939443986, -0.40584515137739716691,
                                                0.0,
                                                0.40584515137739716691, 0.74153118559939443986, 0.94910791234275852453}};
      static constexpr std::array<double, 7> w{{0.12948496616886969327, 0.27970539148927666790, 0.38183005050511894495,
                                                0.41795918367346938776,
                                                0.38183005050511894495, 0.27970539148927666790, 0.12948496616886969327}};
   };

   template <>
   struct GL<8>
   {
      static constexpr std::array<double, 8> x{{-0.96028985649753623168, -0.79666647741362673959, -0.52553240991632898582, -0.18343464249564980494,
                                                0.18343464249564980494, 0.52553240991632898582, 0.79666647741362673959, 0.96028985649753623168}};
      static constexpr std::array<double, 8> w{{0.10122853629037625915, 0.22238103445337447054, 0.31370664587788728734, 0.36268378337836198297,
                                                0.36268378337836198297, 0.31370664587788728734, 0.22238103445337447054, 0.10122853629037625915}};
   };

   // -----------------------------
   // 1D Gauss–Legendre (generic N)
   // -----------------------------
   template <int N, class F, class El>
   inline double gl_quad_1d(F &&f, El &el, double a, double b)
   {
      static_assert(N >= 2 && N <= 8, "Supported N in [2..8]");
      const auto &xs = GL<N>::x;
      const auto &ws = GL<N>::w;
      const double m = 0.5 * (b - a);
      const double c = 0.5 * (b + a);
      double acc = 0.0;
      for (size_t i = 0; i < xs.size(); ++i)
      {
         const double x = m * xs[i] + c;
         acc += ws[i] * f(el, x);
      }
      return m * acc;
   }

   // -----------------------------
   // Adaptive 1D using fixed N
   // -----------------------------
   struct AdaptTol
   {
      double rtol = 3e-8;
      double atol = 1e-12;
   };

   template <int N, class F, class El>
   double quad_adapt(F &&f, El &el, double a, double b, const AdaptTol tol = {})
   {
      const auto whole = gl_quad_1d<N>(f, el, a, b);
      const double mid = 0.5 * (a + b);
      const auto left = gl_quad_1d<N>(f, el, a, mid);
      const auto right = gl_quad_1d<N>(f, el, mid, b);
      const auto better = left + right;
      const double scale = std::fmax(std::fabs(better), 1.0);
      if (std::fabs(better - whole) <= std::fmax(tol.rtol * scale, tol.atol))
      {
         return better;
      }
      // Note: caller's El can track recursion if desired
      return quad_adapt<N>(f, el, a, mid, tol) + quad_adapt<N>(f, el, mid, b, tol);
   }

   // -----------------------------
   // 2D integration over one triangle
   // Notebook-style: build a rotation
   // from AB, then integrate along AB,
   // using two inner 1D integrals:
   //   Sfda : integrate in A-direction
   //   Sfdb : integrate in B-direction
   // Outer (edge) integral uses GL-3
   // Inner uses GL-8 (like the nb).
   // -----------------------------
   struct EdgeABContext
   {
      // Optional user fields (mimic elparams)
      long number_of_calls = 0;
      long sub_quad = 0;
      double A0 = 0.0; // starting A for inner A-integration
      double B0 = 0.0; // starting B for inner B-integration
   };

   // edge integral of the two inner contributions
   template <class F, class El>
   inline double GQ3S(F &&Eqn, El &el,
                      const std::array<double, 2> &n_ab, // normal in AB-frame (components along +A and +B)
                      const Mat2 &R,
                      const std::array<double, 2> &a0, // point A0B0 in AB-frame
                      const std::array<double, 2> &a1)
   {                             // point A1B1 in AB-frame
      const auto &xs = GL<8>::x; // we’ll use 8-pt along edge for robustness
      const auto &ws = GL<8>::w;
      const double mA = 0.5 * (a1[0] - a0[0]);
      const double cA = 0.5 * (a1[0] + a0[0]);
      const double mB = 0.5 * (a1[1] - a0[1]);
      const double cB = 0.5 * (a1[1] + a0[1]);

      double accA = 0.0;
      double accB = 0.0;
      for (size_t i = 0; i < xs.size(); ++i)
      {
         const double A = mA * xs[i] + cA;
         const double B = mB * xs[i] + cB;
         accA += ws[i] * Sfda(Eqn, el, R, A, B);
         accB += ws[i] * Sfdb(Eqn, el, R, A, B);
         //   accA += ws[i] * Sfda(Eqn, el, R, a0[0], A, B);
         //   accB += ws[i] * Sfdb(Eqn, el, R, a0[1], A, B);
      }
      const double scale = 0.5; // from mapping [-1,1] to edge
      return scale * (n_ab[0] * accA + n_ab[1] * accB);
   }

   // =====================
   // Utilities (3D / 4D)
   // =====================
   struct Vec2
   {
      double x, y;
   };
   inline Vec2 operator+(const Vec2 &a, const Vec2 &b) { return {a.x + b.x, a.y + b.y}; }
   inline Vec2 operator-(const Vec2 &a, const Vec2 &b) { return {a.x - b.x, a.y - b.y}; }
   inline Vec2 operator*(double s, const Vec2 &v) { return {s * v.x, s * v.y}; }

   inline double det2(const Vec2 &a, const Vec2 &b)
   {
      return a.x * b.y - a.y * b.x; // parallelogram area (signed)
   }

   struct Vec3
   {
      double x, y, z;
   };
   struct Vec4
   {
      double x, y, z, w;
   };

   inline Vec3 operator+(const Vec3 &a, const Vec3 &b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
   inline Vec3 operator-(const Vec3 &a, const Vec3 &b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
   inline Vec3 operator*(double s, const Vec3 &v) { return {s * v.x, s * v.y, s * v.z}; }

   inline Vec4 operator+(const Vec4 &a, const Vec4 &b) { return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w}; }
   inline Vec4 operator-(const Vec4 &a, const Vec4 &b) { return {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w}; }
   inline Vec4 operator*(double s, const Vec4 &v) { return {s * v.x, s * v.y, s * v.z, s * v.w}; }

   // det3 = det([a b c]) where a,b,c are column vectors
   inline double det3(const Vec3 &a, const Vec3 &b, const Vec3 &c)
   {
      return a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x) + a.z * (b.x * c.y - b.y * c.x);
   }

   // det4 = det of 4x4 whose columns are a,b,c,d (vectors from P0 to others)
   inline double det4(const Vec4 &a, const Vec4 &b, const Vec4 &c, const Vec4 &d)
   {
      // Laplace expansion (explicit 4x4 det)
      const double m11 = a.x, m12 = b.x, m13 = c.x, m14 = d.x;
      const double m21 = a.y, m22 = b.y, m23 = c.y, m24 = d.y;
      const double m31 = a.z, m32 = b.z, m33 = c.z, m34 = d.z;
      const double m41 = a.w, m42 = b.w, m43 = c.w, m44 = d.w;

      const double A = m11 * (m22 * (m33 * m44 - m34 * m43) - m23 * (m32 * m44 - m34 * m42) + m24 * (m32 * m43 - m33 * m42));
      const double B = -m12 * (m21 * (m33 * m44 - m34 * m43) - m23 * (m31 * m44 - m34 * m41) + m24 * (m31 * m43 - m33 * m41));
      const double C = m13 * (m21 * (m32 * m44 - m34 * m42) - m22 * (m31 * m44 - m34 * m41) + m24 * (m31 * m42 - m32 * m41));
      const double D = -m14 * (m21 * (m32 * m43 - m33 * m42) - m22 * (m31 * m43 - m33 * m41) + m23 * (m31 * m42 - m32 * m41));
      return A + B + C + D;
   }

   // ------------------------------
   // 1D: Gauss–Legendre on [a,b]
   // ------------------------------
   template <int N, class F, class El>
   double gl_integrate_1d(F &&f, El &el, double a, double b)
   {
      const auto &xs = quad::GL<N>::x;
      const auto &ws = quad::GL<N>::w;
      const double m = 0.5 * (b - a);
      const double c = 0.5 * (b + a);
      double acc = 0.0;
      for (size_t i = 0; i < xs.size(); ++i)
      {
         const double x = m * xs[i] + c;
         acc += ws[i] * f(el, x);
      }
      return m * acc;
   }

   // ================================
   // Gauss–Legendre on [0,1] (order N)
   // ================================
   template <int N>
   struct GL01
   {
      static constexpr auto &x01() { return GL<N>::x; } // reuse nodes in [-1,1]
      static constexpr auto &w01() { return GL<N>::w; }
      static double node(double xi) { return 0.5 * (xi + 1.0); } // map to [0,1]
      static double weight(double wi) { return 0.5 * wi; }       // scale weights
   };

   template <int N, class F, class El>
   double gl01_square(F &&f, El &el)
   {
      const auto &xs = GL<N>::x;
      const auto &ws = GL<N>::w;
      double acc = 0.0;
      for (size_t i = 0; i < xs.size(); ++i)
      {
         const double u = 0.5 * (xs[i] + 1.0);
         const double wu = 0.5 * ws[i];
         for (size_t j = 0; j < xs.size(); ++j)
         {
            const double v = 0.5 * (xs[j] + 1.0);
            const double wv = 0.5 * ws[j];
            acc += wu * wv * f(el, u, v);
         }
      }
      return acc;
   }

   // ================================
   // Tensor product GL on [0,1]^k
   // (3D & 4D specializations below)
   // ================================
   template <int N, class F, class El>
   double gl01_cube3(F &&f, El &el)
   {
      const auto &xs = GL<N>::x;
      const auto &ws = GL<N>::w;
      double acc = 0.0;
      for (size_t i = 0; i < xs.size(); ++i)
      {
         const double u = GL01<N>::node(xs[i]);
         const double wu = GL01<N>::weight(ws[i]);
         for (size_t j = 0; j < xs.size(); ++j)
         {
            const double v = GL01<N>::node(xs[j]);
            const double wv = GL01<N>::weight(ws[j]);
            for (size_t k = 0; k < xs.size(); ++k)
            {
               const double w = GL01<N>::node(xs[k]);
               const double ww = GL01<N>::weight(ws[k]);
               acc += wu * wv * ww * f(el, u, v, w);
            }
         }
      }
      return acc;
   }

   template <int N, class F, class El>
   double gl01_hyper4(F &&f, El &el)
   {
      const auto &xs = GL<N>::x;
      const auto &ws = GL<N>::w;
      double acc = 0.0;
      for (size_t i = 0; i < xs.size(); ++i)
      {
         const double u = GL01<N>::node(xs[i]);
         const double wu = GL01<N>::weight(ws[i]);
         for (size_t j = 0; j < xs.size(); ++j)
         {
            const double v = GL01<N>::node(xs[j]);
            const double wv = GL01<N>::weight(ws[j]);
            for (size_t k = 0; k < xs.size(); ++k)
            {
               const double w = GL01<N>::node(xs[k]);
               const double ww = GL01<N>::weight(ws[k]);
               for (size_t m = 0; m < xs.size(); ++m)
               {
                  const double z = GL01<N>::node(xs[m]);
                  const double wz = GL01<N>::weight(ws[m]);
                  acc += wu * wv * ww * wz * f(el, u, v, w, z);
               }
            }
         }
      }
      return acc;
   }

   // ========================================
   // Duffy transform for simplex integration
   // ========================================

   // Triangle integral by Duffy transform
   // Triangle vertices P0,P1,P2 (xy space)
   // Eqn signature: double Eqn(El&, double x, double y)
   template <int N, class F, class El>
   double integrate_triangle_duffy(F &&Eqn, El &el,
                                   const Vec2 &P0,
                                   const Vec2 &P1,
                                   const Vec2 &P2)
   {
      // Affine map from reference (r,s): X = P0 + r*(P1-P0) + s*(P2-P0),
      // with 0<=r, 0<=s, r+s<=1
      const Vec2 A = {P1.x - P0.x, P1.y - P0.y};
      const Vec2 B = {P2.x - P0.x, P2.y - P0.y};
      const double Jaff = std::fabs(det2(A, B)); // = 2 * area(triangle)

      // Duffy map from (u,v) in [0,1]^2:
      //   r = u
      //   s = (1 - u) * v
      // Jacobian of (r,s) wrt (u,v) is (1 - u)
      auto f_duffy = [&](El &e, double u, double v) -> double
      {
         const double r = u;
         const double s = (1.0 - u) * v;
         const Vec2 X = {P0.x + r * A.x + s * B.x,
                         P0.y + r * A.y + s * B.y};
         const double Jduf = (1.0 - u);
         return Eqn(e, X.x, X.y) * Jduf;
      };

      // ∫_Δ f = |det2(A,B)| * ∫_[0,1]^2 f(X(u,v)) (1-u) du dv
      return Jaff * gl01_square<N>(f_duffy, el);
   }

   // 3D Tetrahedron integral by Duffy
   // Tet(P0,P1,P2,P3); Eqn: (El&, x,y,z) -> double
   template <int N, class F, class El>
   double integrate_tetrahedron(F &&Eqn, El &el,
                                const Vec3 &P0, const Vec3 &P1,
                                const Vec3 &P2, const Vec3 &P3)
   {
      // Affine map: X = P0 + a*(P1-P0) + b*(P2-P0) + c*(P3-P0),  a,b,c >=0, a+b+c<=1
      const Vec3 A = P1 - P0, B = P2 - P0, C = P3 - P0;
      const double detM = std::fabs(det3(A, B, C)); // 6 * volume of the tet

      // Duffy map: (u,v,w) in [0,1]^3
      //   a = u
      //   b = (1-u) v
      //   c = (1-u)(1-v) w
      // Jacobian = (1-u)^2 (1-v)
      auto f_duffy = [&](El &e, double u, double v, double w) -> double
      {
         const double a = u;
         const double b = (1.0 - u) * v;
         const double c = (1.0 - u) * (1.0 - v) * w;
         const Vec3 X = P0 + a * A + b * B + c * C;
         const double J = (1.0 - u) * (1.0 - u) * (1.0 - v);
         return Eqn(e, X.x, X.y, X.z) * J;
      };

      // ∫_Tet f dV = ∫_[0,1]^3 f(X(u,v,w)) * J(u,v,w) du dv dw * |detM|
      return detM * gl01_cube3<N>(f_duffy, el);
   }

   // 4D 4-simplex integral by Duffy
   // 4-simplex(P0..P4); Eqn: (El&, x,y,z,w) -> double
   template <int N, class F, class El>
   double integrate_4simplex(F &&Eqn, El &el,
                             const Vec4 &P0, const Vec4 &P1,
                             const Vec4 &P2, const Vec4 &P3,
                             const Vec4 &P4)
   {
      // Affine map in 4D: X = P0 + a*(P1-P0) + b*(P2-P0) + c*(P3-P0) + d*(P4-P0)
      const Vec4 A = P1 - P0, B = P2 - P0, C = P3 - P0, D = P4 - P0;
      const double detM = std::fabs(det4(A, B, C, D)); // 24 * 4D hypervolume

      // Duffy map in 4D (u,v,w,z in [0,1]):
      //   a = u
      //   b = (1-u) v
      //   c = (1-u)(1-v) w
      //   d = (1-u)(1-v)(1-w) z
      // Jacobian = (1-u)^3 (1-v)^2 (1-w)
      auto f_duffy = [&](El &e, double u, double v, double w, double z) -> double
      {
         const double a = u;
         const double b = (1.0 - u) * v;
         const double c = (1.0 - u) * (1.0 - v) * w;
         const double d = (1.0 - u) * (1.0 - v) * (1.0 - w) * z;
         const Vec4 X = P0 + a * A + b * B + c * C + d * D;
         const double J = (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - v) * (1.0 - v) * (1.0 - w);
         return Eqn(e, X.x, X.y, X.z, X.w) * J;
      };

      return detM * gl01_hyper4<N>(f_duffy, el);
   }

   namespace quad_detail
   {

      // Rotation that aligns v0->v1 with +A in the CURRENT frame.
      // Returns (R_align, edge_len)
      inline std::pair<Mat2, double> align_in_frame(const std::array<double, 2> &v0,
                                                    const std::array<double, 2> &v1)
      {
         return normVec2D(v0, v1); // same math: rotation + length
      }

      // The notebook’s inner 1D integrals along A and B (GL-8)
      template <class F, class El>
      inline double Sfda(F &&Eqn, El &el, const Mat2 &R, double A, double B)
      {
         const auto &xs = GL<8>::x;
         const auto &ws = GL<8>::w;
         const double A0 = get_A0(el); // <- baseline (0.0 if none)
         double acc = 0.0;
         for (size_t i = 0; i < xs.size(); ++i)
         {
            const double a = 0.5 * ((1.0 + xs[i]) * A + (1.0 - xs[i]) * A0);
            const std::array<double, 2> v{a, B};
            const auto xy = R.mul(v);
            acc += ws[i] * Eqn(el, xy[0], xy[1]);
         }
         return 0.5 * (A - A0) * acc;
      }

      template <class F, class El>
      inline double Sfdb(F &&Eqn, El &el, const Mat2 &R, double A, double B)
      {
         const auto &xs = GL<8>::x;
         const auto &ws = GL<8>::w;
         const double B0 = get_B0(el); // <- baseline (0.0 if none)
         double acc = 0.0;
         for (size_t i = 0; i < xs.size(); ++i)
         {
            const double b = 0.5 * ((1.0 + xs[i]) * B + (1.0 - xs[i]) * B0);
            const std::array<double, 2> v{A, b};
            const auto xy = R.mul(v);
            acc += ws[i] * Eqn(el, xy[0], xy[1]);
         }
         return 0.5 * (B - B0) * acc;
      }

      // Outer edge integral (along AB) combining A/B inner parts.
      // NOTE: NO edge-length here (the nb multiplies that outside).
      template <class F, class El>
      inline double GQ_edge(F &&Eqn, El &el,
                            const std::array<double, 2> &n_ab, // rotated normal in AB frame
                            const Mat2 &R,
                            const std::array<double, 2> &a0, // AB-frame start
                            const std::array<double, 2> &a1)
      { // AB-frame end
         // We use GL-8 along the edge for robustness
         const auto &xs = GL<8>::x;
         const auto &ws = GL<8>::w;
         const double mA = 0.5 * (a1[0] - a0[0]);
         const double cA = 0.5 * (a1[0] + a0[0]);
         const double mB = 0.5 * (a1[1] - a0[1]);
         const double cB = 0.5 * (a1[1] + a0[1]);

         double accA = 0.0, accB = 0.0;
         for (size_t i = 0; i < xs.size(); ++i)
         {
            const double A = mA * xs[i] + cA;
            const double B = mB * xs[i] + cB;
            accA += ws[i] * Sfda(Eqn, el, R, A, B);
            accB += ws[i] * Sfdb(Eqn, el, R, A, B);

            // accA += ws[i] * Sfda(Eqn, el, R, a0[0], A, B);
            // accB += ws[i] * Sfdb(Eqn, el, R, a0[1], A, B);
         }
         // map factor [-1,1] -> edge parameter
         return 0.5 * (n_ab[0] * accA + n_ab[1] * accB);
      }

      // One sub-triangle contribution following the notebook:
      //  - Build R from the FIRST TWO xy points (P0->P1)
      //  - Transform those two into AB frame
      //  - Inside AB, build a SECOND rotation RAS that aligns that (still) edge with +A
      //  - normal_ab = RAS * [0, -1]
      //  - Multiply the edge integral by the edge length in AB (l_NormS)
      template <class F, class El>
      inline double sub_integrate(F &&Eqn, El &el,
                                  const Mat2 &R_base, // <— fixed outer R
                                  const std::array<double, 2> &P0,
                                  const std::array<double, 2> &P1,
                                  const std::array<double, 2> & /*P2*/)
      {
         // Endpoints transformed to AB frame using the fixed R
         const auto A0 = R_base.Tmul(P0);
         const auto A1 = R_base.Tmul(P1);

         // Secondary alignment INSIDE the AB frame, plus edge length
         auto RasPair = normVec2D(A0, A1);
         const Mat2 RAS = RasPair.first;
         const double L_edge = RasPair.second;

         // Notebook’s normal components in AB frame:  n_ab = RAS * [0, -1]
         const std::array<double, 2> n_ab = RAS.mul(std::array<double, 2>{0.0, -1.0});

         // Edge integral in AB frame; Eqn maps AB -> XY via the fixed outer R
         const double edge_val = GQ_edge(Eqn, el, n_ab, R_base, {A0[0], A0[1]}, {A1[0], A1[1]});

         return L_edge * edge_val; // multiply by line measure (like l_NormS in nb)
      }

   } // namespace quad_detail

   // Public: triangle integral = sum of the three notebook-style sub-triangles
   template <class F, class El>
   double integrate_triangle(F &&Eqn, El &el,
                             const std::array<double, 2> &p0,
                             const std::array<double, 2> &p1,
                             const std::array<double, 2> &p2)
   {
      Vec2 P0{p0[0], p0[1]}, P1{p1[0], p1[1]}, P2{p2[0], p2[1]};
      return integrate_triangle_duffy<8>(std::forward<F>(Eqn), el, P0, P1, P2);
   }
   // {
   //    using namespace quad_detail;

   //    // Outer rotation from the FIRST edge only (matches the notebook)
   //    const Mat2 R_base = normVec2D(p0, p1).first;

   //    const double c01 = sub_integrate(Eqn, el, R_base, p0, p1, p2);
   //    const double c12 = sub_integrate(Eqn, el, R_base, p1, p2, p0);
   //    const double c20 = sub_integrate(Eqn, el, R_base, p2, p0, p1);

   //    return 0.5 * (c01 + c12 + c20);
   // }

   // // Integral over a triangle (p0,p1,p2) by the notebook scheme
   // // We build R from (p0->p1). Then:
   // //   - transform first two vertices into AB frame
   // //   - outer integral uses edge along AB (with inner Sfda/Sfdb)
   // //   - sum three cyclic sub-triangles as in the nb
   // template<class F, class El>
   // double integrate_triangle(F&& Eqn, El& el,
   //                           const std::array<double,2>& p0,
   //                           const std::array<double,2>& p1,
   //                           const std::array<double,2>& p2)
   // {
   //     // Rotation aligning AB with +A axis
   //     auto Rpair = normVec2D(p0, p1);
   //     const Mat2 R = Rpair.first;

   //     // AB-frame points
   //     const auto pa = R.Tmul(p0);
   //     const auto pb = R.Tmul(p1);
   //     const auto pc = R.Tmul(p2);

   //     // normal (in AB frame) pointing "down" is (0, -1):
   //     // but we need components of that normal in AB basis for weighting:
   //     // n_ab = (na, nb) = R_AB * (0, -1) = columns of R multiplied by (0,-1).
   //     // In AB frame R_AB is identity; however, notebook mixes an extra rotation from
   //     // subdividing edges. Here, use the AB normal directly for the base edge:
   //     const std::array<double,2> n_ab{ 0.0, -1.0 };

   //     // Contribution from edge (p0,p1) with third point p2:
   //     const double c01 =
   //         GQ3S(Eqn, el, n_ab, R, /*a0*/ {pa[0], pa[1]}, /*a1*/ {pb[0], pb[1]});

   //     // Rotate the triangle cyclically to mimic the notebook’s 3 sums:
   //     // Build new frames for (p1,p2) and (p2,p0)
   //     auto Rpair12 = normVec2D(p1, p2);
   //     const Mat2 R12 = Rpair12.first;
   //     const auto p1_12 = R12.Tmul(p1);
   //     const auto p2_12 = R12.Tmul(p2);
   //     const double c12 =
   //         GQ3S(Eqn, el, {0.0, -1.0}, R12, {p1_12[0], p1_12[1]}, {p2_12[0], p2_12[1]});

   //     auto Rpair20 = normVec2D(p2, p0);
   //     const Mat2 R20 = Rpair20.first;
   //     const auto p2_20 = R20.Tmul(p2);
   //     const auto p0_20 = R20.Tmul(p0);
   //     const double c20 =
   //         GQ3S(Eqn, el, {0.0, -1.0}, R20, {p2_20[0], p2_20[1]}, {p0_20[0], p0_20[1]});

   //     return c01 + c12 + c20;
   // }

} // namespace quad
