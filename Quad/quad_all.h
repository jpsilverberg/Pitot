#pragma once
#include <array>
#include <cmath>
#include <functional>
#include <utility>
#include <type_traits>
#include <dstl/num/types/Float64.h>
#include <dstl/num/types/Accumulator.h>
#include <dstl/math/vec.h>

namespace quad
{
   using namespace dstl::math;
   struct AdaptParams
   {
      double abs_tol = 1e-10;
      double rel_tol = 1e-8;
      int max_depth = 20;
   };
   // // // //=========================//
   // // // // 0) Small math utilities //
   // // // //=========================//
   // // // struct Vec2
   // // // {
   // // //    double x, y;
   // // // };
   // // // struct Vec3
   // // // {
   // // //    double x, y, z;
   // // // };
   // // // struct Vec4
   // // // {
   // // //    double x, y, z, w;
   // // // };

   // // // inline Vec2 operator+(Vec2 a, Vec2 b) { return {a.x + b.x, a.y + b.y}; }
   // // // inline Vec2 operator-(Vec2 a, Vec2 b) { return {a.x - b.x, a.y - b.y}; }
   // // // inline Vec2 operator*(double s, Vec2 v) { return {s * v.x, s * v.y}; }
   // // // inline Vec2 operator*(Vec2 v, double s) { return {s * v.x, s * v.y}; }
   // // // inline double det2(Vec2 a, Vec2 b) { return a.x * b.y - a.y * b.x; }
   // // // inline double norm2(Vec2 v) { return std::sqrt(v.x * v.x + v.y * v.y); }

   // // // inline Vec3 operator+(Vec3 a, Vec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
   // // // inline Vec3 operator-(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
   // // // inline Vec3 operator*(double s, Vec3 v) { return {s * v.x, s * v.y, s * v.z}; }
   // // // inline Vec3 operator*(Vec3 v, double s) { return {s * v.x, s * v.y, s * v.z}; }
   // // // inline double dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
   // // // inline Vec3 cross(Vec3 a, Vec3 b)
   // // // {
   // // //    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
   // // // }
   // // // inline double norm(Vec3 v) { return std::sqrt(dot(v, v)); }

   // // // inline Vec4 operator+(Vec4 a, Vec4 b) { return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w}; }
   // // // inline Vec4 operator-(Vec4 a, Vec4 b) { return {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w}; }
   // // // inline Vec4 operator*(double s, Vec4 v) { return {s * v.x, s * v.y, s * v.z, s * v.w}; }
   // // // inline Vec4 operator*(Vec4 v, double s) { return {s * v.x, s * v.y, s * v.z, s * v.w}; }
   // // // inline double dot4(Vec4 a, Vec4 b) { return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w; }

   // // // // 4D area-normal of a 3D tetra face spanned by (A,B,C) (columns)
   // // // inline Vec4 area_normal4(Vec4 A, Vec4 B, Vec4 C)
   // // // {
   // // //    const double n0 = A.y * (B.z * C.w - B.w * C.z) - A.z * (B.y * C.w - B.w * C.y) + A.w * (B.y * C.z - B.z * C.y);
   // // //    const double n1 = -(A.x * (B.z * C.w - B.w * C.z) - A.z * (B.x * C.w - B.w * C.x) + A.w * (B.x * C.z - B.z * C.x));
   // // //    const double n2 = A.x * (B.y * C.w - B.w * C.y) - A.y * (B.x * C.w - B.w * C.x) + A.w * (B.x * C.y - B.y * C.x);
   // // //    const double n3 = -(A.x * (B.y * C.z - B.z * C.y) - A.y * (B.x * C.z - B.z * C.x) + A.z * (B.x * C.y - B.y * C.x));
   // // //    return {n0, n1, n2, n3};
   // // // }

   //===============================//
   // 1) Gauss–Legendre tables GL<N>//
   //===============================//
   // Provide N=8 and N=10 specializations (extend as needed).
   template <int N>
   struct GL; // primary template (no def)
   template <>
   struct GL<2>
   {
      static constexpr std::array<double, 2> x{
          -0.57735026918962576451, 0.57735026918962576451};
      static constexpr std::array<double, 2> w{
          1.0, 1.0};
   };

   template <>
   struct GL<3>
   {
      static constexpr std::array<double, 3> x{
          -0.77459666924148337704, 0.0, 0.77459666924148337704};
      static constexpr std::array<double, 3> w{
          0.55555555555555555556, 0.88888888888888888889, 0.55555555555555555556};
   };

   template <>
   struct GL<4>
   {
      static constexpr std::array<double, 4> x{
          -0.86113631159405257522, -0.33998104358485626480,
          0.33998104358485626480, 0.86113631159405257522};
      static constexpr std::array<double, 4> w{
          0.34785484513745385737, 0.65214515486254614263,
          0.65214515486254614263, 0.34785484513745385737};
   };

   template <>
   struct GL<5>
   {
      static constexpr std::array<double, 5> x{
          -0.90617984593866399280, -0.53846931010568309104, 0.0,
          0.53846931010568309104, 0.90617984593866399280};
      static constexpr std::array<double, 5> w{
          0.23692688505618908751, 0.47862867049936646804, 0.56888888888888888889,
          0.47862867049936646804, 0.23692688505618908751};
   };

   template <>
   struct GL<6>
   {
      static constexpr std::array<double, 6> x{
          -0.93246951420315202781, -0.66120938646626451366, -0.23861918608319690863,
          0.23861918608319690863, 0.66120938646626451366, 0.93246951420315202781};
      static constexpr std::array<double, 6> w{
          0.17132449237917034504, 0.36076157304813860757, 0.46791393457269104739,
          0.46791393457269104739, 0.36076157304813860757, 0.17132449237917034504};
   };

   template <>
   struct GL<7>
   {
      static constexpr std::array<double, 7> x{
          -0.94910791234275852453, -0.74153118559939443986, -0.40584515137739716691, 0.0,
          0.40584515137739716691, 0.74153118559939443986, 0.94910791234275852453};
      static constexpr std::array<double, 7> w{
          0.12948496616886969327, 0.27970539148927666790, 0.38183005050511894495, 0.41795918367346938776,
          0.38183005050511894495, 0.27970539148927666790, 0.12948496616886969327};
   };

   template <>
   struct GL<8>
   {
      static constexpr std::array<double, 8> x = {
          -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
          0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363};
      static constexpr std::array<double, 8> w = {
          0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
          0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};
   };

   template <>
   struct GL<10>
   {
      static constexpr std::array<double, 10> x = {
          -0.9739065285171717, -0.8650633666889845, -0.6794095682990244, -0.4333953941292472, -0.1488743389816312,
          0.1488743389816312, 0.4333953941292472, 0.6794095682990244, 0.8650633666889845, 0.9739065285171717};
      static constexpr std::array<double, 10> w = {
          0.0666713443086881, 0.1494513491505806, 0.2190863625159820, 0.2692667193099963, 0.2955242247147529,
          0.2955242247147529, 0.2692667193099963, 0.2190863625159820, 0.1494513491505806, 0.0666713443086881};
   };

   //==============================//
   // 2) Entity reference tags     //
   //==============================//
   enum class EntKind
   {
      Point = 0,
      Edge = 1,
      Face = 2,
      Cell3D = 3
   };
   struct EdgeRef2D
   {
      int id;
   }; // 0:(P0,P1),1:(P1,P2),2:(P2,P0)
   struct FaceRef3D
   {
      int id;
   }; // 0..3 faces opposite P0..P3
   struct FaceRef4D
   {
      int id;
   }; // 0..4 faces opposite P0..P4

   //==============================//
   // 3) Core 1D integration       //
   //==============================//

   // Volume on [a,b]: GL-N
   template <int N, class F, class El, class Accumulator = Accum>
   inline double integrate_segment(F &&f, El &el, double a, double b)
   {
      const auto &xs = GL<N>::x;
      const auto &ws = GL<N>::w;
      const double m = 0.5 * (b - a), c = 0.5 * (b + a);
      Accumulator acc{};
      for (size_t i = 0; i < xs.size(); ++i)
         acc += Accum(ws[i] * f(el, m * xs[i] + c));
      return m * acc.Double();
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
      const auto whole = integrate_segment<N>(f, el, a, b);
      const double mid = 0.5 * (a + b);
      const auto left = integrate_segment<N>(f, el, a, mid);
      const auto right = integrate_segment<N>(f, el, mid, b);
      const auto better = left + right;
      const double scale = std::fmax(std::fabs(better), 1.0);
      if (std::fabs(better - whole) <= std::fmax(tol.rtol * scale, tol.atol))
      {
         return better;
      }
      // Note: caller's El can track recursion if desired
      return quad_adapt<N>(f, el, a, mid, tol) + quad_adapt<N>(f, el, mid, b, tol);
   }

   // Boundary (points): separate laws or dispatcher
   template <class G0, class G1, class El>
   inline double boundary_points_segment(G0 &&gL, G1 &&gR, El &el, double a, double b)
   {
      return gL(el, a) + gR(el, b);
   }
   template <class G, class El>
   inline double boundary_points_segment_dispatch(G &&g, El &el, double a, double b)
   {
      return g(el, 0, a) + g(el, 1, b);
   }

   // Boundary (edge as 1D “surface”): GL-N on [a,b]
   template <int N, class G, class El>
   inline double boundary_edge_segment(G &&g, El &el, double a, double b)
   {
      return integrate_segment<N>(std::forward<G>(g), el, a, b);
   }

   //=================================//
   // 4) 2D triangle (volume + bndry) //
   //=================================//

   // Area (volume) via Duffy: r=u, s=(1-u)v; X = P0 + r*A + s*B
   template <int N, class F, class El, class Accumulator = Accum>
   inline double integrate_triangle_duffy(F &&Eqn, El &el, const Vec2 &P0, const Vec2 &P1, const Vec2 &P2)
   {
      const Vec2 A{P1.x - P0.x, P1.y - P0.y};
      const Vec2 B{P2.x - P0.x, P2.y - P0.y};
      const double Jaff = std::fabs(det2(A, B)); // = 2*area
      const auto &xs = GL<N>::x;
      const auto &ws = GL<N>::w;
      Accumulator acc{};
      for (size_t i = 0; i < xs.size(); ++i)
      {
         const double u = 0.5 * (xs[i] + 1.0);
         const double wu = 0.5 * ws[i];
         for (size_t j = 0; j < xs.size(); ++j)
         {
            const double v = 0.5 * (xs[j] + 1.0);
            const double wv = 0.5 * ws[j];
            const double r = u, s = (1.0 - u) * v, Jduf = (1.0 - u);
            const Vec2 X{P0.x + r * A.x + s * B.x, P0.y + r * A.y + s * B.y};
            acc += Accum(wu * wv * Eqn(el, X.x, X.y) * Jduf);
         }
      }
      return Jaff * acc.Double();
   }

   //================================================//
   // 2D: triangle adaptive via Duffy (unit square)  //
   //================================================//
   template <int N, class F, class El, class Acc = Accum>
   inline double integrate_triangle_adapt(F &&Eqn, El &el,
                                          const Vec2 &P0, const Vec2 &P1, const Vec2 &P2,
                                          AdaptParams p = {})
   {
      const auto &xs = GL<N>::x;
      const auto &ws = GL<N>::w;

      // Affine (2x2) Jacobian magnitude (== 2*area)
      const Vec2 A{P1.x - P0.x, P1.y - P0.y};
      const Vec2 B{P2.x - P0.x, P2.y - P0.y};
      const double Jaff = std::fabs(det2(A, B));

      auto box_eval = [&](double u0, double u1, double v0, double v1) -> double
      {
         const double mu = 0.5 * (u1 - u0), cu = 0.5 * (u1 + u0);
         const double mv = 0.5 * (v1 - v0), cv = 0.5 * (v1 + v0);
         Acc acc{};
         for (size_t i = 0; i < xs.size(); ++i)
         {
            const double u = mu * xs[i] + cu, wu = 0.5 * ws[i] * (u1 - u0);
            for (size_t j = 0; j < xs.size(); ++j)
            {
               const double v = mv * xs[j] + cv, wv = 0.5 * ws[j] * (v1 - v0);
               const double r = u, s = (1.0 - u) * v;
               const double Jduf = (1.0 - u);
               const Vec2 X{P0.x + r * A.x + s * B.x, P0.y + r * A.y + s * B.y};
               acc += Acc((wu * wv) * Eqn(el, X.x, X.y) * Jduf);
            }
         }
         return acc.Double();
      };

      std::function<double(double, double, double, double, int)> rec =
          [&](double u0, double u1, double v0, double v1, int depth) -> double
      {
         const double whole = box_eval(u0, u1, v0, v1);

         // split along larger side
         const double du = u1 - u0, dv = v1 - v0;
         if (du >= dv)
         {
            const double um = 0.5 * (u0 + u1);
            const double left = box_eval(u0, um, v0, v1);
            const double right = box_eval(um, u1, v0, v1);
            const double sum_children = left + right;
            const double err = std::abs(sum_children - whole);
            const double thr = p.abs_tol + p.rel_tol * std::max(std::abs(sum_children), std::abs(whole));
            if (err <= thr || depth >= p.max_depth)
               return sum_children;
            return rec(u0, um, v0, v1, depth + 1) + rec(um, u1, v0, v1, depth + 1);
         }
         else
         {
            const double vm = 0.5 * (v0 + v1);
            const double bot = box_eval(u0, u1, v0, vm);
            const double top = box_eval(u0, u1, vm, v1);
            const double sum_children = bot + top;
            const double err = std::abs(sum_children - whole);
            const double thr = p.abs_tol + p.rel_tol * std::max(std::abs(sum_children), std::abs(whole));
            if (err <= thr || depth >= p.max_depth)
               return sum_children;
            return rec(u0, u1, v0, vm, depth + 1) + rec(u0, u1, vm, v1, depth + 1);
         }
      };

      const double unit_int = rec(0.0, 1.0, 0.0, 1.0, 0);
      return Jaff * unit_int;
   }

   // Boundary: points
   template <class G0, class G1, class G2, class El>
   inline double boundary_points_triangle(G0 &&g0, G1 &&g1, G2 &&g2, El &el,
                                          const Vec2 &P0, const Vec2 &P1, const Vec2 &P2)
   {
      return g0(el, P0.x, P0.y) + g1(el, P1.x, P1.y) + g2(el, P2.x, P2.y);
   }
   template <class G, class El>
   inline double boundary_points_triangle_dispatch(G &&g, El &el,
                                                   const Vec2 &P0, const Vec2 &P1, const Vec2 &P2)
   {
      return g(el, 0, P0.x, P0.y) + g(el, 1, P1.x, P1.y) + g(el, 2, P2.x, P2.y);
   }

   // Boundary: per-edge scalar (array)
   template <int N, class G, class El, class Accumulator = Accum>
   inline double boundary_edges_triangle(
       const std::array<G, 3> &g, El &el, const Vec2 &P0, const Vec2 &P1, const Vec2 &P2)
   {
      auto edge = [&](int eid, Vec2 A, Vec2 B)
      {
         const Vec2 AB{B.x - A.x, B.y - A.y};
         const double L = norm2(AB);
         const auto &xs = GL<N>::x;
         const auto &ws = GL<N>::w;
         Accumulator acc{};
         for (size_t i = 0; i < xs.size(); ++i)
         {
            const double t = 0.5 * (xs[i] + 1.0);
            const Vec2 X{A.x + t * AB.x, A.y + t * AB.y};
            acc += Accum(ws[i] * g[eid](el, X.x, X.y));
         }

         return 0.5 * L * acc.Double();
      };
      return edge(0, P0, P1) + edge(1, P1, P2) + edge(2, P2, P0);
   }

   // Boundary: per-edge scalar (dispatcher)
   template <int N, class G, class El, class Accumulator = Accum>
   inline double boundary_edges_triangle_dispatch(G &&g, El &el,
                                                  const Vec2 &P0, const Vec2 &P1, const Vec2 &P2)
   {
      auto edge = [&](int eid, Vec2 A, Vec2 B)
      {
         const Vec2 AB{B.x - A.x, B.y - A.y};
         const double L = norm2(AB);
         const auto &xs = GL<N>::x;
         const auto &ws = GL<N>::w;
         Accumulator acc{};
         for (size_t i = 0; i < xs.size(); ++i)
         {
            const double t = 0.5 * (xs[i] + 1.0);
            const Vec2 X{A.x + t * AB.x, A.y + t * AB.y};
            acc += Accum(ws[i] * g(el, EdgeRef2D{eid}, X.x, X.y));
         }
         return 0.5 * L * acc.Double();
      };
      return edge(0, P0, P1) + edge(1, P1, P2) + edge(2, P2, P0);
   }

   // Boundary: per-edge flux (vector field · outward normal)
   template <int N, class Fv, class El, class Vec2Like, class Accumulator = Accum>
   inline double boundary_edges_triangle_flux_dispatch(Fv &&F, El &el,
                                                       const Vec2 &P0, const Vec2 &P1, const Vec2 &P2)
   {
      auto edge_flux = [&](int eid, Vec2 A, Vec2 B, Vec2 Opp)
      {
         const Vec2 t = {B.x - A.x, B.y - A.y};
         const double L = norm2(t);
         if (L < 1e-300)
            return 0.0;
         Vec2 n = {t.y / L, -t.x / L}; // right-hand outward candidate
         const Vec2 mid = {(A.x + B.x) * 0.5, (A.y + B.y) * 0.5};
         const Vec2 toOpp = {Opp.x - mid.x, Opp.y - mid.y};
         if (n.x * toOpp.x + n.y * toOpp.y > 0.0)
         {
            n.x = -n.x;
            n.y = -n.y;
         }
         const auto &xs = GL<N>::x;
         const auto &ws = GL<N>::w;
         Accumulator acc{};
         for (size_t i = 0; i < xs.size(); ++i)
         {
            const double tpar = 0.5 * (xs[i] + 1.0);
            const Vec2 X{A.x + tpar * t.x, A.y + tpar * t.y};
            const auto v = F(el, EdgeRef2D{eid}, X.x, X.y); // expects .x,.y
            const double fdotn = v.x * n.x + v.y * n.y;
            acc += Accum(ws[i] * fdotn);
         }
         return 0.5 * L * acc.Double();
      };
      return edge_flux(0, P0, P1, P2) + edge_flux(1, P1, P2, P0) + edge_flux(2, P2, P0, P1);
   }

   //==================================//
   // 5) 3D tetrahedron (vol + bndry) //
   //==================================//

   // Volume via 3D Duffy: (u,v,w) in [0,1]^3 -> (r,s,t)=(u,(1-u)v,(1-u)(1-v)w)
   template <int N, class F, class El>
   inline double integrate_tetrahedron(F &&Eqn, El &el,
                                       const Vec3 &P0, const Vec3 &P1, const Vec3 &P2, const Vec3 &P3)
   {
      const Vec3 A{P1.x - P0.x, P1.y - P0.y, P1.z - P0.z};
      const Vec3 B{P2.x - P0.x, P2.y - P0.y, P2.z - P0.z};
      const Vec3 C{P3.x - P0.x, P3.y - P0.y, P3.z - P0.z};
      const double Jaff = std::fabs(dot(cross(A, B), C)); // = 6 * volume
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
            for (size_t k = 0; k < xs.size(); ++k)
            {
               const double w = 0.5 * (xs[k] + 1.0);
               const double ww = 0.5 * ws[k];
               const double r = u, s = (1.0 - u) * v, t = (1.0 - u) * (1.0 - v) * w;
               const Vec3 X{P0.x + r * A.x + s * B.x + t * C.x,
                            P0.y + r * A.y + s * B.y + t * C.y,
                            P0.z + r * A.z + s * B.z + t * C.z};
               const double Jduf = (1.0 - u) * (1.0 - u) * (1.0 - v);
               acc += wu * wv * ww * Eqn(el, X.x, X.y, X.z) * Jduf;
            }
         }
      }
      return Jaff * acc;
   }

   // Boundary: points
   template <class G0, class G1, class G2, class G3, class El>
   inline double boundary_points_tetra(G0 &&g0, G1 &&g1, G2 &&g2, G3 &&g3, El &el,
                                       const Vec3 &P0, const Vec3 &P1, const Vec3 &P2, const Vec3 &P3)
   {
      return g0(el, P0.x, P0.y, P0.z) + g1(el, P1.x, P1.y, P1.z) + g2(el, P2.x, P2.y, P2.z) + g3(el, P3.x, P3.y, P3.z);
   }
   template <class G, class El>
   inline double boundary_points_tetra_dispatch(G &&g, El &el,
                                                const Vec3 &P0, const Vec3 &P1, const Vec3 &P2, const Vec3 &P3)
   {
      return g(el, 0, P0.x, P0.y, P0.z) + g(el, 1, P1.x, P1.y, P1.z) + g(el, 2, P2.x, P2.y, P2.z) + g(el, 3, P3.x, P3.y, P3.z);
   }

   // Boundary: edges (6 of them), array or dispatcher
   template <int N, class G, class El>
   inline double boundary_edges_tetra(const std::array<G, 6> &g, El &el,
                                      const Vec3 &P0, const Vec3 &P1, const Vec3 &P2, const Vec3 &P3)
   {
      auto edge = [&](int eid, Vec3 A, Vec3 B)
      {
         const Vec3 AB{B.x - A.x, B.y - A.y, B.z - A.z};
         const double L = norm(AB);
         const auto &xs = GL<N>::x;
         const auto &ws = GL<N>::w;
         double acc = 0.0;
         for (size_t i = 0; i < xs.size(); ++i)
         {
            const double t = 0.5 * (xs[i] + 1.0);
            const Vec3 X{A.x + t * AB.x, A.y + t * AB.y, A.z + t * AB.z};
            acc += ws[i] * g[eid](el, X.x, X.y, X.z);
         }
         return 0.5 * L * acc;
      };
      return edge(0, P0, P1) + edge(1, P1, P2) + edge(2, P2, P0) + edge(3, P0, P3) + edge(4, P1, P3) + edge(5, P2, P3);
   }
   template <int N, class G, class El>
   inline double boundary_edges_tetra_dispatch(G &&g, El &el,
                                               const Vec3 &P0, const Vec3 &P1, const Vec3 &P2, const Vec3 &P3)
   {
      auto edge = [&](int eid, Vec3 A, Vec3 B)
      {
         const Vec3 AB{B.x - A.x, B.y - A.y, B.z - A.z};
         const double L = norm(AB);
         const auto &xs = GL<N>::x;
         const auto &ws = GL<N>::w;
         double acc = 0.0;
         for (size_t i = 0; i < xs.size(); ++i)
         {
            const double t = 0.5 * (xs[i] + 1.0);
            const Vec3 X{A.x + t * AB.x, A.y + t * AB.y, A.z + t * AB.z};
            acc += ws[i] * g(el, /*edge_id*/ eid, X.x, X.y, X.z);
         }
         return 0.5 * L * acc;
      };
      return edge(0, P0, P1) + edge(1, P1, P2) + edge(2, P2, P0) + edge(3, P0, P3) + edge(4, P1, P3) + edge(5, P2, P3);
   }

   // Boundary: faces (4 triangles), scalar per-face (array)
   template <int N, class G, class El>
   inline double boundary_faces_tetra(const std::array<G, 4> &g, El &el,
                                      const Vec3 &P0, const Vec3 &P1, const Vec3 &P2, const Vec3 &P3)
   {
      auto face2D = [&](int fid, Vec3 A, Vec3 B, Vec3 C)
      {
         const Vec3 a{B.x - A.x, B.y - A.y, B.z - A.z};
         const Vec3 b{C.x - A.x, C.y - A.y, C.z - A.z};
         const double area2 = norm(cross(a, b)); // = 2*area
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
               const double r = u, s = (1.0 - u) * v, J = (1.0 - u);
               const Vec3 X{A.x + r * a.x + s * b.x, A.y + r * a.y + s * b.y, A.z + r * a.z + s * b.z};
               acc += wu * wv * g[fid](el, X.x, X.y, X.z) * J;
            }
         }
         return area2 * acc;
      };
      return face2D(0, P1, P2, P3) + face2D(1, P0, P3, P2) + face2D(2, P0, P1, P3) + face2D(3, P0, P2, P1);
   }

   // Boundary: faces (dispatcher) & flux (F·n dS)
   template <int N, class G, class El>
   inline double boundary_faces_tetra_dispatch(G &&g, El &el,
                                               const Vec3 &P0, const Vec3 &P1, const Vec3 &P2, const Vec3 &P3)
   {
      auto face2D = [&](int fid, Vec3 A, Vec3 B, Vec3 C)
      {
         const Vec3 a{B.x - A.x, B.y - A.y, B.z - A.z};
         const Vec3 b{C.x - A.x, C.y - A.y, C.z - A.z};
         const double area2 = norm(cross(a, b));
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
               const double r = u, s = (1.0 - u) * v, J = (1.0 - u);
               const Vec3 X{A.x + r * a.x + s * b.x, A.y + r * a.y + s * b.y, A.z + r * a.z + s * b.z};
               acc += wu * wv * g(el, FaceRef3D{fid}, X.x, X.y, X.z) * J;
            }
         }
         return area2 * acc;
      };
      return face2D(0, P1, P2, P3) + face2D(1, P0, P3, P2) + face2D(2, P0, P1, P3) + face2D(3, P0, P2, P1);
   }

   // Flux: F returns Vec3 (fx,fy,fz); we integrate F·N (N = area-normal)
   template <int N, class Fv, class El>
   inline double boundary_faces_tetra_flux_dispatch(Fv &&F, El &el,
                                                    const Vec3 &P0, const Vec3 &P1, const Vec3 &P2, const Vec3 &P3)
   {
      auto face_flux = [&](int fid, Vec3 A, Vec3 B, Vec3 C, Vec3 Opp)
      {
         Vec3 a{B.x - A.x, B.y - A.y, B.z - A.z};
         Vec3 b{C.x - A.x, C.y - A.y, C.z - A.z};
         Vec3 NN = cross(a, b); // area-normal
         if (dot(NN, Opp - A) > 0.0)
            NN = NN * (-1.0); // ensure outward
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
               const double r = u, s = (1.0 - u) * v, J = (1.0 - u);
               const Vec3 X{A.x + r * a.x + s * b.x, A.y + r * a.y + s * b.y, A.z + r * a.z + s * b.z};
               const auto V = F(el, FaceRef3D{fid}, X.x, X.y, X.z);
               const double fdotn = V.x * NN.x + V.y * NN.y + V.z * NN.z;
               acc += wu * wv * fdotn * J;
            }
         }
         return acc; // N already contains area scaling
      };
      return face_flux(0, P1, P2, P3, P0) + face_flux(1, P0, P3, P2, P1) + face_flux(2, P0, P1, P3, P2) + face_flux(3, P0, P2, P1, P3);
   }

   //====================================//
   // 6) 4D 4-simplex (vol + 3D boundary)//
   //====================================//

   // Volume via 4D Duffy: (u,v,w,t)->(r,s,q,p) = (u,(1-u)v,(1-u)(1-v)w,(1-u)(1-v)(1-w)t)
   template <int N, class F, class El>
   inline double integrate_4simplex(F &&Eqn, El &el,
                                    const Vec4 &P0, const Vec4 &P1, const Vec4 &P2, const Vec4 &P3, const Vec4 &P4)
   {
      const Vec4 A{P1.x - P0.x, P1.y - P0.y, P1.z - P0.z, P1.w - P0.w};
      const Vec4 B{P2.x - P0.x, P2.y - P0.y, P2.z - P0.z, P2.w - P0.w};
      const Vec4 C{P3.x - P0.x, P3.y - P0.y, P3.z - P0.z, P3.w - P0.w};
      const Vec4 D{P4.x - P0.x, P4.y - P0.y, P4.z - P0.z, P4.w - P0.w};

      // |det([A B C D])| = 24 * hypervolume. Compute via Gram determinant sqrt(det(M^T M))
      // Explicit 4x4 det helper (Laplace expansion) for speed & simplicity:
      auto det4 = [](const Vec4 &a, const Vec4 &b, const Vec4 &c, const Vec4 &d) -> double
      {
         // build matrix columns a,b,c,d; expand (compact but fine for header use)
         const double m00 = a.x, m01 = b.x, m02 = c.x, m03 = d.x;
         const double m10 = a.y, m11 = b.y, m12 = c.y, m13 = d.y;
         const double m20 = a.z, m21 = b.z, m22 = c.z, m23 = d.z;
         const double m30 = a.w, m31 = b.w, m32 = c.w, m33 = d.w;
         const double c00 = m11 * (m22 * m33 - m23 * m32) - m12 * (m21 * m33 - m23 * m31) + m13 * (m21 * m32 - m22 * m31);
         const double c01 = m10 * (m22 * m33 - m23 * m32) - m12 * (m20 * m33 - m23 * m30) + m13 * (m20 * m32 - m22 * m30);
         const double c02 = m10 * (m21 * m33 - m23 * m31) - m11 * (m20 * m33 - m23 * m30) + m13 * (m20 * m31 - m21 * m30);
         const double c03 = m10 * (m21 * m32 - m22 * m31) - m11 * (m20 * m32 - m22 * m30) + m12 * (m20 * m31 - m21 * m30);
         return m00 * c00 - m01 * c01 + m02 * c02 - m03 * c03;
      };
      const double Jaff = std::fabs(det4(A, B, C, D)); // = 24 * hypervolume

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
            for (size_t k = 0; k < xs.size(); ++k)
            {
               const double w = 0.5 * (xs[k] + 1.0);
               const double ww = 0.5 * ws[k];
               for (size_t m = 0; m < xs.size(); ++m)
               {
                  const double t = 0.5 * (xs[m] + 1.0);
                  const double wt = 0.5 * ws[m];
                  const double r = u;
                  const double s = (1.0 - u) * v;
                  const double q = (1.0 - u) * (1.0 - v) * w;
                  const double p = (1.0 - u) * (1.0 - v) * (1.0 - w) * t;
                  const Vec4 X{P0.x + r * A.x + s * B.x + q * C.x + p * D.x,
                               P0.y + r * A.y + s * B.y + q * C.y + p * D.y,
                               P0.z + r * A.z + s * B.z + q * C.z + p * D.z,
                               P0.w + r * A.w + s * B.w + q * C.w + p * D.w};
                  const double Jduf = (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - v) * (1.0 - v) * (1.0 - w);
                  acc += wu * wv * ww * wt * Eqn(el, X.x, X.y, X.z, X.w) * Jduf;
               }
            }
         }
      }
      return Jaff * acc;
   }

   // Boundary: points (5)
   template <class G0, class G1, class G2, class G3, class G4, class El>
   inline double boundary_points_4simplex(G0 &&g0, G1 &&g1, G2 &&g2, G3 &&g3, G4 &&g4, El &el,
                                          const Vec4 &P0, const Vec4 &P1, const Vec4 &P2, const Vec4 &P3, const Vec4 &P4)
   {
      return g0(el, P0.x, P0.y, P0.z, P0.w) + g1(el, P1.x, P1.y, P1.z, P1.w) + g2(el, P2.x, P2.y, P2.z, P2.w) + g3(el, P3.x, P3.y, P3.z, P3.w) + g4(el, P4.x, P4.y, P4.z, P4.w);
   }

   // Boundary: 3D faces (5 tets), scalar per-face (array)
   template <int N, class G, class El>
   inline double boundary_faces_4simplex(const std::array<G, 5> &g, El &el,
                                         const Vec4 &P0, const Vec4 &P1, const Vec4 &P2, const Vec4 &P3, const Vec4 &P4)
   {
      auto face3D = [&](int fid, Vec4 A, Vec4 B, Vec4 C, Vec4 D)
      {
         const Vec4 a{B.x - A.x, B.y - A.y, B.z - A.z, B.w - A.w};
         const Vec4 b{C.x - A.x, C.y - A.y, C.z - A.z, C.w - A.w};
         const Vec4 c{D.x - A.x, D.y - A.y, D.z - A.z, D.w - A.w};
         const Vec4 NN = area_normal4(a, b, c);
         const double area3 = std::sqrt(dot(NN, NN)); // 3-volume scale
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
               for (size_t k = 0; k < xs.size(); ++k)
               {
                  const double w = 0.5 * (xs[k] + 1.0);
                  const double ww = 0.5 * ws[k];
                  const double r = u, s = (1.0 - u) * v, t = (1.0 - u) * (1.0 - v) * w;
                  const Vec4 X{A.x + r * a.x + s * b.x + t * c.x,
                               A.y + r * a.y + s * b.y + t * c.y,
                               A.z + r * a.z + s * b.z + t * c.z,
                               A.w + r * a.w + s * b.w + t * c.w};
                  const double J = (1.0 - u) * (1.0 - u) * (1.0 - v); // 3D Duffy
                  acc += wu * wv * ww * g[fid](el, X.x, X.y, X.z, X.w) * J;
               }
            }
         }
         return area3 * acc;
      };
      return face3D(0, P1, P2, P3, P4) + face3D(1, P0, P3, P2, P4) + face3D(2, P0, P1, P3, P4) + face3D(3, P0, P2, P1, P4) + face3D(4, P0, P1, P2, P3);
   }

   // (Optional) 4D flux variant would dot a Vec4 field with the oriented 4D area-normal;
   // wire analogously to 3D flux using area_normal4 and orientation against the opposite vertex.

} // namespace quad
