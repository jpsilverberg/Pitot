#include <gtest/gtest.h>
#include <thread>

#include "Numbers.h"
#include <dstl/dlog.h>
#include <iostream>
#include <array>

using namespace dstl;
using namespace quad;

struct Ctx
{
   long NumberOfCalls = 0;
   long SubQuad = 0;
};

template <typename T>
inline T sqr(T x) { return x * x; }

// 1D test: ∫_0^1 x^15 dx = 1/16
double f1(Ctx &el, double x)
{
   (void)el;
   return std::pow(x, 15);
}

// 2D test: f(x,y)=x+y on triangle [(2,2),(5,2),(2,5)] → 27
double f2(Ctx &el, double x, double y)
{
   ++el.NumberOfCalls;
   return x + y;
}

TEST(Quad, Test1)
{

   using namespace quad;

   // 1D adaptive
   {
      Ctx el1;
      double i1 = quad_adapt<8>(f1, el1, 0.0, 1.0);
      std::cout << "I1 ~= " << i1 << " (expected 1/16 = " << 1.0 / 16.0 << ")\n";
   }

   // 2D triangle
   {
      Ctx el2;
      quad::Vec2 p0{2.0, 2.0}, p1{5.0, 2.0}, p2{2.0, 5.0};
      auto wrapper = [](Ctx &e, double x)
      { return f1(e, x); }; // just shows the 1D signature exists
      (void)wrapper;

      auto Fxy = [](Ctx &e, double x, double y) -> double
      {
         ++e.NumberOfCalls;
         return x + y;
      };

      double Itri = integrate_triangle_duffy<8>(Fxy, el2, p0, p1, p2);
      std::cout << "I_triangle ~= " << Itri << "  calls=" << el2.NumberOfCalls << "\n";
   }
   {
      Ctx el2;
      quad::Vec2 p0{2.0, 2.0}, p1{5.0, 2.0}, p2{2.0, 5.0};
      auto F1 = [](Ctx &, double, double)
      { return 1.0; };
      double A = quad::integrate_triangle_duffy<8>(F1, el2, p0, p1, p2); // should be 4.5 (area)
      std::cout << "Area triangle = " << A << "\n";
   }

   {
      // Define vertices
      Vec3 p0{2, 2, 2}, p1{5, 2, 2}, p2{2, 5, 2}, p3{2, 2, 5};

      // Context (your existing Ctx works)
      Ctx el{};

      // Integrand examples
      auto one3 = [](Ctx &, double, double, double)
      { return 1.0; };
      auto lin3 = [](Ctx &, double x, double y, double z)
      { return x + y + z; };

      // Choose order (N = 6..8 is very accurate)
      double vol = integrate_tetrahedron<6>(one3, el, p0, p1, p2, p3); // ~ volume
      VAR(vol);
      double Ilin = integrate_tetrahedron<6>(lin3, el, p0, p1, p2, p3);
      VAR(Ilin);
   }

   {
      Ctx el{};

      // Right triangle (2,2)-(5,2)-(2,5): area 4.5
      Vec2 t0{2, 2}, t1{5, 2}, t2{2, 5};

      auto one2 = [](Ctx &, double, double)
      { return 1.0; };
      auto lin2 = [](Ctx &, double x, double y)
      { return x + y; };

      double A = integrate_triangle_duffy<8>(one2, el, t0, t1, t2); // 4.5
      VAR(A);
      double Ixy = integrate_triangle_duffy<8>(lin2, el, t0, t1, t2); // 27.0
      VAR(Ixy);

      // If you kept the alias:
      double A2 = integrate_triangle_duffy<8>(one2, el, {2, 2}, {5, 2}, {2, 5}); // 4.5
      VAR(A2);
      double I2 = integrate_triangle_duffy<8>(lin2, el, {2, 2}, {5, 2}, {2, 5}); // 27.0
      VAR(I2);
   }
}

// ======================================================
// 1D TESTS
// ======================================================

TEST(Quad_1D, PolynomialOnUnitInterval)
{
   Ctx el{};
   auto p0 = [](Ctx &e, double)
   { ++e.NumberOfCalls; return 1.0; };
   auto p1 = [](Ctx &e, double x)
   { ++e.NumberOfCalls; return x; };
   auto p5 = [](Ctx &e, double x)
   { ++e.NumberOfCalls; return std::pow(x, 5); };

   constexpr double a = 0.0, b = 1.0;
   // Exact integrals
   const double I0 = 1.0;
   const double I1 = 0.5;
   const double I5 = 1.0 / 6.0;

   auto I0accum = integrate_segment<8>(p0, el, a, b);
   VAR(I0accum);
   EXPECT_NEAR(I0accum, I0, 1e-14);
   auto I1accum = integrate_segment<8>(p1, el, a, b);
   VAR(I1accum);
   EXPECT_NEAR(I1accum, I1, 1e-14);
   auto I5accum = integrate_segment<8>(p5, el, a, b);
   VAR(I5accum);
   EXPECT_NEAR(I5accum, I5, 1e-14);
}

TEST(Quad_1D, PolynomialOnShiftedInterval)
{
   Ctx el{};
   auto p3 = [](Ctx &e, double x)
   { ++e.NumberOfCalls; return x * x * x; };
   const double a = -2.0, b = 3.0;
   const double exact = (std::pow(b, 4) - std::pow(a, 4)) / 4.0; // ∫ x^3 dx = x^4/4
   EXPECT_NEAR(integrate_segment<8>(p3, el, a, b), exact, 1e-13);
}

TEST(Quad_1D, SinOn0Pi)
{
   Ctx el{};
   auto s = [](Ctx &e, double x)
   { ++e.NumberOfCalls; return std::sin(x); };
   const double a = 0.0, b = M_PI;
   const double exact = 2.0; // ∫_0^π sin x dx = 2
   auto Isin = quad_adapt<8>(s, el, a, b);
   VAR(Isin);
   EXPECT_NEAR(Isin, exact, 1e-13);
}

// ======================================================
// 2D TRIANGLE (Duffy)
// ======================================================

TEST(Quad_2D, ReferenceTriangle_ConstantsAndLinear)
{
   Ctx el{};
   // Reference triangle: (0,0), (1,0), (0,1) -> area = 1/2
   Vec2 p0{0, 0}, p1{1, 0}, p2{0, 1};
   auto one = [](Ctx &e, double , double )
   { ++e.NumberOfCalls; return 1.0; };
   auto xfn = [](Ctx &e, double x, double )
   { ++e.NumberOfCalls; return x; };
   auto yfn = [](Ctx &e, double , double y)
   { ++e.NumberOfCalls; return y; };
   auto sfn = [](Ctx &e, double x, double y)
   { ++e.NumberOfCalls; return x + y; };

   const double A = 0.5;
   EXPECT_NEAR(quad::integrate_triangle_adapt<8>(one, el, p0, p1, p2), A, 1e-14);
   EXPECT_NEAR(quad::integrate_triangle_adapt<8>(xfn, el, p0, p1, p2), 1.0 / 6.0, 1e-14);
   EXPECT_NEAR(quad::integrate_triangle_adapt<8>(yfn, el, p0, p1, p2), 1.0 / 6.0, 1e-14);
   EXPECT_NEAR(quad::integrate_triangle_adapt<8>(sfn, el, p0, p1, p2), 1.0 / 3.0, 1e-14);
}

TEST(Quad_2D, SkewedTriangle_ConstantsAndLinear)
{
   Ctx el{};
   // Triangle from your earlier test: (2,2), (5,2), (2,5)
   Vec2 p0{2, 2}, p1{5, 2}, p2{2, 5};
   auto one = [](Ctx &e, double , double )
   { ++e.NumberOfCalls; return 1.0; };
   auto sfn = [](Ctx &e, double x, double y)
   { ++e.NumberOfCalls; return x + y; };

   // Area = |det([p1-p0, p2-p0])|/2 = |det([3,0; 0,3])|/2 = 9/2 = 4.5
   EXPECT_NEAR(quad::integrate_triangle_duffy<8>(one, el, p0, p1, p2), 4.5, 1e-13);
   // ∫(x+y) dA = (x+y at centroid) * area
   // centroid = ( (2+5+2)/3, (2+2+5)/3 ) = (3, 3)  => x+y = 6;  6 * 4.5 = 27
   EXPECT_NEAR(quad::integrate_triangle_duffy<8>(sfn, el, p0, p1, p2), 27.0, 1e-12);
}

// ======================================================
// 3D TETRAHEDRON (Duffy)
// ======================================================

TEST(Quad_3D, ReferenceTetra_ConstantsAndLinear)
{
   Ctx el{};
   Vec3 p0{0, 0, 0}, p1{1, 0, 0}, p2{0, 1, 0}, p3{0, 0, 1};
   auto one = [](Ctx &e, double , double , double )
   { ++e.NumberOfCalls; return 1.0; };
   auto sfn = [](Ctx &e, double x, double y, double z)
   { ++e.NumberOfCalls; return x+y+z; };

   // Volume = 1/6
   EXPECT_NEAR(quad::integrate_tetrahedron<8>(one, el, p0, p1, p2, p3), 1.0 / 6.0, 1e-14);
   // Centroid sum = 3*(1/4) = 3/4, integral = (3/4) * (1/6) = 1/8
   EXPECT_NEAR(quad::integrate_tetrahedron<8>(sfn, el, p0, p1, p2, p3), 1.0 / 8.0, 1e-14);
}

TEST(Quad_3D, AxisAlignedTetra_ConstantsAndLinear)
{
   Ctx el{};
   // Your earlier 3x3x3 along axes:
   Vec3 p0{2, 2, 2}, p1{5, 2, 2}, p2{2, 5, 2}, p3{2, 2, 5};
   auto one = [](Ctx &e, double , double , double )
   { ++e.NumberOfCalls; return 1.0; };
   auto sfn = [](Ctx &e, double x, double y, double z)
   { ++e.NumberOfCalls; return x+y+z; };

   // Volume = |det([3,0,0; 0,3,0; 0,0,3])|/6 = 27/6 = 4.5
   EXPECT_NEAR(quad::integrate_tetrahedron<8>(one, el, p0, p1, p2, p3), 4.5, 1e-12);
   // Centroid = (2.75, 2.75, 2.75) -> sum = 8.25; integral = 8.25 * 4.5 = 37.125
   EXPECT_NEAR(quad::integrate_tetrahedron<8>(sfn, el, p0, p1, p2, p3), 37.125, 1e-12);
}

// ======================================================
// 4D 4-SIMPLEX (Duffy)
// ======================================================

TEST(Quad_4D, Unit4Simplex_ConstantsAndLinear)
{
   Ctx el{};
   Vec4 q0{0, 0, 0, 0}, q1{1, 0, 0, 0}, q2{0, 1, 0, 0}, q3{0, 0, 1, 0}, q4{0, 0, 0, 1};
   auto one = [](Ctx &e, double , double , double , double )
   { ++e.NumberOfCalls; return 1.0; };
   auto sfn = [](Ctx &e, double x, double y, double z, double w)
   { ++e.NumberOfCalls; return x+y+z+w; };

   // Hypervolume = 1/24
   EXPECT_NEAR(quad::integrate_4simplex<8>(one, el, q0, q1, q2, q3, q4), 1.0 / 24.0, 1e-14);
   // Centroid = (1/5,...,1/5), sum = 4/5; integral = (4/5) * (1/24) = 1/30
   EXPECT_NEAR(quad::integrate_4simplex<8>(sfn, el, q0, q1, q2, q3, q4), 1.0 / 30.0, 1e-14);
}
