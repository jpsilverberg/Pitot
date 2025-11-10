#include <gtest/gtest.h>
#include <dstl/dlog.h>

#include "../include/dstl/LinAlg.h"

using namespace dstl::linalg;

TEST(Mat, LU_SolveAndDeterminant)
{
   using M = Mat<3, 3, double>;

   M A{4, 2, 3,
       3, 1, 2,
       2, 1, 3};

   dstl::linalg::Vec<3, double> x1;
   dstl::linalg::Vec<3, double> x2;
   {
      auto F = LU<3, double>::factor(A);
      VAR(F);
      x1 = F.solve(Vec3{1, 2, 3});
      VAR(x1);
      auto detA = F.det();
      VAR(detA);
      auto invA = F.inverse();
      VAR(invA);
   }

   {
      auto F = LU_Left<3, double>::factor(A);
      VAR(F);
      x2 = F.solve(Vec3{1, 2, 3});
      VAR(x2);
      auto detA = F.det();
      VAR(detA);
      auto invA = F.inverse();
      VAR(invA);
   }

   auto dx = x1 - x2;
   VAR(dx);
}

TEST(Mat, StiffMatrix)
{

   using M = Mat<3, 3, double>;
   using V = Vec<3, double>;

   auto residual = [](const M &A, const V &x, const V &b)
   {
      V r = A * x - b;
      return norm(r);
   };

   auto reconstruct_err = [](const M &A, const M &LU, const Vec<3,uint32_t> &piv)
   {
      // build P*A
      M PA{};
      for (size_t i = 0; i < 3; ++i)
         for (size_t j = 0; j < 3; ++j)
            PA(i, j) = A(piv[i], j);
      // extract L,U
      M L = M::eye(), U{};
      for (size_t i = 0; i < 3; ++i)
      {
         for (size_t j = 0; j < 3; ++j)
         {
            if (i > j)
               L(i, j) = LU(i, j);
            else
               U(i, j) = LU(i, j);
         }
      }
      M R = L * U - PA;
      double e = 0.0;
      for (double x : R)
         e = std::max(e, std::abs(x));
      return e;
   };

   M A{1, 2, 3,
       4, 5, 6,
       7, 8, 9.001};

   dstl::linalg::Vec<3, double> x1;
   dstl::linalg::Vec<3, double> x2;
   dstl::linalg::Vec<3, double> rhs{1,2,3};

   auto F1 = LU<3, double>::factor(A);
   VAR(F1);
   x1 = F1.solve(rhs);
   VAR(x1);

   auto F2 = LU_Left<3, double>::factor(A);
   VAR(F2);
   x2 = F2.solve(rhs);
   VAR(x2);

   auto dx = x1 - x2;
   VAR(dx);

   auto r1 = residual(A, x1, rhs);
   EXPECT_NEAR(r1, 0.0, 1e-12);
   VAR(r1);
   auto r2 = residual(A, x2, rhs);
   VAR(r2);
   EXPECT_NEAR(r2, 0.0, 1e-12);

   EXPECT_NEAR(norm(x1 - x2), 0.0, 1e-12);

   // If you expose pivots from the factor (they are in your dumps), also:
   auto rr1 = reconstruct_err(A, F1.lu, F1.piv);
   VAR(rr1);
   EXPECT_LT(rr1, 1e-12);
   auto rr2 = reconstruct_err(A, F2.lu, F2.piv);
   VAR(rr2);
   EXPECT_LT(rr2, 1e-12);
}