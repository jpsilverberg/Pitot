#include <gtest/gtest.h>

#include "../include/dstl/LinAlg.h"

using namespace dstl::linalg;

TEST(VecBasics, AddScale) {
  Vec3 a{1.0,2.0,3.0};
  Vec3 b{4.0,5.0,6.0};
  auto c = a + b;            // (5,7,9)
  auto d = 2.0 * a;          // (2,4,6)
  EXPECT_DOUBLE_EQ(c[0], 5.0); EXPECT_DOUBLE_EQ(c[1], 7.0); EXPECT_DOUBLE_EQ(c[2], 9.0);
  EXPECT_DOUBLE_EQ(d[0], 2.0); EXPECT_DOUBLE_EQ(d[1], 4.0); EXPECT_DOUBLE_EQ(d[2], 6.0);
}

TEST(VecDotNorm, DotAndNorm) {
  Vec3 x{3.0, 4.0, 12.0};
  EXPECT_DOUBLE_EQ(x.sqnorm(), 3*3 + 4*4 + 12*12);
  EXPECT_DOUBLE_EQ(x.dot(x), x.sqnorm());
  EXPECT_DOUBLE_EQ(x.norm(), std::sqrt(x.sqnorm()));
}

TEST(Vec2, DetAndCrossZ) {
  Vec2 a{1.0, 0.0}, b{0.0, 2.0};
  EXPECT_DOUBLE_EQ(det2(a,b), 2.0);
  EXPECT_DOUBLE_EQ(crossz(a,b), 2.0);
}

TEST(Vec3, CrossRightHanded) {
  Vec3 ex{1,0,0}, ey{0,1,0}, ez{0,0,1};
  auto r = cross(ex, ey);
  EXPECT_DOUBLE_EQ(r[0], ez[0]);
  EXPECT_DOUBLE_EQ(r[1], ez[1]);
  EXPECT_DOUBLE_EQ(r[2], ez[2]);
}

TEST(Vec4, AreaNormal4_BasisFace) {
  // Face spanned by columns A=(1,0,0,0), B=(0,1,0,0), C=(0,0,1,0) in R^4 -> normal should be along +w
  Vec4 A{1,0,0,0}, B{0,1,0,0}, C{0,0,1,0};
  auto n = area_normal4(A,B,C);
  // (0,0,0, - (1*1 - 0))? For this orientation, formula yields n=(0,0,0,-1)
  EXPECT_DOUBLE_EQ(n[0], 0.0);
  EXPECT_DOUBLE_EQ(n[1], 0.0);
  EXPECT_DOUBLE_EQ(n[2], 0.0);
  EXPECT_DOUBLE_EQ(n[3], -1.0);
}
