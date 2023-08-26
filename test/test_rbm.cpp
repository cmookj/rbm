//===------------------------------------------------------------*- C++ -*-===//
//
//  TestRigidBodyMotion.cpp
//  TestRigidBodyMotionLibrary
//
//  Created by Changmook Chun on 11/23/22.
//
//===----------------------------------------------------------------------===//

#include "rbm.hpp"
#include <gtest/gtest.h>

using namespace std;
using namespace gpw::blat;
using namespace gpw::geometry;

TEST(SO3, TrivialConstruction) {
  SO3 I{{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};

  vec3 w{0., 0., 0.};
  SO3 R1{w, 0.};

  for (std::size_t i = 1; i < 4; ++i)
    for (std::size_t j = 1; j < 4; ++j)
      EXPECT_NEAR(R1(i, j), I(i, j), 0.0001);

  SO3 R2{w};

  for (std::size_t i = 1; i < 4; ++i)
    for (std::size_t j = 1; j < 4; ++j)
      EXPECT_NEAR(R2(i, j), I(i, j), 0.0001);

  auto R3 = SO3{1., 2., 3.};
  SO3 R_approx{
      {-.6949, .7135, .0893}, {-.1920, -.3038, .9332}, {.6930, .6313, .3481}};

  for (std::size_t i = 1; i < 4; ++i)
    for (std::size_t j = 1; j < 4; ++j)
      EXPECT_NEAR(R3(i, j), R_approx(i, j), 0.0001);
}

TEST(SO3, TransposeInverse) {
  SO3 I{{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
  SO3 I_inv = inv(I);

  for (std::size_t i = 1; i < 4; ++i)
    for (std::size_t j = 1; j < 4; ++j)
      EXPECT_NEAR(I_inv(i, j), I(i, j), 0.0001);

  SO3 R{1., 2., 3.};
  // -.6949, .7135, .0893
  // -.1920, -.3038, .9332
  // .6930, .6313, .3481

  SO3 R1 = R.transpose();
  SO3 R2 = transpose(R);

  for (std::size_t i = 1; i < 4; ++i)
    for (std::size_t j = 1; j < 4; ++j) {
      EXPECT_NEAR(R1(i, j), R(j, i), 0.0001);
      EXPECT_NEAR(R1(i, j), R2(i, j), 0.0001);
    }

  SO3 I_3 = R * inv(R);

  for (std::size_t i = 1; i < 4; ++i)
    for (std::size_t j = 1; j < 4; ++j)
      EXPECT_NEAR(I_3(i, j), I(i, j), 0.0001);
}

TEST(SO3, ExponentialLogarithm) {
  auto zero = mat3{};
  vec3 w{1., 2., 3.};
  auto W = skew(w);
  EXPECT_EQ(W + transpose(W), zero);

  auto R = expm(w);
  SO3 R_approx{
      {-.6949, .7135, .0893}, {-.1920, -.3038, .9332}, {.6930, .6313, .3481}};
  for (std::size_t i = 1; i < 4; ++i)
    for (std::size_t j = 1; j < 4; ++j)
      EXPECT_NEAR(R(i, j), R_approx(i, j), 0.0001);

  auto log_R = logm(R);
  vec3 log_R_approx = {-.6793, -1.3585, -2.0378};
  for (std::size_t i = 1; i < 4; ++i)
    EXPECT_NEAR(log_R(i), log_R_approx(i), 0.0001);

  auto exp_log_R = expm(log_R);
  for (std::size_t i = 1; i < 4; ++i)
    for (std::size_t j = 1; j < 4; ++j)
      EXPECT_FLOAT_EQ(R(i, j), exp_log_R(i, j));
}

TEST(SE3, Creation) {
  auto T = SE3{vec3{1., 2., 3.}, vec3{3., 2., 1.}};
  SE3 T_expected{SO3{{-.6949, .7135, .0893},
                     {-.1920, -.3038, .9332},
                     {.6930, .6313, .3481}},
                 vec3{-.1522, 2.3854, 1.7938}};

  for (std::size_t i = 1; i < 4; ++i) {
    EXPECT_NEAR(T.p(i), T_expected.p(i), 0.0001);
    for (std::size_t j = 1; j < 4; ++j)
      EXPECT_NEAR(T.R(i, j), T_expected.R(i, j), 0.0001);
  }
}

TEST(SE3, Inverse) {
  auto T = SE3{vec3{1., 2., 3.}, vec3{3., 2., 1.}};
  auto Tinv = inv(T);
  auto I_expected = T * Tinv;
  auto I3x3 = SO3{{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1}};

  for (std::size_t i = 1; i < 4; ++i) {
    for (std::size_t j = 1; j < 4; ++j) {
      EXPECT_NEAR(I_expected.R(i, j), I3x3(i, j), 0.0001);
    }
    EXPECT_NEAR(I_expected.p(i), 0., 0.0001);
  }
}

TEST(SO3, Projection) {
  mat3 M{};
  mat3 error{{0.0001, 0.0002, -0.0001},
             {-0.0001, -0.0002, 0.0002},
             {0.0002, 0.0001, -0.0001}};
  M += error;

  EXPECT_NE(det(M), 1.);

  SO3 R{M};
  EXPECT_FLOAT_EQ(det(R), 1.);
}
