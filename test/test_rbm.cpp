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

TEST (SO3, TrivialConstruction) {
    SO3 I {
        {1., 0., 0.},
        {0., 1., 0.},
        {0., 0., 1.}
    };

    vec3 w {0., 0., 0.};
    SO3 R1 {w, 0.};

    for (std::size_t i = 1; i < 4; ++i)
        for (std::size_t j = 1; j < 4; ++j)
            EXPECT_NEAR (R1 (i, j), I (i, j), 0.0001);

    SO3 R2 {w};

    for (std::size_t i = 1; i < 4; ++i)
        for (std::size_t j = 1; j < 4; ++j)
            EXPECT_NEAR (R2 (i, j), I (i, j), 0.0001);

    auto R3 = SO3 {1., 2., 3.};
    SO3 R_approx {
        {-.6949,  .7135, .0893},
        {-.1920, -.3038, .9332},
        { .6930,  .6313, .3481}
    };

    for (std::size_t i = 1; i < 4; ++i)
        for (std::size_t j = 1; j < 4; ++j)
            EXPECT_NEAR (R3 (i, j), R_approx (i, j), 0.0001);
}

TEST (SO3, Construction) {
    auto R1 = SO3::euler_zyx (90.0, 0.0, 0.0);
    double tol = 1e-12;

    EXPECT_NEAR (R1 (1, 1), 0., tol);
    EXPECT_NEAR (R1 (2, 1), 1., tol);
    EXPECT_NEAR (R1 (3, 1), 0., tol);

    EXPECT_NEAR (R1 (1, 2), -1., tol);
    EXPECT_NEAR (R1 (2, 2), 0., tol);
    EXPECT_NEAR (R1 (3, 2), 0., tol);

    EXPECT_NEAR (R1 (1, 3), 0., tol);
    EXPECT_NEAR (R1 (2, 3), 0., tol);
    EXPECT_NEAR (R1 (3, 3), 1., tol);
}

TEST (SO3, Multiplication) {
    SO3 R1 {1, 3, 5};

    EXPECT_NEAR (R1 (1, 1), 0.93527384, 0.000001);
    EXPECT_NEAR (R1 (1, 2), 0.30904994, 0.000001);
    EXPECT_NEAR (R1 (1, 3), -0.17248473, 0.000001);
    EXPECT_NEAR (R1 (2, 1), -0.29762768, 0.000001);
    EXPECT_NEAR (R1 (2, 2), 0.95050352, 0.000001);
    EXPECT_NEAR (R1 (2, 3), 0.08922342, 0.000001);
    EXPECT_NEAR (R1 (3, 1), 0.19152184, 0.000001);
    EXPECT_NEAR (R1 (3, 2), -0.0321121, 0.000001);
    EXPECT_NEAR (R1 (3, 3), 0.98096289, 0.000001);

    SO3 R2 {-1, -2, -3};

    EXPECT_NEAR (R2 (1, 1), -0.69492056, 0.000001);
    EXPECT_NEAR (R2 (1, 2), -0.19200697, 0.000001);
    EXPECT_NEAR (R2 (1, 3), 0.69297817, 0.000001);
    EXPECT_NEAR (R2 (2, 1), 0.71352099, 0.000001);
    EXPECT_NEAR (R2 (2, 2), -0.30378504, 0.000001);
    EXPECT_NEAR (R2 (2, 3), 0.6313497, 0.000001);
    EXPECT_NEAR (R2 (3, 1), 0.08929286, 0.000001);
    EXPECT_NEAR (R2 (3, 2), 0.93319235, 0.000001);
    EXPECT_NEAR (R2 (3, 3), 0.34810748, 0.000001);

    SO3 R12 = R1 * R2;

    EXPECT_NEAR (R12 (1, 1), -0.44482905, 0.000001);
    EXPECT_NEAR (R12 (1, 2), -0.43442528, 0.000001);
    EXPECT_NEAR (R12 (1, 3), 0.78319971, 0.000001);

    EXPECT_NEAR (R12 (2, 1), 0.89299882, 0.000001);
    EXPECT_NEAR (R12 (2, 2), -0.14833955, 0.000001);
    EXPECT_NEAR (R12 (2, 3), 0.42490997, 0.000001);

    EXPECT_NEAR (R12 (3, 1), -0.06841214, 0.000001);
    EXPECT_NEAR (R12 (3, 2), 0.88840872, 0.000001);
    EXPECT_NEAR (R12 (3, 3), 0.45392701, 0.000001);
}

TEST (SO3, TransposeInverse) {
    SO3 I {
        {1., 0., 0.},
        {0., 1., 0.},
        {0., 0., 1.}
    };
    SO3 I_inv = inv (I);

    for (std::size_t i = 1; i < 4; ++i)
        for (std::size_t j = 1; j < 4; ++j)
            EXPECT_NEAR (I_inv (i, j), I (i, j), 0.0001);

    SO3 R {1., 2., 3.};
    // -.6949, .7135, .0893
    // -.1920, -.3038, .9332
    // .6930, .6313, .3481

    SO3 R1 = R.transpose();
    SO3 R2 = transpose (R);

    for (std::size_t i = 1; i < 4; ++i)
        for (std::size_t j = 1; j < 4; ++j) {
            EXPECT_NEAR (R1 (i, j), R (j, i), 0.0001);
            EXPECT_NEAR (R1 (i, j), R2 (i, j), 0.0001);
        }

    SO3 R_inv = inv (R);
    for (std::size_t i = 1; i < 4; ++i)
        for (std::size_t j = 1; j < 4; ++j)
            EXPECT_NEAR (R (i, j), R_inv (j, i), 0.0001);

    SO3 I_3 = R * R_inv;

    for (std::size_t i = 1; i < 4; ++i)
        for (std::size_t j = 1; j < 4; ++j)
            EXPECT_NEAR (I_3 (i, j), I (i, j), 0.0001);
}

TEST (SO3, RelativeOrientation) {
    SO3 R1 {1, 3, 5};
    SO3 R {-1, -2, -3};
    SO3 R2 = R1 * R;
    SO3 R_calc = relative (R1, R2);

    for (std::size_t i = 1; i < 4; ++i)
        for (std::size_t j = 1; j < 4; ++j)
            EXPECT_NEAR (R (i, j), R_calc (i, j), 0.00001);
}

TEST (SO3, Projection) {
    SO3 R1_org {1, 3, 5};
    mat3 M1 {};
    mat3 error1 {
        { 0.0000001,  0.0000002, -0.0000001},
        {-0.0000001, -0.0000002,  0.0000002},
        { 0.0000002,  0.0000001, -0.0000001}
    };
    for (std::size_t i = 1; i < 4; ++i)
        for (std::size_t j = 1; j < 4; ++j)
            M1 (i, j) = R1_org (i, j) + error1 (i, j);

    EXPECT_NE (det (M1), 1.);

    SO3 R1 {M1};
    EXPECT_TRUE (std::abs (det (R1) - 1.) < std::abs (det (M1) - 1.));

    for (std::size_t i = 1; i < 4; ++i)
        for (std::size_t j = 1; j < 4; ++j)
            EXPECT_NEAR (R1_org (i, j), R1 (i, j), 0.00001);

    SO3 R2_org {3, 5, 7};
    // -0.73944516  0.11502181  0.66331806
    //  0.59015865 -0.36334891  0.72089551
    //  0.3239346   0.92452558  0.20079547
    mat3 M2 {};
    mat3 error2 {
        { 0.0000001, -0.0000001, -0.0000001},
        {-0.0000001,  0.0000001,  0.0000001},
        { 0.0000001,  0.0000001, -0.0000001}
    };
    for (std::size_t i = 1; i < 4; ++i)
        for (std::size_t j = 1; j < 4; ++j)
            M2 (i, j) = R2_org (i, j) + error2 (i, j);

    EXPECT_NE (det (M2), 1.);
    SO3 R2 {M2};
    EXPECT_TRUE (std::abs (det (R2) - 1.) < std::abs (det (M2) - 1.));

    for (std::size_t i = 1; i < 4; ++i)
        for (std::size_t j = 1; j < 4; ++j)
            EXPECT_NEAR (R2_org (i, j), R2 (i, j), 0.00001);
}

TEST (SO3, ExponentialLogarithm) {
    auto zero = mat3 {};
    vec3 w {1., 2., 3.};
    auto W = skew (w);
    EXPECT_EQ (W + transpose (W), zero);

    auto R = expm (w);
    SO3 R_approx {
        {-.6949,  .7135, .0893},
        {-.1920, -.3038, .9332},
        { .6930,  .6313, .3481}
    };
    for (std::size_t i = 1; i < 4; ++i)
        for (std::size_t j = 1; j < 4; ++j)
            EXPECT_NEAR (R (i, j), R_approx (i, j), 0.0001);

    auto log_R = logm (R);
    vec3 log_R_approx = {-.6793, -1.3585, -2.0378};
    for (std::size_t i = 1; i < 4; ++i)
        EXPECT_NEAR (log_R (i), log_R_approx (i), 0.0001);

    auto exp_log_R = expm (log_R);
    for (std::size_t i = 1; i < 4; ++i)
        for (std::size_t j = 1; j < 4; ++j)
            EXPECT_FLOAT_EQ (R (i, j), exp_log_R (i, j));

    SO3 I {}; // Identity
    auto log_I = logm (I);
    EXPECT_FLOAT_EQ (norm (log_I), 0.0);
    for (std::size_t i = 1; i < 4; ++i)
        EXPECT_FLOAT_EQ (log_I (i), 0.);
}

TEST (SE3, Creation) {
    auto T = SE3 {
        vec3 {1., 2., 3.},
         vec3 {3., 2., 1.}
    };
    SE3 T_expected {
        SO3 {
             {-.6949, .7135, .0893},
             {-.1920, -.3038, .9332},
             {.6930, .6313, .3481}                 },
        vec3 {                -.1522, 2.3854, 1.7938}
    };

    for (std::size_t i = 1; i < 4; ++i) {
        EXPECT_NEAR (T.p (i), T_expected.p (i), 0.0001);
        for (std::size_t j = 1; j < 4; ++j)
            EXPECT_NEAR (T.R (i, j), T_expected.R (i, j), 0.0001);
    }
<<<<<<< HEAD

  SO3 I_3 = R * inv(R);

  for (std::size_t i = 1; i < 4; ++i)
    for (std::size_t j = 1; j < 4; ++j)
      EXPECT_NEAR(I_3(i, j), I(i, j), 0.0001);
=======
>>>>>>> 61ddbbf24edf19187979778f31bf58d29af5656c
}

TEST (SE3, Inverse) {
    auto T = SE3 {
        vec3 {1., 2., 3.},
         vec3 {3., 2., 1.}
    };
    auto Tinv = inv (T);
    auto I_expected = T * Tinv;
    auto I3x3 = SO3 {
        {1., 0., 0.},
        {0., 1., 0.},
        {0., 0.,  1}
    };

    for (std::size_t i = 1; i < 4; ++i) {
        for (std::size_t j = 1; j < 4; ++j) {
            EXPECT_NEAR (I_expected.R (i, j), I3x3 (i, j), 0.0001);
        }
        EXPECT_NEAR (I_expected.p (i), 0., 0.0001);
    }
}

TEST (SE3, Multiplication) {
    // First test
    vec3 w1 {1, 2, 3};
    SO3 R1 = SO3 {w1};
    vec3 p1 {6, 7, 8};
    SE3 T1 {R1, p1};

    EXPECT_NEAR (T1.R (1, 1), -0.69492056, 0.000001);
    EXPECT_NEAR (T1.R (1, 2), 0.71352099, 0.000001);
    EXPECT_NEAR (T1.R (1, 3), 0.08929286, 0.000001);
    EXPECT_NEAR (T1.R (2, 1), -0.19200697, 0.000001);
    EXPECT_NEAR (T1.R (2, 2), -0.30378504, 0.000001);
    EXPECT_NEAR (T1.R (2, 3), 0.93319235, 0.000001);
    EXPECT_NEAR (T1.R (3, 1), 0.69297817, 0.000001);
    EXPECT_NEAR (T1.R (3, 2), 0.6313497, 0.000001);
    EXPECT_NEAR (T1.R (3, 3), 0.34810748, 0.000001);

    vec3 w2 {-3, -2, -1};
    SO3 R2 = SO3 {w2};
    vec3 p2 {5, 4, 3};
    SE3 T2 {R2, p2};

    EXPECT_NEAR (T2.R (1, 1), 0.34810748, 0.000001);
    EXPECT_NEAR (T2.R (1, 2), 0.6313497, 0.000001);
    EXPECT_NEAR (T2.R (1, 3), 0.69297817, 0.000001);
    EXPECT_NEAR (T2.R (2, 1), 0.93319235, 0.000001);
    EXPECT_NEAR (T2.R (2, 2), -0.30378504, 0.000001);
    EXPECT_NEAR (T2.R (2, 3), -0.19200697, 0.000001);
    EXPECT_NEAR (T2.R (3, 1), 0.08929286, 0.000001);
    EXPECT_NEAR (T2.R (3, 2), 0.71352099, 0.000001);
    EXPECT_NEAR (T2.R (3, 3), -0.69492056, 0.000001);

    SE3 T = T1 * T2;

    EXPECT_NEAR (T.R (1, 1), 0.4319185, 0.000001);
    EXPECT_NEAR (T.R (1, 2), -0.59178256, 0.000001);
    EXPECT_NEAR (T.R (1, 3), -0.68061722, 0.000001);
    EXPECT_NEAR (T.R (2, 1), -0.26700153, 0.000001);
    EXPECT_NEAR (T.R (2, 2), 0.63691414, 0.000001);
    EXPECT_NEAR (T.R (2, 3), -0.72322234, 0.000001);
    EXPECT_NEAR (T.R (3, 1), 0.86148511, 0.000001);
    EXPECT_NEAR (T.R (3, 2), 0.49409895, 0.000001);
    EXPECT_NEAR (T.R (3, 3), 0.11708815, 0.000001);

    EXPECT_NEAR (T.p (1), 5.64735975, 0.000001);
    EXPECT_NEAR (T.p (2), 7.62440202, 0.000001);
    EXPECT_NEAR (T.p (3), 15.03461207, 0.000001);

    // Second test
    vec3 w3 {2, -6, 4};
    SO3 R3 = SO3 {w3};
    vec3 p3 {3, 5, 7};
    SE3 T3 {R3, p3};

    EXPECT_NEAR (T3.R (1, 1), 0.40779158, 0.000001);
    EXPECT_NEAR (T3.R (1, 2), -0.6348844, 0.000001);
    EXPECT_NEAR (T3.R (1, 3), -0.65622239, 0.000001);
    EXPECT_NEAR (T3.R (2, 1), 0.36155744, 0.000001);
    EXPECT_NEAR (T3.R (2, 2), 0.77222753, 0.000001);
    EXPECT_NEAR (T3.R (2, 3), -0.52243742, 0.000001);
    EXPECT_NEAR (T3.R (3, 1), 0.83844037, 0.000001);
    EXPECT_NEAR (T3.R (3, 2), -0.0242165, 0.000001);
    EXPECT_NEAR (T3.R (3, 3), 0.54445506, 0.000001);

    vec3 w4 {-9, 7, 5};
    SO3 R4 = SO3 {w4};
    vec3 p4 {-3, 6, -9};
    SE3 T4 {R4, p4};

    EXPECT_NEAR (T4.R (1, 1), 0.99676544, 0.000001);
    EXPECT_NEAR (T4.R (1, 2), 0.04391646, 0.000001);
    EXPECT_NEAR (T4.R (1, 3), -0.06730524, 0.000001);
    EXPECT_NEAR (T4.R (2, 1), -0.04942395, 0.000001);
    EXPECT_NEAR (T4.R (2, 2), 0.99536672, 0.000001);
    EXPECT_NEAR (T4.R (2, 3), -0.0824765, 0.000001);
    EXPECT_NEAR (T4.R (3, 1), 0.06337132, 0.000001);
    EXPECT_NEAR (T4.R (3, 2), 0.08553622, 0.000001);
    EXPECT_NEAR (T4.R (3, 3), 0.99431767, 0.000001);

    T = T3 * T4;

    EXPECT_NEAR (T.R (1, 1), 0.39626536, 0.000001);
    EXPECT_NEAR (T.R (1, 2), -0.67016482, 0.000001);
    EXPECT_NEAR (T.R (1, 3), -0.62757698, 0.000001);
    EXPECT_NEAR (T.R (2, 1), 0.28911388, 0.000001);
    EXPECT_NEAR (T.R (2, 2), 0.73984058, 0.000001);
    EXPECT_NEAR (T.R (2, 3), -0.6074941, 0.000001);
    EXPECT_NEAR (T.R (3, 1), 0.8714281, 0.000001);
    EXPECT_NEAR (T.R (3, 2), 0.05928766, 0.000001);
    EXPECT_NEAR (T.R (3, 3), 0.48692715, 0.000001);

    EXPECT_NEAR (T.p (1), 3.87332038, 0.000001);
    EXPECT_NEAR (T.p (2), 13.25062968, 0.000001);
    EXPECT_NEAR (T.p (3), -0.56071567, 0.000001);
}

TEST (SE3, CoordinateTransformation) {
    auto T_01 = SE3 {
        SO3 {},
         vec3 {1.0, 2.0, 3.0}
    };
    auto p_1 = vec3 {1.0, 2.0, 3.0};
    auto p_0 = T_01 * p_1;

    EXPECT_FLOAT_EQ (p_0 (1), 2.0);
    EXPECT_FLOAT_EQ (p_0 (2), 4.0);
    EXPECT_FLOAT_EQ (p_0 (3), 6.0);

    auto T1 = SE3 {
        SO3 {2, -6, 4},
         vec3 {3,  5, 7}
    };

    auto T2 = SE3 {
        SO3 {-9, 7,  5},
         vec3 {-3, 6, -9}
    };

    auto T12 = relative (T1, T2);
    auto T22 = T1 * T12;

    EXPECT_TRUE (similar (T2, T22, gpw::geometry::EPS * 1e4));
    std::cout << "T2:\n";
    std::cout << to_string (T2) << '\n';

    std::cout << "T22:\n";
    std::cout << to_string (T22) << '\n';

    // clang-format off
    // 0.44173603   0.44950802   0.77640957 -15.5002379 
    // -0.67253199   0.73869632  -0.04503848   4.96899801
    // -0.59377605  -0.50226516   0.62861731  -5.29638404
    // 
    // 0.99676544  0.04391646 -0.06730524 -3.
    // -0.04942395  0.99536672 -0.0824765   6.
    // 0.06337132  0.08553622  0.99431767 -9.
    // 
    // 0.99676544  0.04391646 -0.06730524 -3.
    // -0.04942395  0.99536672 -0.0824765   6.
    // 0.06337132  0.08553622  0.99431767 -9.
    // clang-format on
}

TEST (SE3, RandomPose) {
    auto I3 = SO3 {};
    auto p3 = vec3 {};

    auto T0 = SE3 {
        vec3 {0., 0., 0.},
         vec3 {0., 0., 0.}
    };

    EXPECT_TRUE (similar (T0.R(), I3));
    EXPECT_TRUE (similar (T0.p(), p3));
}
