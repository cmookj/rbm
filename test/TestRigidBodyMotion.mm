//===---------------------------------------------------------*- objC++ -*-===//
//
//  TestRigidBodyMotion.mm
//  TestRigidBodyMotionLibrary
//
//  Created by Changmook Chun on 11/23/22.
//
//===----------------------------------------------------------------------===//

#import <XCTest/XCTest.h>

#include "rbm.hpp"

using namespace std;
using namespace gpw::blat;
using namespace gpw::geometry;

@interface TestRigidBodyMotion : XCTestCase

@end

@implementation TestRigidBodyMotion

- (void)setUp {
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
}

- (void)testSO3 {
    auto zero = mat3 {};
    vec3 w {1., 2., 3.};
    auto W = skew(w);
    XCTAssert(W + transpose(W) == zero);
    
    auto R = expm(w);
    SO3 R_approx {{-.6949, .7135, .0893},
        {-.1920, -.3038, .9332},
        {.6930, .6313, .3481}};
    for (std::size_t i = 1; i < 4; ++i)
        for (std::size_t j = 1; j < 4; ++j)
            XCTAssert(std::fabs(R(i, j) - R_approx(i, j)) < 0.0001);
    
    auto log_R = logm(R);
    vec3 log_R_approx = {-.6793, -1.3585, -2.0378};
    for (std::size_t i = 1; i < 4; ++i)
        XCTAssert(std::fabs(log_R(i) - log_R_approx(i)) < 0.0001);
    
    auto exp_log_R = expm(log_R);
    for (std::size_t i = 1; i < 4; ++i)
        for (std::size_t j = 1; j < 4; ++j)
            XCTAssert(std::fabs(R(i, j) - exp_log_R(i, j)) < 0.0001);
}

- (void)testSE3 {
    auto T = SE3{vec3{1., 2., 3.}, vec3{3., 2., 1.}};
    SE3 T_expected {
        SO3{
            {-.6949,  .7135,  .0893},
            {-.1920, -.3038,  .9332},
            { .6930,  .6313,  .3481}},
        vec3{-.1522, 2.3854, 1.7938}
    };
    
    for (std::size_t i = 1; i < 4; ++i) {
        XCTAssert(std::fabs(T.p(i) - T_expected.p(i) < 0.0001));
        for (std::size_t j = 1; j < 4; ++j)
            XCTAssert(std::fabs(T.R(i, j) - T_expected.R(i, j)) < 0.0001);
    }
    
}

- (void)testProjection {
    mat3 M {};
    mat3 error {
        { 0.0001,  0.0002, -0.0001},
        {-0.0001, -0.0002,  0.0002},
        { 0.0002,  0.0001,  -0.0001}
    };
    M += error;
    
    XCTAssert(det(M) != 1.);
    
    SO3 R {M};
    XCTAssert(std::fabs(det(R) - 1.) < 0.00000001);
}

- (void)testPerformanceExample {
    // This is an example of a performance test case.
    [self measureBlock:^{
        // Put the code you want to measure the time of here.
    }];
}

@end
