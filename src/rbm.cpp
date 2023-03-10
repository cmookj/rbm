//===------------------------------------------------------------*- C++ -*-===//
//
//  rbm.cpp
//  RigidBodyMotionLibrary
//
//  Created by Changmook Chun on 11/23/22.
//
//===----------------------------------------------------------------------===//

#include "rbm.hpp"

#include <algorithm>
#include <cstring>
#include <exception>
#include <functional>
#include <limits>
#include <numeric>

#if defined(_WIN32) || defined(_WIN64)

#define NOMINMAX
const double M_PI   = 3.14159265358979323846;
const double M_PI_2 = 1.57079632679489661923;
const double M_PI_4 = 0.785398163397448309616;

#endif

using namespace std;
using namespace gpw::geometry;
using namespace gpw::blat;


#pragma mark - SO(3)

SO3::SO3()
: mat3 {{{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}}} {}

SO3::SO3(const std::initializer_list<std::initializer_list<double>>& l)
    : mat3 {l} {}

SO3::SO3(const array<double, 9>& el)
: mat3 {{el[0], el[1], el[2]}, {el[3], el[4], el[5]}, {el[6], el[7], el[8]}} {}

SO3::SO3(const double& w1, const double& w2, const double& w3) {
    _expm(w1, w2, w3);
}

SO3::SO3(const vec3& v, const double& q) {
    vec3 axis {v / norm(v) * q};
    _expm(axis(1), axis(2), axis(3));
}

SO3::SO3(const vec3& v) {
    _expm(v(1), v(2), v(3));
}

/**
 @brief Finds the closest SO(3) to a 3 x 3 matrix
 
 @details
 Projection of a general 3 x 3 matrix M onto 3-dimensional special orthogonal
 group is to find a matrix R in SO(3) which is the closest to M in the sense of
 Frobenius norm.  That is we want to find R s.t.
   J(R) = tr{(R - M)'(R - M)}
 is minimized, where M' is the transpose of M.  Simple manipulation of the
 objective function results in
   J(R) = -2 tr{RM'}.
 
 When M is invertible, the objective function has optimal solution in closed
 form:
   R = V diag(1, 1, 1) U',  if det(M) > 0;
       V diag(1, 1, -1) U', if det(M) < 0,
 where M = U S V' (singular value decomposition) and
                   [ a 0 0 ]
   diag(a, b, c) = [ 0 b 0 ].
                   [ 0 0 c ]
                  
 */
SO3::SO3(const mat3& M)
: mat3 {{{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}}} {
    double dtm = det(M);
    if (dtm != 0) {
        auto SVD = svd(M);

        mat3 I = gpw::blat::identity<3>();
        if (dtm < 0.)
            I(3, 3) = -1.;
        
        mat3 projection = SVD.V * I * gpw::blat::transpose(SVD.U);
        _elem = projection.elem();
    }
}

SO3 SO3::euler_zyx(double q1, double q2, double q3) {
    q1 *= M_PI / 180.;
    q2 *= M_PI / 180.;
    q3 *= M_PI / 180.;

    double s1 = sin(q1);
    double s2 = sin(q2);
    double s3 = sin(q3);

    double c1 = cos(q1);
    double c2 = cos(q2);
    double c3 = cos(q3);

    SO3 R;
    R._elem[0] = c1 * c2;
    R._elem[1] = c2 * s1;
    R._elem[2] = -s2;

    R._elem[3] = c1 * s2 * s3 - c3 * s1;
    R._elem[4] = c1 * c3 + s1 * s2 * s3;
    R._elem[5] = c2 * s3;

    R._elem[6] = s1 * s3 + c1 * c3 * s2;
    R._elem[7] = c3 * s1 * s2 - c1 * s3;
    R._elem[8] = c2 * c3;

    return R;
}

SO3 SO3::euler_zyz(double q1, double q2, double q3) {
    q1 *= M_PI / 180.;
    q2 *= M_PI / 180.;
    q3 *= M_PI / 180.;

    double s1 = sin(q1);
    double s2 = sin(q2);
    double s3 = sin(q3);

    double c1 = cos(q1);
    double c2 = cos(q2);
    double c3 = cos(q3);

    SO3 R;
    R._elem[0] = c1 * c2 * c3 - s1 * s3;
    R._elem[1] = c1 * s3 + c2 * c3 * s1;
    R._elem[2] = -c3 * s2;

    R._elem[3] = -c3 * s1 - c1 * c2 * s3;
    R._elem[4] = c1 * c3 - c2 * s1 * s3;
    R._elem[5] = s2 * s3;

    R._elem[6] = c1 * s2;
    R._elem[7] = s1 * s2;
    R._elem[8] = c2;

    return R;
}

SO3 SO3::transpose() const {
    SO3 R {{_elem[0], _elem[3], _elem[6]},
        {_elem[1], _elem[4], _elem[7]},
        {_elem[2], _elem[5], _elem[8]}};
    return R;
}

SO3 SO3::inv() const {
    SO3 R {{_elem[0], _elem[3], _elem[6]},
        {_elem[1], _elem[4], _elem[7]},
        {_elem[2], _elem[5], _elem[8]}};
    return R;
}

void SO3::_expm(double w1, double w2, double w3) {
    double theta = sqrt(w1 * w1 + w2 * w2 + w3 * w3);
    if (theta == 0.) {
        _elem[0] = 1.;
        _elem[1] = 0.;
        _elem[2] = 0.;
        _elem[3] = 0.;
        _elem[4] = 1.;
        _elem[5] = 0.;
        _elem[6] = 0.;
        _elem[7] = 0.;
        _elem[8] = 1.;
    } else {
        w1 /= theta;
        w2 /= theta;
        w3 /= theta;
        
        double a = sin(theta);
        double b = 1. - cos(theta);

        _elem[0] = 1. - b * (w2 * w2 + w3 * w3);
        _elem[1] = b * w1 * w2 + w3 * a;
        _elem[2] = b * w3 * w1 - w2 * a;

        _elem[3] = b * w1 * w2 - w3 * a;
        _elem[4] = 1. - b * (w3 * w3 + w1 * w1);
        _elem[5] = b * w2 * w3 + w1 * a;

        _elem[6] = b * w3 * w1 + w2 * a;
        _elem[7] = b * w2 * w3 - w1 * a;
        _elem[8] = 1. - b * (w1 * w1 + w2 * w2);
    }
}


#pragma mark - SE(3)

SE3::SE3() {
    attitude = SO3();
    position = vec3();
}

/**
 @brief Constructs SE3 from a rotation matrix and a position vector
 */
SE3::SE3(const SO3& R, const vec3& p) {
    attitude = R;
    position = p;
}

/**
 @brief Constructs SE3 using the matrix exponential of se3
 
 @details Both the screw [w] in so3 and position are given as 3-dimensional vectors.
 */
SE3::SE3(const vec3& w, const vec3& v) {
    _expm(w, v);
}

/**
 @brief Constructs SE3 using the matrix exponential of se3
 
 @details The screw and position (both of them are 3-dimensional vectors) are
 given as 6 scalars.
 */
SE3::SE3(
    const double& w1, const double& w2, const double& w3, const double& v1,
    const double& v2, const double& v3
) {
    _expm(vec3 {{w1, w2, w3}}, vec3 {{v1, v2, v3}});
}

const SO3& SE3::R() const {
    return attitude;
}

void SE3::set_R(const SO3& R) {
    attitude = R;
}

const vec3& SE3::p() const {
    return position;
}

void SE3::set_p(const vec3& p) {
    position = p;
}

const double& SE3::R(const std::size_t i, const std::size_t j) const {
    if (3 < i || 3 < j)
        throw std::runtime_error("Out of index");
    
    return attitude(i, j);
}

double& SE3::R(const std::size_t i, const std::size_t j) {
    return const_cast<double &>(static_cast<const SE3&>(*this).R(i, j));
}

const double& SE3::p(const std::size_t i) const {
    if (3 < i)
        throw std::runtime_error("Out of index");
    
    return position(i);
}

double& SE3::p(const std::size_t i) {
    return const_cast<double&>(static_cast<const SE3&>(*this).p(i));
}

SE3& SE3::operator*=(const SE3& T) {
    position += attitude * T.p();
    attitude *= T.R();

    return *this;
}

vec3 SE3::operator*(const vec3& v) const {
    return attitude * v + position;
}

void SE3::print(ostream& strm) const {
    strm.setf(ios_base::fixed, ios_base::floatfield);
    strm.precision(4);

    for (size_t i = 0; i < 3; ++i) {
        strm << "[  ";
        for (size_t j = 0; j < 3; ++j) {
            strm.width(10);
            strm << attitude(i + 1, j + 1) << "  ";
        }
        strm << "  ][  ";
        strm.width(10);
        strm << position(i + 1) << "  ]" << endl;
    }
    strm << endl;

    strm.setf(ios_base::fmtflags(0), ios_base::floatfield);
    strm.precision(6);
    strm.width(0);
}

void SE3::_expm(vec3 w, vec3 v) {
    double theta = norm(w);
    if (theta != 0.) {
        attitude = SO3 {w};
        w /= theta;
        v /= theta;
        double a = 1. - cos(theta);
        double b = theta - sin(theta);

        position = vec3 {
            {(theta - b * (w(2) * w(2) + w(3) * w(3))) * v(1) +
                 (b * w(1) * w(2) - a * w(3)) * v(2) +
                 (b * w(3) * w(1) + a * w(2)) * v(3),
             (b * w(1) * w(2) + a * w(3)) * v(1) +
                 (theta - b * (w(3) * w(3) + w(1) * w(1))) * v(2) +
                 (b * w(2) * w(3) - a * w(1)) * v(3),
             (b * w(3) * w(1) - a * w(2)) * v(1) +
                 (b * w(2) * w(3) + a * w(1)) * v(2) +
                 (theta - b * (w(1) * w(1) + w(2) * w(2))) * v(3)}};
    }
}

vector<double> SE3::elem() const {
    vector<double> elm(16, 0);

    elm[ 0] = attitude.elem()[0];
    elm[ 1] = attitude.elem()[1];
    elm[ 2] = attitude.elem()[2];
    elm[ 3] = 0.;
    
    elm[ 4] = attitude.elem()[3];
    elm[ 5] = attitude.elem()[4];
    elm[ 6] = attitude.elem()[5];
    elm[ 7] = 0.;
    
    elm[ 8] = attitude.elem()[6];
    elm[ 9] = attitude.elem()[7];
    elm[10] = attitude.elem()[8];
    elm[11] = 0.;
    
    elm[12] = position.elem()[0];
    elm[13] = position.elem()[1];
    elm[14] = position.elem()[2];
    elm[15] = 1.;
    
    return elm;
}


#pragma mark - OTHER FUNCTIONS

template <typename T>
bool gpw::geometry::similar(const T& a, const T& b) {
    if (fabs(a - b) < numeric_limits<T>::epsilon())
        return true;
    else
        return false;
}

template bool gpw::geometry::similar<float>(const float&, const float&);
template bool gpw::geometry::similar<double>(const double&, const double&);

vec3 gpw::geometry::cross(const vec3& a, const vec3& b) {
    return vec3 {
        {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
         a[0] * b[1] - a[1] * b[0]}};
}

mat3 gpw::geometry::skew(
    const double& v1, const double& v2, const double& v3
) {
    mat3 m {{0., v3, -v2}, {-v3, 0., v1}, {v2, -v1, 0.}};

    return m;
}

mat3 gpw::geometry::skew(const vec3& v) {
    return skew(v(1), v(2), v(3));
}

SO3 gpw::geometry::operator*(const SO3& R1, const SO3& R2) {
    auto R = gpw::blat::operator*(R1, R2);
    return SO3(move(R));
}

vec3 gpw::geometry::operator*(const SO3& R, const vec3& v) {
    auto Rv = gpw::blat::operator*(R, v);
    return vec3(move(Rv.elem()));
}

SO3 gpw::geometry::expm(const double w1, const double w2, const double w3) {
    return SO3(w1, w2, w3);
}

SO3 gpw::geometry::expm(const vec3& w) {
    return SO3(w(1), w(2), w(3));
}

vec3 gpw::geometry::logm(const SO3& R) {
    vec3   w;
    double theta = 0.;
    double trace = tr(R);
    if (trace == 3)
        return w;

    vec3 zero;
    if (trace == -1) {
        theta = M_PI;
        w(1)  = R(1, 3);
        w(2)  = R(2, 3);
        w(3)  = 1 + R(3, 3);
        if (w != zero) {
            w *= theta / sqrt(2 * (1 + R(3, 3)));
            return w;
        }

        w(1) = R(1, 2);
        w(2) = 1 + R(2, 2);
        w(3) = R(3, 2);
        if (w != zero) {
            w *= theta / sqrt(2 * (1 + R(2, 2)));
            return w;
        }

        w(1) = 1 + R(1, 1);
        w(2) = R(2, 1);
        w(3) = R(3, 1);
        w *= theta / sqrt(2 * (1 + R(1, 1)));
        return w;
    }

    theta = acos((trace - 1.) / 2.);
    w(1)  = -R(2, 3) + R(3, 2);
    w(2)  = R(1, 3) - R(3, 1);
    w(3)  = -R(1, 2) + R(2, 1);
    w *= theta / (2. * sin(theta));
    return w;
}

SO3 gpw::geometry::transpose(const SO3& R) {
    mat3 M {R};
    return SO3 {gpw::blat::transpose(M)};
}

SO3 gpw::geometry::inv(const SO3& R) {
    return transpose(R);
}

vec3 gpw::geometry::geodesic(const SO3& R1, const SO3& R2) {
    return logm(transpose(R1) * R2);
}

vector<vector<vec3>> gpw::geometry::interpolate(
    const vector<double>& t, const vector<SO3>& R, const vec3& omega,
    const vec3& alpha
) {
    size_t n = t.size() - 1;

    vector<vector<vec3>> abc;
    for (size_t i = 0; i != n + 1; i++) {
        vector<vec3> abc_i;
        abc_i.push_back(vec3());
        abc_i.push_back(vec3());
        abc_i.push_back(vec3());
        abc.push_back(abc_i);
    }

    for (size_t i = 1; i != n + 1; i++) {
        if (gpw::blat::tr(transpose(R[i - 1]) * R[i]) == -1)
            return abc;
    }

    vector<vec3> r;
    r.push_back(vec3());
    vector<mat3> A;
    A.push_back(mat3());
    for (size_t i = 1; i != n + 1; i++) {
        mat3 I   = gpw::blat::identity<3>();
        vec3 r_i = logm(transpose(R[i - 1]) * R[i]);
        r.push_back(r_i);

        double abs_r_i   = norm(r_i);
        double abs_r_i_2 = abs_r_i * abs_r_i;
        double abs_r_i_3 = abs_r_i_2 * abs_r_i;
        mat3   brk_r_i   = skew(r_i);
        mat3   A_i =
            I - (1. - cos(abs_r_i)) / abs_r_i_2 * brk_r_i +
            (abs_r_i - sin(abs_r_i)) / abs_r_i_3 * brk_r_i * brk_r_i;
        A.push_back(A_i);
    }

    abc[1][2] = omega;
    abc[1][1] = alpha / 2.;
    abc[1][0] = r[1] - abc[1][1] - abc[1][2];

    for (size_t i = 2; i != n + 1; i++) {
        vec3 s = r[i];
        vec3 t = 3 * abc[i - 1][0] + 2 * abc[i - 1][1] + abc[i - 1][2];
        vec3 u = 6 * abc[i - 1][0] + 2 * abc[i - 1][1];

        vec3 c_i  = A[i - 1] * abc[i - 1][2];
        abc[i][2] = c_i;

        double abs_s   = norm(s);
        double abs_s_2 = abs_s * abs_s;
        double abs_s_3 = abs_s_2 * abs_s;
        double abs_s_4 = abs_s_3 * abs_s;
        double abs_s_5 = abs_s_4 * abs_s;
        double sin_s   = sin(abs_s);
        double cos_s   = cos(abs_s);

        vec3 b_i =
            (u -
             inner(s, t) / abs_s_4 * (2 * cos_s + abs_s * sin_s - 2.) *
                 cross(s, t) -
             (1. - cos_s) / abs_s_2 * cross(s, u) +
             inner(s, t) / abs_s_5 * (3 * sin_s - abs_s * cos_s - 2. * abs_s) *
                 cross(s, cross(s, t)) +
             (abs_s - sin_s) / abs_s_3 *
                 (cross(t, cross(s, t)) + cross(s, cross(s, u)))) /
            2.;
        abc[i][1] = b_i;

        abc[i][0] = s - b_i - c_i;
    }

    return abc;
}

SO3 gpw::geometry::rotation_interpolated(
    const double t, const vector<double>& ts, const vector<SO3>& R,
    const vector<vector<vec3>>& abc
) {
    if ((t < ts[0]) || (ts.back() < t))
        return SO3();

    size_t i {0};

    for (size_t j = 1; j != ts.size() + 1; j++) {
        if ((ts[j - 1] <= t) && (t <= ts[j])) {
            i = j;
            break;
        }
    }

    double tau = (t - ts[i - 1]) / (ts[i] - ts[i - 1]);
    return R[i - 1] * expm(
                          abc[i][0] * tau * tau * tau + abc[i][1] * tau * tau +
                          abc[i][2] * tau
                      );
}

SE3 gpw::geometry::expm(const double w1, const double w2, const double w3,
                        const double v1, const double v2, const double v3) {
    return SE3 {w1, w2, w3, v1, v2, v3};
}

SE3 gpw::geometry::expm(const vec3& w, const vec3& v) {
    return SE3 {w, v};
}

SE3 gpw::geometry::operator*(SE3 a, SE3 b) {
    return a *= b;
}

bool gpw::geometry::operator==(const SE3& T, const gpw::blat::mat<4, 4>& M) {
    return (T.elem() == M.elem());
}

bool gpw::geometry::operator==(const gpw::blat::mat<4, 4>& M, const SE3& T) {
    return T == M;
}

bool gpw::geometry::operator!=(const SE3& T, const gpw::blat::mat<4, 4>& M) {
    return !(T == M);
}

bool gpw::geometry::operator!=(const gpw::blat::mat<4, 4>& M, const SE3& T) {
    return !(M == T);
}
