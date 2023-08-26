//===----------------------------------------------------------*- C++ -*-===//
//
//  rbm.hpp
//  RigidBodyMotionLibrary
//
//  Created by Changmook Chun on 11/23/22.
//
//===--------------------------------------------------------------------===//

#ifndef RIGID_BODY_MOTION
#define RIGID_BODY_MOTION

#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include "blat/blat.hpp"

namespace gpw::geometry {

const double EPS = 2.2204e-16;

using vec3 = gpw::blat::vec<3>;
using mat3 = gpw::blat::mat<3, 3>;

enum class output_fmt { sht, nml, ext, sci, scx };
int
set_format (std::stringstream& strm, output_fmt fmt = output_fmt::nml);

#pragma mark - SO(3)

/**
 @brief SO(3): Special Orthogonal Matrix in 3-Dimensional Space

 @details Note that the base type of this class stores the elements of a Matrix
 in column-major order.
 */
class SO3 : public mat3 {
public:
    // Default constructor and destructor
    SO3();
    virtual ~SO3() = default;

    // Copy constructor and assignment
    SO3 (const SO3&) = default;
    SO3&
    operator= (const SO3&) = default;

    // Move constructor and assignment
    SO3 (SO3&&) noexcept = default;
    SO3&
    operator= (SO3&&) noexcept = default;

    // Other constructors
    SO3 (const std::initializer_list<std::initializer_list<double>>&);
    SO3 (const double&, const double&, const double&);
    SO3 (const vec3&, const double&);
    SO3 (const vec3&);
    SO3 (const mat3&);

    static SO3
    euler_zyx (double, double, double);
    static SO3
    euler_zyz (double, double, double);

    SO3
    transpose () const;
    SO3
    inv () const;

private:
    void
    _expm (double, double, double);
};

#pragma mark - SE(3)

/**
 @brief SE(3): Special Euclidean Matrix in 3-Dimensional Space
 */
class SE3 {
private:
    SO3 attitude;
    vec3 position;

public:
    // Default constructor and destructor
    SE3();
    virtual ~SE3() = default;

    // Copy constructor and assignment
    SE3 (const SE3&) = default;
    SE3&
    operator= (const SE3&) = default;

    // Move constructor and assignment
    SE3 (SE3&&) noexcept = default;
    SE3&
    operator= (SE3&&) noexcept = default;

    // Other constructors
    SE3 (const SO3&, const vec3&);
    SE3 (const vec3&, const vec3&);
    SE3 (const double&, const double&, const double&, const double&, const double&, const double&);

    const SO3&
    R () const;
    void
    set_R (const SO3&);
    const vec3&
    p () const;
    void
    set_p (const vec3&);

    const double&
    R (const std::size_t, const std::size_t) const;
    double&
    R (const std::size_t, const std::size_t);

    const double&
    p (const std::size_t) const;
    double&
    p (const std::size_t);

    SE3&
    operator*= (const SE3&);
    vec3
    operator* (const vec3&) const;

    void
    print (std::ostream& strm = std::cout) const;

    std::vector<double>
    elem () const;

private:
    void _expm (vec3, vec3);
};

#pragma mark - OTHER FUNCTIONS

std::string
to_string (const SE3&, output_fmt fmt = output_fmt::nml);

vec3
cross (const vec3&, const vec3&);

mat3
skew (const double&, const double&, const double&);
mat3
skew (const vec3&);

template <typename T>
bool
similar (const T&, const T&, const T tol = std::numeric_limits<T>::epsilon());

SO3
operator* (const SO3&, const SO3&);
vec3
operator* (const SO3&, const vec3&);

SO3
expm (const double, const double, const double);
SO3
expm (const vec3&);
vec3
logm (const SO3&);

SO3
transpose (const SO3&);
SO3
inv (const SO3&);

/*
 @brief Calculate relative orientation of the second w.r.t. the first.
 */
SO3
relative (const SO3&, const SO3&);

vec3
geodesic (const SO3&, const SO3&);

std::vector<std::vector<vec3>>
interpolate (const std::vector<double>&, const std::vector<SO3>&, const vec3&, const vec3&);

SO3
rotation_interpolated (const double, const std::vector<double>&, const std::vector<SO3>&, const std::vector<std::vector<vec3>>&);

SE3
expm (
    const double,
    const double,
    const double,
    const double,
    const double,
    const double
);
SE3
expm (const vec3&, const vec3&);

SE3 operator* (SE3, SE3);
SE3
inv (const SE3&);

/*
 @brief Calculate relative pose of the second w.r.t. the first.
 */
SE3
relative (const SE3&, const SE3&);

bool
operator== (const SE3&, const gpw::blat::mat<4, 4>&);

bool
operator== (const gpw::blat::mat<4, 4>&, const SE3&);

bool
operator!= (const SE3&, const gpw::blat::mat<4, 4>&);

bool
operator!= (const gpw::blat::mat<4, 4>&, const SE3&);

bool
operator== (const SE3&, const SE3&);

bool
operator!= (const SE3&, const SE3&);

bool
similar (
    const SO3&,
    const SO3&,
    const double tol = std::numeric_limits<double>::epsilon()
);

bool
similar (
    const SE3&,
    const SE3&,
    const double tol = std::numeric_limits<double>::epsilon()
);

} // namespace gpw::geometry

#endif
