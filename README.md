# rbm: Rigid Body Motion Library
A library to handle rigid body motion in 3D space using Lie Group.

## Dependancies
* This library depends on [vma](https://github.com/cmookj/vma) for basic linear algebra.  The linear algebra routines depend on 
  - Apple's Accelerate framework on macOS, iOS, etc.
  - libopenblas and libgfortran on other OSes, behind the scene.
* For unit test, this repository automatically fetches Google's GTest toolkit.

## Basic Functionalities
To represent the orientation in 3D space, this library implements SO(3) (_Special Orthogonal_) Group.
The orientation and position are represented as an element in SE(3) (_Special Euclidean_) Group, which can be considered as the combination of the SO(3) and 3-dimensional vector space.

The following functions are available in the library:
* Mapping between so(3) and SO(3), i.e., matrix exponential on so(3) and logarithm on SO(3).
* Generation of the rotation matrix (an element in SO(3) space) from Euler angle representations.
* Inverse of the elements in SO(3).
* Interpolation of orientation using a similar methodology as Bezier curves.
* Mapping between se(3) and SE(3), i.e., matrix exponential on se(3) and logarithm on SE(3).
* Generation of the SE(3) elements from SO(3) and a 3-dimensional vector.
* Combination of mappings in SE(3), inverse mapping, etc.
