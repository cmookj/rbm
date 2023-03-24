cmake_minimum_required(VERSION 3.16..3.24)
project(RigidBodyMotion VERSION 0.0.1
    DESCRIPTION "Lie group class library for rigid body motion"
   LANGUAGES CXX)

add_library(rbm STATIC
    ${PROJECT_SOURCE_DIR}/src/rbm.cpp)

target_include_directories(rbm PUBLIC
    ${PROJECT_SOURCE_DIR}/include
    /usr/local/include/lapack
    /usr/local/include
    )

