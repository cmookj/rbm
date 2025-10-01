cmake_minimum_required(VERSION 3.16)

project(rbm VERSION 0.0.1
    DESCRIPTION "Lie group class library for rigid body motion"
    LANGUAGES CXX)

# GoogleTest requires at least C++17
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Dependencies
include(Dependencies.cmake)

# Generate compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_subdirectory(rbm)

enable_testing()

if(APPLE)
    find_library(ACCELERATE_LIBRARY Accelerate)
    set(extra_libs ${ACCELERATE_LIBRARY})
elseif(UNIX AND NOT APPLE)
    set(extra_libs
        libopenblas.a
        libgfortran.a # libopenblas requires link to _gfortran_concat_string
        )
endif()

add_executable(rbm_test
  test.cpp
)

target_link_libraries(rbm_test PRIVATE
  rbm
  vma
  toolbox
    ${extra_libs}
  GTest::gtest_main
)

target_include_directories(rbm_test PRIVATE
  ${PROJECT_SOURCE_DIR}/rbm/include
)

include(GoogleTest)
gtest_discover_tests(rbm_test)
