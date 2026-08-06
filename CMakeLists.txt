cmake_minimum_required(VERSION 3.25)
project(rbm VERSION 0.0.1
    DESCRIPTION "Lie group class library for rigid body motion"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Generate compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Dependencies
include(Dependencies.cmake)

# Main library target
# ------------------------------------------------------------------------------
add_library(rbm STATIC
    ${PROJECT_SOURCE_DIR}/src/rbm.cpp)

target_include_directories(rbm PUBLIC
    ${PROJECT_SOURCE_DIR}/include)

target_link_libraries(rbm PRIVATE
    vma
    cch-utils
)

# Test target
# ------------------------------------------------------------------------------
add_subdirectory(tests)
