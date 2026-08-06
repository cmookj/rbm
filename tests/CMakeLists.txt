set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

enable_testing()

add_executable(rbm-test
    test.cpp
)

if(APPLE)
    find_library(ACCELERATE_LIBRARY Accelerate)
    set(extra_libs ${ACCELERATE_LIBRARY})
elseif(UNIX AND NOT APPLE)
    set(extra_libs
        libopenblas.a
        libgfortran.a # libopenblas requires link to _gfortran_concat_string
        )
endif()

target_link_libraries(rbm-test PRIVATE
  rbm
  vma
  cch-utils
  ${extra_libs}
  GTest::gtest_main
)

target_compile_options(rbm-test PRIVATE
    $<$<CONFIG:Debug>:-g -Wall -Wextra -Wpedantic>
    $<$<CONFIG:Release>:-O3 -Wall -Wextra -Wpedantic>
)

target_include_directories(rbm-test PRIVATE
  ${PROJECT_SOURCE_DIR}/include
)

include(GoogleTest)
gtest_discover_tests(rbm-test)
