# BUILD_TESTING variable is created by include(CTest)
# It is set to ON by default.
if (BUILD_TESTING)
    add_executable(test_rbm test_rbm.cpp)

    if(APPLE)
        find_library(ACCELERATE_LIBRARY Accelerate)
        set(extra_libs ${ACCELERATE_LIBRARY})
    elseif(UNIX AND NOT APPLE)
        set(extra_libs 
            libopenblas.a
            libgfortran.a # libopenblas requires link to _gfortran_concat_string
            )
    endif()

    target_link_libraries(test_rbm PRIVATE 
        rbm 
        ${extra_libs}
        GTest::gtest_main)

    target_include_directories(test_rbm PRIVATE 
        ${PROJECT_SOURCE_DIR}/rbm/include
        /usr/local/include
        )

    include(GoogleTest)
    # Finds all the Google tests associated with the executable 
    gtest_discover_tests(test_rbm)
endif()
