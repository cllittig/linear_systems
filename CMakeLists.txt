cmake_minimum_required(VERSION 3.14)
project(LinearSystems)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Find BLAS and LAPACK
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Add project libraries
file(GLOB_RECURSE SOURCES_ALGEBRA "src/algebra_linear/*.cpp")
file(GLOB_RECURSE SOURCES_METODOS "src/metodos/*.cpp")
file(GLOB_RECURSE SOURCES_UTILS "src/utils/*.cpp")

# Main executable
file(GLOB_RECURSE SOURCES_MAIN "src/main.cpp")
add_executable(programa_tcc 
    ${SOURCES_MAIN} 
    ${SOURCES_ALGEBRA} 
    ${SOURCES_METODOS} 
    ${SOURCES_UTILS}
    src/algebra_linear/matriz_identity.cpp
)
target_link_libraries(programa_tcc PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} m pthread)

# Google Test configuration
option(BUILD_TESTS "Build tests" ON)

if(BUILD_TESTS)
    # Download and build Google Test
    include(FetchContent)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1
    )
    
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)

    enable_testing()

    # Test executable
    add_executable(programa_testes
        tests/test_matriz.cpp
        tests/test_matriz_advanced.cpp
        tests/test_vector.cpp
        tests/test_vector_advanced.cpp
        tests/test_conjugate_gradient.cpp
        tests/test_compare_with_octave.cpp
        tests/test_compare_all_methods.cpp
        tests/test_numerical_methods.cpp
        src/algebra_linear/matriz_identity.cpp
        ${SOURCES_ALGEBRA}
        ${SOURCES_METODOS}
        ${SOURCES_UTILS}
    )

    target_link_libraries(programa_testes PRIVATE 
        gtest_main 
        gtest 
        ${BLAS_LIBRARIES} 
        ${LAPACK_LIBRARIES} 
        m pthread
    )

    add_test(NAME LinearSystemsTests COMMAND programa_testes)
endif()

add_executable(bench_compare
    src/bench/bench_compare_octave.cpp
    ${SOURCES_ALGEBRA}
    ${SOURCES_METODOS}
    ${SOURCES_UTILS}
)

target_link_libraries(bench_compare PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} m pthread)

add_executable(bench_fair
    src/bench/bench_fair_comparison.cpp
    ${SOURCES_ALGEBRA}
    ${SOURCES_METODOS}
    ${SOURCES_UTILS}
)

target_link_libraries(bench_fair PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} m pthread)
