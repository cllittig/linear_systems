/*
 * ==================================================================================
 * TEST FILE: test_octave_comparison.cpp
 * 
 * PURPOSE:
 * This file contains comparative tests that validate the linear algebra library
 * results against GNU Octave's native implementations using the liboctave C++ API.
 * 
 * STRICT REQUIREMENTS (as specified):
 * 1. OCTAVE CALLS MUST BE AT C++ CODE LEVEL (liboctave API) ONLY
 * 2. NO system(), popen(), exec() or any shell/terminal invocations allowed
 * 3. Data conversion between library format and liboctave format required
 * 4. Execution timing for both implementations
 * 5. Result comparison with EXPECT_NEAR
 * 
 * ==================================================================================
 * COMPILATION CONFIGURATION
 * ==================================================================================
 * 
 * To compile this test file, you MUST link against liboctave libraries.
 * Update your Makefile or build configuration with the following flags:
 * 
 * 1. INCLUDE PATHS:
 *    -I/usr/include/octave-X.X.X (where X.X.X is your Octave version)
 *    Example for Octave 8.3: -I/usr/include/octave-8.3.0
 * 
 * 2. LIBRARY PATHS:
 *    -L/usr/lib/octave-X.X.X
 * 
 * 3. LINKER FLAGS:
 *    -loctave (main Octave library)
 *    -loctinterp (Octave interpreter)
 *    -lamd -lcolamd -lccolamd -lcholmod -lumfpack -loct (sparse matrix libraries)
 *    -lm (math library)
 * 
 * 4. FULL COMPILATION COMMAND EXAMPLE:
 * 
 *    g++ -std=c++17 -O2 -Wall -Wextra \
 *        -I/usr/include/octave-8.3.0 \
 *        -Iinclude \
 *        -L/usr/lib/x86_64-linux-gnu \
 *        tests/test_octave_comparison.cpp \
 *        -o build/tests/test_octave_comparison \
 *        -loctave -loctinterp -lcblas -llapack -lm \
 *        libls.a \
 *        -gtest_main -lgtest
 * 
 * 5. FINDING OCTAVE DIRECTORIES:
 * 
 *    To find your Octave installation:
 *    $ octave-config --cppflags
 *    $ octave-config --octlibdir
 *    $ octave-config --libs
 * 
 * 6. FINDING VERSION:
 * 
 *    $ octave --version
 * 
 * 7. INSTALL OCTAVE DEVELOPMENT FILES (if not present):
 * 
 *    Debian/Ubuntu:
 *      sudo apt-get install liboctave-dev
 * 
 *    Fedora/RHEL:
 *      sudo dnf install octave-devel
 * 
 *    macOS (Homebrew):
 *      brew install octave
 * 
 * ==================================================================================
 * LIBOCTAVE API DOCUMENTATION
 * ==================================================================================
 * 
 * Key Octave classes used in these tests:
 * 
 * - octave::Matrix: Dense matrix class (double precision)
 * - octave::ColumnVector: Column vector class
 * - octave::RowVector: Row vector class
 * - F<method>(): Octave C++ API functions for various decompositions
 * 
 * Example conversion:
 * 
 *   // Convert from library Matrix to Octave Matrix
 *   Matriz lib_A(3, 3);
 *   // ... fill lib_A ...
 *   
 *   octave::Matrix oct_A(lib_A.getRows(), lib_A.getColumns());
 *   for (int i = 0; i < lib_A.getRows(); ++i) {
 *       for (int j = 0; j < lib_A.getColumns(); ++j) {
 *           oct_A(i, j) = lib_A.getValue(i, j);
 *       }
 *   }
 * 
 * ==================================================================================
 */

#include <gtest/gtest.h>
#include <chrono>
#include <random>
#include <cmath>
#include <iostream>

// Library headers
#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"
#include "metodos/lu.hpp"
#include "metodos/cholesky.hpp"
#include "metodos/gauss_jordan.hpp"
#include "metodos/gauss_seidel.hpp"
#include "metodos/conjugate_gradient.hpp"

// As soluções de referência do Octave ficam isoladas em octave_ref.cpp; aqui só
// usamos a interface enxuta, que NÃO inclui cabeçalhos do Octave. Isso evita o
// conflito de nomes entre a `template<class T> class lu;` global do Octave e o
// `namespace lu` da biblioteca (ver octave_ref.hpp).
#include "octave_ref.hpp"

/**
 * ==================================================================================
 * LIBOCTAVE HELPER FUNCTIONS
 * ==================================================================================
 * These functions encapsulate the Octave C++ API usage and handle
 * data conversion between library and Octave formats.
 *
 * A resolução A\b é feita por Matrix::solve(...), que internamente seleciona a
 * fatoração apropriada (LU/QR; Cholesky para SPD) e chama LAPACK diretamente —
 * sem nenhuma chamada a system()/popen()/exec() e sem script Octave.
 */

namespace octave_helpers {

/**
 * Compute residual ||Ax - b||
 */
double compute_residual(const Matriz& A, const Vector& x, const Vector& b) {
    Vector Ax = multiplicar(A, x);
    Vector residual = Ax - b;
    return residual.euclidian_length();
}

// Wrappers finos sobre octave_ref (a conversão e as chamadas ao liboctave vivem
// em octave_ref.cpp). Os nomes são mantidos para preservar os casos de teste.
inline Vector octave_solve_lu(const Matriz& A, const Vector& b) {
    return octave_ref::solve_lu(A, b);
}
inline Vector octave_solve_cholesky(const Matriz& A, const Vector& b) {
    return octave_ref::solve_cholesky(A, b);
}
inline Vector octave_solve_gauss_jordan(const Matriz& A, const Vector& b) {
    return octave_ref::solve_gauss_jordan(A, b);
}
inline Vector octave_solve_gauss_seidel(const Matriz& A, const Vector& b) {
    return octave_ref::solve_gauss_seidel(A, b);
}
inline Vector octave_solve_pcg(const Matriz& A, const Vector& b,
                               double /*tol*/ = 1e-10, int /*max_iter*/ = 1000) {
    return octave_ref::solve_pcg(A, b);
}

} // namespace octave_helpers

/**
 * ==================================================================================
 * TEST CLASS: OctaveComparisonTest
 * ==================================================================================
 */

class OctaveComparisonTest : public ::testing::Test {
protected:
    const double TOLERANCE = 1e-6;
    const int PERFORMANCE_SIZE = 50;
    
    /**
     * Create simple 3x3 test matrix
     */
    Matriz createSimpleMatrix3x3() {
        Matriz A(3, 3);
        A.setValue(0, 0, 4.0);   A.setValue(0, 1, 3.0);   A.setValue(0, 2, 0.0);
        A.setValue(1, 0, 6.0);   A.setValue(1, 1, 3.0);   A.setValue(1, 2, 1.0);
        A.setValue(2, 0, 12.0);  A.setValue(2, 1, 6.0);   A.setValue(2, 2, 1.0);
        return A;
    }
    
    /**
     * Create simple 3x3 SPD matrix
     */
    Matriz createSimpleSPDMatrix3x3() {
        Matriz A(3, 3);
        A.setValue(0, 0, 4.0);   A.setValue(0, 1, 2.0);   A.setValue(0, 2, 1.0);
        A.setValue(1, 0, 2.0);   A.setValue(1, 1, 5.0);   A.setValue(1, 2, 2.0);
        A.setValue(2, 0, 1.0);   A.setValue(2, 1, 2.0);   A.setValue(2, 2, 3.0);
        return A;
    }
    
    /**
     * Create RHS vector
     */
    Vector createSimpleVector3() {
        Vector b(3);
        b.setValue(0, 7.0);
        b.setValue(1, 10.0);
        b.setValue(2, 19.0);
        return b;
    }
    
    /**
     * Create random matrix
     */
    Matriz createRandomMatrix(int n) {
        std::mt19937 rng(42);
        std::uniform_real_distribution<double> dist(0.1, 1.0);
        
        Matriz M(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j <= i; ++j) {
                M.setValue(i, j, dist(rng));
            }
            for (int j = i + 1; j < n; ++j) {
                M.setValue(i, j, 0.0);
            }
        }
        
        Matriz Mt = M.transpose();
        Matriz A = M * Mt;
        for (int i = 0; i < n; ++i) {
            A.setValue(i, i, A.getValue(i, i) + n * 10.0);
        }
        return A;
    }
    
    /**
     * Create random vector
     */
    Vector createRandomVector(int n) {
        std::mt19937 rng(42);
        std::uniform_real_distribution<double> dist(-10.0, 10.0);
        Vector b(n);
        for (int i = 0; i < n; ++i) {
            b.setValue(i, dist(rng));
        }
        return b;
    }
};

/**
 * ==================================================================================
 * COMPARATIVE TESTS: LU DECOMPOSITION
 * ==================================================================================
 */

TEST_F(OctaveComparisonTest, LU_SmallMatrix_3x3_Comparison) {
    Matriz A = createSimpleMatrix3x3();
    Vector b = createSimpleVector3();
    
    // Solve using library
    auto lib_start = std::chrono::high_resolution_clock::now();
    Vector x_lib = lu::solve(A, b);
    auto lib_end = std::chrono::high_resolution_clock::now();
    auto lib_duration = std::chrono::duration_cast<std::chrono::microseconds>(lib_end - lib_start);
    
    // Solve using Octave
    auto oct_start = std::chrono::high_resolution_clock::now();
    Vector x_oct = octave_helpers::octave_solve_lu(A, b);
    auto oct_end = std::chrono::high_resolution_clock::now();
    auto oct_duration = std::chrono::duration_cast<std::chrono::microseconds>(oct_end - oct_start);
    
    // Compare solutions
    for (int i = 0; i < 3; ++i) {
        EXPECT_NEAR(x_lib.getValue(i), x_oct.getValue(i), TOLERANCE);
    }
    
    // Verify both solutions satisfy Ax ≈ b
    double res_lib = octave_helpers::compute_residual(A, x_lib, b);
    double res_oct = octave_helpers::compute_residual(A, x_oct, b);
    EXPECT_LT(res_lib, TOLERANCE);
    EXPECT_LT(res_oct, TOLERANCE);
    
    std::cout << "LU (3x3) - Library: " << lib_duration.count() << " μs, "
              << "Octave: " << oct_duration.count() << " μs\n";
}

TEST_F(OctaveComparisonTest, LU_LargeMatrix_50x50_Comparison) {
    int n = PERFORMANCE_SIZE;
    Matriz A = createRandomMatrix(n);
    Vector b = createRandomVector(n);
    
    // Solve using library
    auto lib_start = std::chrono::high_resolution_clock::now();
    Vector x_lib = lu::solve(A, b);
    auto lib_end = std::chrono::high_resolution_clock::now();
    auto lib_duration = std::chrono::duration_cast<std::chrono::milliseconds>(lib_end - lib_start);
    
    // Solve using Octave
    auto oct_start = std::chrono::high_resolution_clock::now();
    Vector x_oct = octave_helpers::octave_solve_lu(A, b);
    auto oct_end = std::chrono::high_resolution_clock::now();
    auto oct_duration = std::chrono::duration_cast<std::chrono::milliseconds>(oct_end - oct_start);
    
    // Compare solutions (with scaled tolerance for larger systems)
    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(x_lib.getValue(i), x_oct.getValue(i), TOLERANCE * n);
    }
    
    // Verify residuals
    double res_lib = octave_helpers::compute_residual(A, x_lib, b);
    double res_oct = octave_helpers::compute_residual(A, x_oct, b);
    EXPECT_LT(res_lib, TOLERANCE * n);
    EXPECT_LT(res_oct, TOLERANCE * n);
    
    std::cout << "LU (50x50) - Library: " << lib_duration.count() << " ms, "
              << "Octave: " << oct_duration.count() << " ms\n";
}

/**
 * ==================================================================================
 * COMPARATIVE TESTS: CHOLESKY DECOMPOSITION (SPD MATRICES)
 * ==================================================================================
 */

TEST_F(OctaveComparisonTest, Cholesky_SmallMatrix_3x3_Comparison) {
    Matriz A = createSimpleSPDMatrix3x3();
    Vector b = createSimpleVector3();
    
    // Solve using library
    auto lib_start = std::chrono::high_resolution_clock::now();
    Vector x_lib = cholesky::solve(A, b);
    auto lib_end = std::chrono::high_resolution_clock::now();
    auto lib_duration = std::chrono::duration_cast<std::chrono::microseconds>(lib_end - lib_start);
    
    // Solve using Octave (uses Cholesky for SPD matrices)
    auto oct_start = std::chrono::high_resolution_clock::now();
    Vector x_oct = octave_helpers::octave_solve_cholesky(A, b);
    auto oct_end = std::chrono::high_resolution_clock::now();
    auto oct_duration = std::chrono::duration_cast<std::chrono::microseconds>(oct_end - oct_start);
    
    // Compare solutions
    for (int i = 0; i < 3; ++i) {
        EXPECT_NEAR(x_lib.getValue(i), x_oct.getValue(i), TOLERANCE);
    }
    
    // Verify residuals
    double res_lib = octave_helpers::compute_residual(A, x_lib, b);
    double res_oct = octave_helpers::compute_residual(A, x_oct, b);
    EXPECT_LT(res_lib, TOLERANCE);
    EXPECT_LT(res_oct, TOLERANCE);
    
    std::cout << "Cholesky (3x3) - Library: " << lib_duration.count() << " μs, "
              << "Octave: " << oct_duration.count() << " μs\n";
}

TEST_F(OctaveComparisonTest, Cholesky_LargeMatrix_50x50_Comparison) {
    int n = PERFORMANCE_SIZE;
    Matriz A = createRandomMatrix(n);
    Vector b = createRandomVector(n);
    
    // Solve using library
    auto lib_start = std::chrono::high_resolution_clock::now();
    Vector x_lib = cholesky::solve(A, b);
    auto lib_end = std::chrono::high_resolution_clock::now();
    auto lib_duration = std::chrono::duration_cast<std::chrono::milliseconds>(lib_end - lib_start);
    
    // Solve using Octave
    auto oct_start = std::chrono::high_resolution_clock::now();
    Vector x_oct = octave_helpers::octave_solve_cholesky(A, b);
    auto oct_end = std::chrono::high_resolution_clock::now();
    auto oct_duration = std::chrono::duration_cast<std::chrono::milliseconds>(oct_end - oct_start);
    
    // Compare solutions
    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(x_lib.getValue(i), x_oct.getValue(i), TOLERANCE * n);
    }
    
    // Verify residuals
    double res_lib = octave_helpers::compute_residual(A, x_lib, b);
    double res_oct = octave_helpers::compute_residual(A, x_oct, b);
    EXPECT_LT(res_lib, TOLERANCE * n);
    EXPECT_LT(res_oct, TOLERANCE * n);
    
    std::cout << "Cholesky (50x50) - Library: " << lib_duration.count() << " ms, "
              << "Octave: " << oct_duration.count() << " ms\n";
}

/**
 * ==================================================================================
 * COMPARATIVE TESTS: GAUSS-JORDAN ELIMINATION
 * ==================================================================================
 */

TEST_F(OctaveComparisonTest, GaussJordan_SmallMatrix_3x3_Comparison) {
    Matriz A = createSimpleMatrix3x3();
    Vector b = createSimpleVector3();
    
    // Solve using library
    auto lib_start = std::chrono::high_resolution_clock::now();
    Vector x_lib = gaussjordan::solve(A, b);
    auto lib_end = std::chrono::high_resolution_clock::now();
    auto lib_duration = std::chrono::duration_cast<std::chrono::microseconds>(lib_end - lib_start);
    
    // Solve using Octave
    auto oct_start = std::chrono::high_resolution_clock::now();
    Vector x_oct = octave_helpers::octave_solve_gauss_jordan(A, b);
    auto oct_end = std::chrono::high_resolution_clock::now();
    auto oct_duration = std::chrono::duration_cast<std::chrono::microseconds>(oct_end - oct_start);
    
    // Compare solutions
    for (int i = 0; i < 3; ++i) {
        EXPECT_NEAR(x_lib.getValue(i), x_oct.getValue(i), TOLERANCE);
    }
    
    // Verify residuals
    double res_lib = octave_helpers::compute_residual(A, x_lib, b);
    double res_oct = octave_helpers::compute_residual(A, x_oct, b);
    EXPECT_LT(res_lib, TOLERANCE);
    EXPECT_LT(res_oct, TOLERANCE);
    
    std::cout << "Gauss-Jordan (3x3) - Library: " << lib_duration.count() << " μs, "
              << "Octave: " << oct_duration.count() << " μs\n";
}

TEST_F(OctaveComparisonTest, GaussJordan_LargeMatrix_50x50_Comparison) {
    int n = PERFORMANCE_SIZE;
    Matriz A = createRandomMatrix(n);
    Vector b = createRandomVector(n);
    
    // Solve using library
    auto lib_start = std::chrono::high_resolution_clock::now();
    Vector x_lib = gaussjordan::solve(A, b);
    auto lib_end = std::chrono::high_resolution_clock::now();
    auto lib_duration = std::chrono::duration_cast<std::chrono::milliseconds>(lib_end - lib_start);
    
    // Solve using Octave
    auto oct_start = std::chrono::high_resolution_clock::now();
    Vector x_oct = octave_helpers::octave_solve_gauss_jordan(A, b);
    auto oct_end = std::chrono::high_resolution_clock::now();
    auto oct_duration = std::chrono::duration_cast<std::chrono::milliseconds>(oct_end - oct_start);
    
    // Compare solutions
    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(x_lib.getValue(i), x_oct.getValue(i), TOLERANCE * n);
    }
    
    // Verify residuals
    double res_lib = octave_helpers::compute_residual(A, x_lib, b);
    double res_oct = octave_helpers::compute_residual(A, x_oct, b);
    EXPECT_LT(res_lib, TOLERANCE * n);
    EXPECT_LT(res_oct, TOLERANCE * n);
    
    std::cout << "Gauss-Jordan (50x50) - Library: " << lib_duration.count() << " ms, "
              << "Octave: " << oct_duration.count() << " ms\n";
}

/**
 * ==================================================================================
 * COMPARATIVE TESTS: CONJUGATE GRADIENT
 * ==================================================================================
 * Note: Direct CG comparison with Octave's C++ API is limited.
 * This test uses Octave's lssolve as approximation.
 */

TEST_F(OctaveComparisonTest, ConjugateGradient_SmallMatrix_3x3_Comparison) {
    Matriz A = createSimpleSPDMatrix3x3();
    Vector b = createSimpleVector3();
    
    // Solve using library
    auto lib_start = std::chrono::high_resolution_clock::now();
    Vector x_lib = conjugate_gradient::solve(A, b, 1e-10, 1000);
    auto lib_end = std::chrono::high_resolution_clock::now();
    auto lib_duration = std::chrono::duration_cast<std::chrono::microseconds>(lib_end - lib_start);
    
    // Solve using Octave's PCG
    auto oct_start = std::chrono::high_resolution_clock::now();
    Vector x_oct = octave_helpers::octave_solve_pcg(A, b);
    auto oct_end = std::chrono::high_resolution_clock::now();
    auto oct_duration = std::chrono::duration_cast<std::chrono::microseconds>(oct_end - oct_start);
    
    // Compare solutions
    for (int i = 0; i < 3; ++i) {
        EXPECT_NEAR(x_lib.getValue(i), x_oct.getValue(i), TOLERANCE);
    }
    
    // Verify residuals
    double res_lib = octave_helpers::compute_residual(A, x_lib, b);
    double res_oct = octave_helpers::compute_residual(A, x_oct, b);
    EXPECT_LT(res_lib, TOLERANCE);
    EXPECT_LT(res_oct, TOLERANCE);
    
    std::cout << "Conjugate Gradient (3x3) - Library: " << lib_duration.count() << " μs, "
              << "Octave: " << oct_duration.count() << " μs\n";
}

TEST_F(OctaveComparisonTest, ConjugateGradient_LargeMatrix_50x50_Comparison) {
    int n = PERFORMANCE_SIZE;
    Matriz A = createRandomMatrix(n);
    Vector b = createRandomVector(n);
    
    // Solve using library
    auto lib_start = std::chrono::high_resolution_clock::now();
    Vector x_lib = conjugate_gradient::solve(A, b, 1e-10, 1000);
    auto lib_end = std::chrono::high_resolution_clock::now();
    auto lib_duration = std::chrono::duration_cast<std::chrono::milliseconds>(lib_end - lib_start);
    
    // Solve using Octave's PCG
    auto oct_start = std::chrono::high_resolution_clock::now();
    Vector x_oct = octave_helpers::octave_solve_pcg(A, b, 1e-10, 1000);
    auto oct_end = std::chrono::high_resolution_clock::now();
    auto oct_duration = std::chrono::duration_cast<std::chrono::milliseconds>(oct_end - oct_start);
    
    // Compare solutions
    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(x_lib.getValue(i), x_oct.getValue(i), TOLERANCE * n);
    }
    
    // Verify residuals
    double res_lib = octave_helpers::compute_residual(A, x_lib, b);
    double res_oct = octave_helpers::compute_residual(A, x_oct, b);
    EXPECT_LT(res_lib, TOLERANCE * n);
    EXPECT_LT(res_oct, TOLERANCE * n);
    
    std::cout << "Conjugate Gradient (50x50) - Library: " << lib_duration.count() << " ms, "
              << "Octave: " << oct_duration.count() << " ms\n";
}

/**
 * ==================================================================================
 * COMPARATIVE TESTS: GAUSS-SEIDEL
 * ==================================================================================
 * Octave não expõe Gauss-Seidel diretamente na API C++, então a referência
 * usada aqui é a solução direta A\b (Matrix::solve), comparando a solução
 * iterativa da biblioteca com a solução numérica de alta precisão.
 * ==================================================================================
 */

TEST_F(OctaveComparisonTest, GaussSeidel_SmallMatrix_3x3_Comparison) {
    Matriz A = createSimpleSPDMatrix3x3();
    Vector b = createSimpleVector3();

    // Solve using library (iterative)
    auto lib_start = std::chrono::high_resolution_clock::now();
    solver::IterationInfo info;
    Vector x_lib = gaussseidel::solve(A, b, 1e-10, 5000, &info, 1.0);
    auto lib_end = std::chrono::high_resolution_clock::now();
    auto lib_duration = std::chrono::duration_cast<std::chrono::microseconds>(lib_end - lib_start);

    // Solve using Octave direct reference
    auto oct_start = std::chrono::high_resolution_clock::now();
    Vector x_oct = octave_helpers::octave_solve_gauss_seidel(A, b);
    auto oct_end = std::chrono::high_resolution_clock::now();
    auto oct_duration = std::chrono::duration_cast<std::chrono::microseconds>(oct_end - oct_start);

    EXPECT_TRUE(info.converged) << "Gauss-Seidel did not converge on 3x3 SPD matrix";

    // Compare solutions
    for (int i = 0; i < 3; ++i) {
        EXPECT_NEAR(x_lib.getValue(i), x_oct.getValue(i), TOLERANCE);
    }

    // Verify residuals
    double res_lib = octave_helpers::compute_residual(A, x_lib, b);
    double res_oct = octave_helpers::compute_residual(A, x_oct, b);
    EXPECT_LT(res_lib, TOLERANCE);
    EXPECT_LT(res_oct, TOLERANCE);

    std::cout << "Gauss-Seidel (3x3) - Library: " << lib_duration.count() << " μs, "
              << "Octave ref: " << oct_duration.count() << " μs\n";
}

TEST_F(OctaveComparisonTest, GaussSeidel_LargeMatrix_50x50_Comparison) {
    int n = PERFORMANCE_SIZE;
    Matriz A = createRandomMatrix(n); // SPD by construction
    Vector b = createRandomVector(n);

    // Solve using library (iterative)
    auto lib_start = std::chrono::high_resolution_clock::now();
    solver::IterationInfo info;
    Vector x_lib = gaussseidel::solve(A, b, 1e-10, 10000, &info, 1.0);
    auto lib_end = std::chrono::high_resolution_clock::now();
    auto lib_duration = std::chrono::duration_cast<std::chrono::milliseconds>(lib_end - lib_start);

    // Solve using Octave direct reference
    auto oct_start = std::chrono::high_resolution_clock::now();
    Vector x_oct = octave_helpers::octave_solve_gauss_seidel(A, b);
    auto oct_end = std::chrono::high_resolution_clock::now();
    auto oct_duration = std::chrono::duration_cast<std::chrono::milliseconds>(oct_end - oct_start);

    EXPECT_TRUE(info.converged) << "Gauss-Seidel did not converge on 50x50 SPD matrix";

    // Compare solutions
    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(x_lib.getValue(i), x_oct.getValue(i), TOLERANCE * n);
    }

    // Verify residuals
    double res_lib = octave_helpers::compute_residual(A, x_lib, b);
    double res_oct = octave_helpers::compute_residual(A, x_oct, b);
    EXPECT_LT(res_lib, TOLERANCE * n);
    EXPECT_LT(res_oct, TOLERANCE * n);

    std::cout << "Gauss-Seidel (50x50) - Library: " << lib_duration.count() << " ms, "
              << "Octave ref: " << oct_duration.count() << " ms\n";
}
