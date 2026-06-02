#include <gtest/gtest.h>
#include <chrono>
#include <random>
#include <cmath>
#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"
#include "metodos/gauss_jordan.hpp"

/**
 * @class GaussJordanTest
 * @brief Test suite for Gauss-Jordan elimination
 * 
 * Gauss-Jordan is a direct elimination method (like LU) but eliminates
 * above and below the pivot to obtain the reduced row echelon form.
 * This suite validates:
 * - Small matrices (3x3) with known solutions
 * - Numerical stability with floating-point arithmetic
 * - Larger systems (50x50) for stress testing
 * - Floating-point tolerance (1e-6)
 */

class GaussJordanTest : public ::testing::Test {
protected:
    const double TOLERANCE = 1e-6;
    const int PERFORMANCE_SIZE = 50;
    
    /**
     * Create a simple 3x3 test matrix with known solution
     */
    Matriz createSimpleMatrix3x3() {
        Matriz A(3, 3);
        // System: 4*x + 3*y + 0*z = 7
        //         6*x + 3*y + 1*z = 10
        //        12*x + 6*y + 1*z = 19
        A.setValue(0, 0, 4.0);   A.setValue(0, 1, 3.0);   A.setValue(0, 2, 0.0);
        A.setValue(1, 0, 6.0);   A.setValue(1, 1, 3.0);   A.setValue(1, 2, 1.0);
        A.setValue(2, 0, 12.0);  A.setValue(2, 1, 6.0);   A.setValue(2, 2, 1.0);
        return A;
    }
    
    /**
     * Create RHS vector for 3x3 case
     */
    Vector createSimpleVector3() {
        Vector b(3);
        b.setValue(0, 7.0);
        b.setValue(1, 10.0);
        b.setValue(2, 19.0);
        return b;
    }
    
    /**
     * Expected solution
     */
    Vector getExpectedSolution3() {
        Vector x(3);
        x.setValue(0, 1.0);
        x.setValue(1, 1.0);
        x.setValue(2, 1.0);
        return x;
    }
    
    /**
     * Helper to compute residual ||Ax - b||
     */
    double computeResidual(const Matriz& A, const Vector& x, const Vector& b) {
        Vector Ax = multiplicar(A, x);
        Vector residual = Ax - b;
        return residual.euclidian_length();
    }
    
    /**
     * Create random matrix with good conditioning
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
        
        // Create A = M * M^T for better conditioning
        Matriz Mt = M.transpose();
        Matriz A = M * Mt;
        
        // Add diagonal dominance
        for (int i = 0; i < n; ++i) {
            A.setValue(i, i, A.getValue(i, i) + n * 10.0);
        }
        
        return A;
    }
    
    /**
     * Create random RHS vector
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
 * @test GaussJordan_SmallMatrix_3x3
 * Test Gauss-Jordan on 3x3 system with known solution
 */
TEST_F(GaussJordanTest, GaussJordan_SmallMatrix_3x3) {
    Matriz A = createSimpleMatrix3x3();
    Vector b = createSimpleVector3();
    Vector expectedX = getExpectedSolution3();
    
    // Solve using Gauss-Jordan
    Vector x = gaussjordan::solve(A, b);
    
    // Verify solution
    double residual = computeResidual(A, x, b);
    EXPECT_NEAR(residual, 0.0, TOLERANCE);
    
    // Check individual components
    for (int i = 0; i < 3; ++i) {
        EXPECT_NEAR(x.getValue(i), expectedX.getValue(i), TOLERANCE * 10);
    }
}

/**
 * @test GaussJordan_SmallMatrix_Residual
 * Verify residual minimization
 */
TEST_F(GaussJordanTest, GaussJordan_SmallMatrix_Residual) {
    Matriz A = createSimpleMatrix3x3();
    Vector b = createSimpleVector3();
    
    Vector x = gaussjordan::solve(A, b);
    double residual = computeResidual(A, x, b);
    
    EXPECT_LT(residual, TOLERANCE);
}

/**
 * @test GaussJordan_TwoByTwo
 * Simple 2x2 system
 */
TEST_F(GaussJordanTest, GaussJordan_TwoByTwo) {
    Matriz A(2, 2);
    A.setValue(0, 0, 2.0);   A.setValue(0, 1, 1.0);
    A.setValue(1, 0, 1.0);   A.setValue(1, 1, 2.0);
    
    Vector b(2);
    b.setValue(0, 5.0);
    b.setValue(1, 5.0);
    
    Vector x = gaussjordan::solve(A, b);
    
    // Solution should be approximately [5/3, 5/3]
    double expected = 5.0 / 3.0;
    EXPECT_NEAR(x.getValue(0), expected, TOLERANCE);
    EXPECT_NEAR(x.getValue(1), expected, TOLERANCE);
    
    double residual = computeResidual(A, x, b);
    EXPECT_LT(residual, TOLERANCE);
}

/**
 * @test GaussJordan_FourByFour
 * 4x4 system
 */
TEST_F(GaussJordanTest, GaussJordan_FourByFour) {
    Matriz A(4, 4);
    // Simple diagonal system scaled by 2
    A.setValue(0, 0, 2.0);   A.setValue(1, 1, 2.0);
    A.setValue(2, 2, 2.0);   A.setValue(3, 3, 2.0);
    // Off-diagonal: small values
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (i != j) {
                A.setValue(i, j, 0.1);
            }
        }
    }
    
    Vector b(4);
    for (int i = 0; i < 4; ++i) {
        b.setValue(i, 2.4);  // Each component contributes 0.1*4 + 2.0*x = 2.4
    }
    
    Vector x = gaussjordan::solve(A, b);
    
    // Verify solution
    double residual = computeResidual(A, x, b);
    EXPECT_LT(residual, TOLERANCE);
}

/**
 * @test GaussJordan_LargeMatrix_50x50_Correctness
 * Stress test: Gauss-Jordan on 50x50 matrix
 */
TEST_F(GaussJordanTest, GaussJordan_LargeMatrix_50x50_Correctness) {
    int n = PERFORMANCE_SIZE;
    Matriz A = createRandomMatrix(n);
    Vector b = createRandomVector(n);
    
    Vector x = gaussjordan::solve(A, b);
    
    // Verify solution
    double residual = computeResidual(A, x, b);
    EXPECT_LT(residual, TOLERANCE * n);
}

/**
 * @test GaussJordan_LargeMatrix_50x50_Performance
 * Performance benchmark for Gauss-Jordan on 50x50
 */
TEST_F(GaussJordanTest, GaussJordan_LargeMatrix_50x50_Performance) {
    int n = PERFORMANCE_SIZE;
    Matriz A = createRandomMatrix(n);
    Vector b = createRandomVector(n);
    
    auto start = std::chrono::high_resolution_clock::now();
    Vector x = gaussjordan::solve(A, b);
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    // Verify solution
    double residual = computeResidual(A, x, b);
    EXPECT_LT(residual, TOLERANCE * n);
    
    // Gauss-Jordan is typically slower than LU due to operations above diagonal
    EXPECT_LT(duration.count(), 2000) << "Gauss-Jordan (50x50) took " << duration.count() << " ms";
    
    std::cout << "Gauss-Jordan elimination (50x50) took " << duration.count() << " ms\n";
}

/**
 * @test GaussJordan_MultipleRHS
 * Test with multiple right-hand sides
 */
TEST_F(GaussJordanTest, GaussJordan_MultipleRHS) {
    Matriz A(3, 3);
    A.setValue(0, 0, 2.0);   A.setValue(0, 1, 1.0);   A.setValue(0, 2, 0.0);
    A.setValue(1, 0, 1.0);   A.setValue(1, 1, 3.0);   A.setValue(1, 2, 1.0);
    A.setValue(2, 0, 0.0);   A.setValue(2, 1, 1.0);   A.setValue(2, 2, 2.0);
    
    Matriz B(3, 2);
    B.setValue(0, 0, 5.0);   B.setValue(0, 1, 6.0);
    B.setValue(1, 0, 8.0);   B.setValue(1, 1, 9.0);
    B.setValue(2, 0, 4.0);   B.setValue(2, 1, 5.0);
    
    // Solve for multiple RHS
    Matriz X = gaussjordan::gauss_jordan(A, B);
    
    // Verify each column
    for (int col = 0; col < 2; ++col) {
        Vector x(3);
        x.setValue(0, X.getValue(0, col));
        x.setValue(1, X.getValue(1, col));
        x.setValue(2, X.getValue(2, col));
        
        Vector b(3);
        b.setValue(0, B.getValue(0, col));
        b.setValue(1, B.getValue(1, col));
        b.setValue(2, B.getValue(2, col));
        
        double residual = computeResidual(A, x, b);
        EXPECT_LT(residual, TOLERANCE);
    }
}

/**
 * @test GaussJordan_SingularMatrix_ThrowsException
 * Verify exception handling for singular matrix
 */
TEST_F(GaussJordanTest, GaussJordan_SingularMatrix_ThrowsException) {
    Matriz singular(2, 2);
    singular.setValue(0, 0, 1.0);   singular.setValue(0, 1, 2.0);
    singular.setValue(1, 0, 2.0);   singular.setValue(1, 1, 4.0);
    
    Vector b(2);
    b.setValue(0, 3.0);
    b.setValue(1, 6.0);
    
    EXPECT_THROW(gaussjordan::solve(singular, b), std::runtime_error);
}

/**
 * @test GaussJordan_AlmostSingularMatrix
 * Test behavior with ill-conditioned matrix
 */
TEST_F(GaussJordanTest, GaussJordan_AlmostSingularMatrix) {
    Matriz A(2, 2);
    // Almost singular: very small determinant
    double eps = 1e-10;
    A.setValue(0, 0, 1.0);       A.setValue(0, 1, 1.0 - eps);
    A.setValue(1, 0, 1.0 + eps); A.setValue(1, 1, 1.0);
    
    Vector b(2);
    b.setValue(0, 2.0);
    b.setValue(1, 2.0);
    
    // May throw or give poor results
    try {
        Vector x = gaussjordan::solve(A, b);
        // If it doesn't throw, residual may be large
        double residual = computeResidual(A, x, b);
        EXPECT_LT(residual, 1e-3);  // Loose tolerance for ill-conditioned
    } catch (const std::runtime_error& e) {
        // Expected for nearly singular matrix
        EXPECT_TRUE(true);
    }
}

/**
 * @test GaussJordan_IdentityMatrix
 * Test with identity matrix (trivial but important)
 */
TEST_F(GaussJordanTest, GaussJordan_IdentityMatrix) {
    int n = 3;
    Matriz I(n, n);
    for (int i = 0; i < n; ++i) {
        I.setValue(i, i, 1.0);
    }
    
    Vector b(n);
    b.setValue(0, 5.0);
    b.setValue(1, 3.0);
    b.setValue(2, 7.0);
    
    Vector x = gaussjordan::solve(I, b);
    
    // For I*x = b, solution should be x = b
    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(x.getValue(i), b.getValue(i), TOLERANCE);
    }
}
