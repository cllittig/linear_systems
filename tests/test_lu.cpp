#include <gtest/gtest.h>
#include <chrono>
#include <random>
#include <cmath>
#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"
#include "metodos/lu.hpp"

/**
 * @class LUDecompositionTest
 * @brief Test suite for LU decomposition method
 * 
 * This suite validates the LU decomposition implementation with:
 * - Small matrices (3x3) with known values computed by hand
 * - Stress tests with larger matrices (50x50) for performance validation
 * - Floating-point tolerance using EXPECT_NEAR (1e-6)
 */

class LUDecompositionTest : public ::testing::Test {
protected:
    const double TOLERANCE = 1e-6;
    const int PERFORMANCE_SIZE = 50;
    
    // Helper function to create a simple symmetric positive definite matrix
    Matriz createSimpleMatrix3x3() {
        Matriz A(3, 3);
        // Known 3x3 matrix for manual verification
        // [  4  3  0 ]
        // [  6  3  1 ]
        // [ 12  6  1 ]
        A.setValue(0, 0, 4.0);   A.setValue(0, 1, 3.0);   A.setValue(0, 2, 0.0);
        A.setValue(1, 0, 6.0);   A.setValue(1, 1, 3.0);   A.setValue(1, 2, 1.0);
        A.setValue(2, 0, 12.0);  A.setValue(2, 1, 6.0);   A.setValue(2, 2, 1.0);
        return A;
    }
    
    // Helper function to create a simple RHS vector for 3x3 case
    Vector createSimpleVector3() {
        Vector b(3);
        // [ 7  ]
        // [ 10 ]
        // [ 19 ]
        b.setValue(0, 7.0);
        b.setValue(1, 10.0);
        b.setValue(2, 19.0);
        return b;
    }
    
    // Expected solution for 3x3 case (computed manually or with external tool)
    // Solution: x = [1, 1, 1]
    Vector getExpectedSolution3() {
        Vector x(3);
        x.setValue(0, 1.0);
        x.setValue(1, 1.0);
        x.setValue(2, 1.0);
        return x;
    }
    
    // Helper to verify solution by computing residual: ||Ax - b||
    double computeResidual(const Matriz& A, const Vector& x, const Vector& b) {
        Vector Ax = multiplicar(A, x);
        Vector residual = Ax - b;
        return residual.euclidian_length();
    }
    
    // Helper to create a large random positive definite matrix
    Matriz createRandomMatrix(int n) {
        std::mt19937 rng(42);  // Fixed seed for reproducibility
        std::uniform_real_distribution<double> dist(0.1, 1.0);
        
        Matriz M(n, n);
        // Create random lower triangular matrix L
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j <= i; ++j) {
                M.setValue(i, j, dist(rng));
            }
            for (int j = i + 1; j < n; ++j) {
                M.setValue(i, j, 0.0);
            }
        }
        
        // Create A = L * L^T to ensure it's positive definite
        // A will be symmetric and positive definite
        Matriz Mt = M.transpose();
        Matriz A = M * Mt;
        
        // Add diagonal dominance for better conditioning
        for (int i = 0; i < n; ++i) {
            A.setValue(i, i, A.getValue(i, i) + n * 10.0);
        }
        
        return A;
    }
    
    // Create random RHS vector
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
 * @test LU_SmallMatrix_3x3
 * Test LU decomposition with a small 3x3 matrix with known solution
 * Verifies that the computed solution satisfies Ax = b
 */
TEST_F(LUDecompositionTest, LU_SmallMatrix_3x3) {
    Matriz A = createSimpleMatrix3x3();
    Vector b = createSimpleVector3();
    Vector expectedX = getExpectedSolution3();
    
    // Solve using LU decomposition
    Vector x = lu::solve(A, b);
    
    // Verify solution: check residual is near zero
    double residual = computeResidual(A, x, b);
    EXPECT_NEAR(residual, 0.0, TOLERANCE);
    
    // Also check each component of the solution
    for (int i = 0; i < 3; ++i) {
        EXPECT_NEAR(x.getValue(i), expectedX.getValue(i), TOLERANCE);
    }
}

/**
 * @test LU_SmallMatrix_Residual
 * Test that the computed solution minimizes the residual ||Ax - b||
 */
TEST_F(LUDecompositionTest, LU_SmallMatrix_Residual) {
    Matriz A = createSimpleMatrix3x3();
    Vector b = createSimpleVector3();
    
    Vector x = lu::solve(A, b);
    double residual = computeResidual(A, x, b);
    
    EXPECT_LT(residual, TOLERANCE);
}

/**
 * @test LU_LargeMatrix_50x50_Correctness
 * Stress test: LU decomposition with 50x50 matrix
 * Validates numerical stability and convergence
 */
TEST_F(LUDecompositionTest, LU_LargeMatrix_50x50_Correctness) {
    int n = PERFORMANCE_SIZE;
    Matriz A = createRandomMatrix(n);
    Vector b = createRandomVector(n);
    
    // Solve using LU decomposition
    Vector x = lu::solve(A, b);
    
    // Verify solution satisfies Ax ≈ b
    double residual = computeResidual(A, x, b);
    EXPECT_LT(residual, TOLERANCE * n);  // Scale tolerance for larger problems
}

/**
 * @test LU_LargeMatrix_50x50_Performance
 * Performance benchmark for LU decomposition on 50x50 matrix
 * Measures execution time and ensures it completes within reasonable time
 */
TEST_F(LUDecompositionTest, LU_LargeMatrix_50x50_Performance) {
    int n = PERFORMANCE_SIZE;
    Matriz A = createRandomMatrix(n);
    Vector b = createRandomVector(n);
    
    auto start = std::chrono::high_resolution_clock::now();
    Vector x = lu::solve(A, b);
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    // Verify solution is correct
    double residual = computeResidual(A, x, b);
    EXPECT_LT(residual, TOLERANCE * n);
    
    // Performance benchmark: should complete in reasonable time (< 1 second for 50x50)
    EXPECT_LT(duration.count(), 1000) << "LU decomposition took " << duration.count() << " ms";
    
    std::cout << "LU decomposition (50x50) took " << duration.count() << " ms\n";
}

/**
 * @test LU_MultipleRHS
 * Test LU decomposition with multiple right-hand sides (matrix RHS)
 */
TEST_F(LUDecompositionTest, LU_MultipleRHS_2x2) {
    Matriz A(2, 2);
    A.setValue(0, 0, 2.0);  A.setValue(0, 1, 1.0);
    A.setValue(1, 0, 1.0);  A.setValue(1, 1, 2.0);
    
    Matriz B(2, 2);
    B.setValue(0, 0, 5.0);  B.setValue(0, 1, 6.0);
    B.setValue(1, 0, 5.0);  B.setValue(1, 1, 8.0);
    
    // Solve using LU
    Matriz X = lu::lu(A, B);
    
    // Verify each column: A * x_i ≈ b_i
    for (int col = 0; col < 2; ++col) {
        Vector x(2);
        x.setValue(0, X.getValue(0, col));
        x.setValue(1, X.getValue(1, col));
        
        Vector b(2);
        b.setValue(0, B.getValue(0, col));
        b.setValue(1, B.getValue(1, col));
        
        double residual = computeResidual(A, x, b);
        EXPECT_NEAR(residual, 0.0, TOLERANCE);
    }
}

/**
 * @test LU_SingularMatrix_ThrowsException
 * Verify that LU decomposition throws exception for singular matrix
 */
TEST_F(LUDecompositionTest, LU_SingularMatrix_ThrowsException) {
    Matriz singular(2, 2);
    // Create a singular matrix (linearly dependent rows)
    singular.setValue(0, 0, 1.0);  singular.setValue(0, 1, 2.0);
    singular.setValue(1, 0, 2.0);  singular.setValue(1, 1, 4.0);
    
    Vector b(2);
    b.setValue(0, 3.0);
    b.setValue(1, 6.0);
    
    EXPECT_THROW(lu::solve(singular, b), std::runtime_error);
}
