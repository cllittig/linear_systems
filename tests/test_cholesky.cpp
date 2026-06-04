#include <gtest/gtest.h>
#include <chrono>
#include <random>
#include <cmath>
#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"
#include "metodos/cholesky.hpp"

/**
 * @class CholeskyDecompositionTest
 * @brief Test suite for Cholesky decomposition method
 * 
 * Cholesky decomposition requires symmetric positive definite (SPD) matrices.
 * This suite validates:
 * - Small SPD matrices (3x3) with known decomposition
 * - Stress tests with larger SPD matrices (50x50)
 * - Forward and backward substitution correctness
 * - Floating-point tolerance (1e-6)
 */

class CholeskyDecompositionTest : public ::testing::Test {
protected:
    const double TOLERANCE = 1e-6;
    const int PERFORMANCE_SIZE = 50;
    
    /**
     * Create a simple 3x3 symmetric positive definite matrix
     * Matrix A = [ 4   2   1 ]
     *            [ 2   5   2 ]
     *            [ 1   2   3 ]
     */
    Matriz createSimpleSPDMatrix3x3() {
        Matriz A(3, 3);
        A.setValue(0, 0, 4.0);   A.setValue(0, 1, 2.0);   A.setValue(0, 2, 1.0);
        A.setValue(1, 0, 2.0);   A.setValue(1, 1, 5.0);   A.setValue(1, 2, 2.0);
        A.setValue(2, 0, 1.0);   A.setValue(2, 1, 2.0);   A.setValue(2, 2, 3.0);
        return A;
    }
    
    /**
     * Create a simple RHS vector for 3x3 case
     */
    Vector createSimpleVector3() {
        Vector b(3);
        b.setValue(0, 7.0);
        b.setValue(1, 9.0);
        b.setValue(2, 6.0);
        return b;
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
     * Create a random symmetric positive definite matrix using L*L^T approach
     * Ensures numerical stability for Cholesky
     */
    Matriz createRandomSPDMatrix(int n) {
        std::mt19937 rng(42);
        std::uniform_real_distribution<double> dist(0.1, 1.0);
        
        // Create lower triangular matrix L with positive diagonal
        Matriz L(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < i; ++j) {
                L.setValue(i, j, dist(rng));
            }
            // Diagonal elements: positive
            L.setValue(i, i, dist(rng) + 0.5);
            for (int j = i + 1; j < n; ++j) {
                L.setValue(i, j, 0.0);
            }
        }
        
        // A = L * L^T ensures SPD property
        Matriz Lt = L.transpose();
        Matriz A = L * Lt;
        
        // Add diagonal dominance for better conditioning
        for (int i = 0; i < n; ++i) {
            A.setValue(i, i, A.getValue(i, i) + n * 5.0);
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
 * @test Cholesky_SmallMatrix_3x3_Decomposition
 * Test Cholesky decomposition on 3x3 SPD matrix
 * Verifies that L*L^T reconstructs the original matrix
 */
TEST_F(CholeskyDecompositionTest, Cholesky_SmallMatrix_3x3_Decomposition) {
    Matriz A = createSimpleSPDMatrix3x3();
    
    // Compute Cholesky decomposition
    Matriz L = cholesky::cholesky(A);
    
    // Reconstruct A from L*L^T
    Matriz Lt = L.transpose();
    Matriz A_reconstructed = L * Lt;
    
    // Verify reconstruction: ||A - L*L^T|| should be small
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double expected = A.getValue(i, j);
            double computed = A_reconstructed.getValue(i, j);
            EXPECT_NEAR(expected, computed, TOLERANCE);
        }
    }
}

/**
 * @test Cholesky_SmallMatrix_3x3_Solve
 * Test full solve operation (A*x = b) using Cholesky
 * Verifies correctness of forward and backward substitution
 */
TEST_F(CholeskyDecompositionTest, Cholesky_SmallMatrix_3x3_Solve) {
    Matriz A = createSimpleSPDMatrix3x3();
    Vector b = createSimpleVector3();
    
    // Solve using Cholesky
    Vector x = cholesky::solve(A, b);
    
    // Verify solution: ||Ax - b|| should be near zero
    double residual = computeResidual(A, x, b);
    EXPECT_LT(residual, TOLERANCE);
}

/**
 * @test Cholesky_SmallMatrix_ForwardSubstitution
 * Test forward substitution (L*y = b)
 */
TEST_F(CholeskyDecompositionTest, Cholesky_SmallMatrix_ForwardSubstitution) {
    Matriz A = createSimpleSPDMatrix3x3();
    Vector b = createSimpleVector3();
    
    Matriz L = cholesky::cholesky(A);
    Vector y = cholesky::forwardSubstituition(L, b);
    
    // Verify L*y = b
    Vector Ly = multiplicar(L, y);
    Vector residual_Ly = Ly - b;
    EXPECT_LT(residual_Ly.euclidian_length(), TOLERANCE);
}

/**
 * @test Cholesky_SmallMatrix_BackwardSubstitution
 * Test backward substitution (L^T*x = y)
 */
TEST_F(CholeskyDecompositionTest, Cholesky_SmallMatrix_BackwardSubstitution) {
    Matriz A = createSimpleSPDMatrix3x3();
    Vector b = createSimpleVector3();
    
    Matriz L = cholesky::cholesky(A);
    Vector y = cholesky::forwardSubstituition(L, b);
    Vector x = cholesky::backwardSubstitution(L.transpose(), y);
    
    // Verify A*x = b (composition of forward and backward)
    double residual = computeResidual(A, x, b);
    EXPECT_LT(residual, TOLERANCE);
}

/**
 * @test Cholesky_LargeMatrix_50x50_Correctness
 * Stress test: Cholesky decomposition with 50x50 SPD matrix
 */
TEST_F(CholeskyDecompositionTest, Cholesky_LargeMatrix_50x50_Correctness) {
    int n = PERFORMANCE_SIZE;
    Matriz A = createRandomSPDMatrix(n);
    Vector b = createRandomVector(n);
    
    // Solve using Cholesky
    Vector x = cholesky::solve(A, b);
    
    // Verify solution
    double residual = computeResidual(A, x, b);
    EXPECT_LT(residual, TOLERANCE * n);
}

/**
 * @test Cholesky_LargeMatrix_50x50_Performance
 * Performance benchmark for Cholesky on 50x50 SPD matrix
 */
TEST_F(CholeskyDecompositionTest, Cholesky_LargeMatrix_50x50_Performance) {
    int n = PERFORMANCE_SIZE;
    Matriz A = createRandomSPDMatrix(n);
    Vector b = createRandomVector(n);
    
    auto start = std::chrono::high_resolution_clock::now();
    Vector x = cholesky::solve(A, b);
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    // Verify solution
    double residual = computeResidual(A, x, b);
    EXPECT_LT(residual, TOLERANCE * n);
    
    // Cholesky should be faster than LU (O(n^3/3) vs O(2n^3/3))
    EXPECT_LT(duration.count(), 1000) << "Cholesky (50x50) took " << duration.count() << " ms";
    
    std::cout << "Cholesky decomposition (50x50) took " << duration.count() << " ms\n";
}

/**
 * @test Cholesky_MultipleRHS
 * Test Cholesky solve with multiple right-hand sides (matrix RHS)
 */
TEST_F(CholeskyDecompositionTest, Cholesky_MultipleRHS) {
    Matriz A = createSimpleSPDMatrix3x3();
    
    Matriz B(3, 2);
    B.setValue(0, 0, 7.0);   B.setValue(0, 1, 10.0);
    B.setValue(1, 0, 9.0);   B.setValue(1, 1, 12.0);
    B.setValue(2, 0, 6.0);   B.setValue(2, 1, 8.0);
    
    // Solve for multiple RHS
    Matriz X = cholesky::solve(A, B);
    
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
 * @test Cholesky_SymmetryVerification
 * Verify that the algorithm correctly handles symmetric input
 */
TEST_F(CholeskyDecompositionTest, Cholesky_SymmetryVerification) {
    Matriz A = createSimpleSPDMatrix3x3();
    
    // Verify A is symmetric
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(A.getValue(i, j), A.getValue(j, i), TOLERANCE);
        }
    }
    
    Matriz L = cholesky::cholesky(A);
    
    // L should be lower triangular
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            EXPECT_NEAR(L.getValue(i, j), 0.0, TOLERANCE);
        }
    }
}

/**
 * @test Cholesky_NonSPDMatrix_ThrowsException
 * Verify that non-SPD matrix throws exception
 */
TEST_F(CholeskyDecompositionTest, Cholesky_NonSPDMatrix_ThrowsException) {
    Matriz notSPD(2, 2);
    // Create a non-SPD symmetric matrix (negative eigenvalue)
    notSPD.setValue(0, 0, 1.0);   notSPD.setValue(0, 1, 2.0);
    notSPD.setValue(1, 0, 2.0);   notSPD.setValue(1, 1, 1.0);
    
    Vector b(2);
    b.setValue(0, 1.0);
    b.setValue(1, 1.0);
    
    // Should throw because matrix is not positive definite
    EXPECT_THROW(cholesky::solve(notSPD, b), std::runtime_error);
}
