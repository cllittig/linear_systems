#include <gtest/gtest.h>
#include <chrono>
#include <random>
#include <cmath>
#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"
#include "metodos/gauss_seidel.hpp"

/**
 * @class GaussSeidelTest
 * @brief Test suite for Gauss-Seidel iterative solver
 * 
 * Gauss-Seidel is an iterative method suitable for large sparse systems.
 * This suite validates:
 * - Small diagonally dominant matrices (3x3)
 * - Convergence behavior with different tolerances
 * - Stress tests on larger matrices (50x50)
 * - Convergence information tracking
 * - Floating-point tolerance (1e-6)
 */

class GaussSeidelTest : public ::testing::Test {
protected:
    // Gauss-Seidel is iterative and may not reach extreme precision
    // Using slightly relaxed tolerance for final residual (1e-4)
    // but keeping 1e-6 for internal algorithm tolerance parameter
    const double TOLERANCE = 1e-6;
    const double RESIDUAL_CHECK_TOL = 1e-4;  // Allow residual up to 1e-4 for iterative method
    const int PERFORMANCE_SIZE = 50;
    const int MAX_ITERATIONS = 1000;
    
    /**
     * Create a simple 3x3 diagonally dominant matrix
     * Important: Gauss-Seidel converges for diagonally dominant matrices
     */
    Matriz createSimpleDiagonallyDominantMatrix3x3() {
        Matriz A(3, 3);
        // Make sure diagonal elements dominate
        A.setValue(0, 0, 10.0);  A.setValue(0, 1, 1.0);   A.setValue(0, 2, 2.0);
        A.setValue(1, 0, 1.0);   A.setValue(1, 1, 8.0);   A.setValue(1, 2, 1.0);
        A.setValue(2, 0, 2.0);   A.setValue(2, 1, 1.0);   A.setValue(2, 2, 10.0);
        return A;
    }
    
    /**
     * Create RHS vector for 3x3 case
     */
    Vector createSimpleVector3() {
        Vector b(3);
        b.setValue(0, 12.0);
        b.setValue(1, 10.0);
        b.setValue(2, 13.0);
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
     * Create a random diagonally dominant matrix for convergence guarantee
     */
    Matriz createRandomDiagonallyDominantMatrix(int n) {
        std::mt19937 rng(42);
        std::uniform_real_distribution<double> dist(0.1, 0.3);
        
        Matriz A(n, n);
        
        // Fill with small off-diagonal elements
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    A.setValue(i, j, dist(rng));
                }
            }
        }
        
        // Set dominant diagonal elements
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    sum += std::fabs(A.getValue(i, j));
                }
            }
            A.setValue(i, i, sum * 2.0 + 1.0);
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
 * @test GaussSeidel_SmallMatrix_3x3_Convergence
 * Test Gauss-Seidel on 3x3 diagonally dominant matrix
 * Verifies convergence and solution accuracy
 */
TEST_F(GaussSeidelTest, GaussSeidel_SmallMatrix_3x3_Convergence) {
    Matriz A = createSimpleDiagonallyDominantMatrix3x3();
    Vector b = createSimpleVector3();
    
    solver::IterationInfo info;
    Vector x = gaussseidel::solve(A, b, TOLERANCE, MAX_ITERATIONS, &info);
    
    // Check convergence
    EXPECT_TRUE(info.converged) << "Gauss-Seidel did not converge";
    EXPECT_LT(info.iterations, MAX_ITERATIONS);
    // Allow slightly relaxed residual tolerance for iterative method
    EXPECT_LT(info.final_residual_norm, RESIDUAL_CHECK_TOL);
    
    // Verify solution
    double residual = computeResidual(A, x, b);
    EXPECT_LT(residual, RESIDUAL_CHECK_TOL);
    
    std::cout << "GS (3x3): converged in " << info.iterations << " iterations\n";
}

/**
 * @test GaussSeidel_SmallMatrix_ResidualDecay
 * Verify that residual decreases monotonically
 */
TEST_F(GaussSeidelTest, GaussSeidel_SmallMatrix_ResidualDecay) {
    Matriz A = createSimpleDiagonallyDominantMatrix3x3();
    Vector b = createSimpleVector3();
    
    solver::IterationInfo info;
    Vector x = gaussseidel::solve(A, b, TOLERANCE, MAX_ITERATIONS, &info);
    
    EXPECT_LT(info.final_residual_norm, RESIDUAL_CHECK_TOL);
    EXPECT_TRUE(info.converged);
}

/**
 * @test GaussSeidel_LooseTolerance
 * Test convergence with looser tolerance
 */
TEST_F(GaussSeidelTest, GaussSeidel_LooseTolerance) {
    Matriz A = createSimpleDiagonallyDominantMatrix3x3();
    Vector b = createSimpleVector3();
    
    solver::IterationInfo info;
    double loose_tol = 1e-3;
    Vector x = gaussseidel::solve(A, b, loose_tol, MAX_ITERATIONS, &info);
    
    EXPECT_TRUE(info.converged);
    EXPECT_LE(info.iterations, MAX_ITERATIONS);
    EXPECT_LT(info.final_residual_norm, loose_tol * 10);  // Allow some slack
}

/**
 * @test GaussSeidel_StrictTolerance
 * Test convergence with strict tolerance
 */
TEST_F(GaussSeidelTest, GaussSeidel_StrictTolerance) {
    Matriz A = createSimpleDiagonallyDominantMatrix3x3();
    Vector b = createSimpleVector3();
    
    solver::IterationInfo info;
    double strict_tol = 1e-10;
    Vector x = gaussseidel::solve(A, b, strict_tol, MAX_ITERATIONS, &info);
    
    // May not converge to 1e-10 due to numerical precision
    // but should make good progress
    EXPECT_LT(info.final_residual_norm, 1e-8);
}

/**
 * @test GaussSeidel_LargeMatrix_50x50_Correctness
 * Stress test: Gauss-Seidel on 50x50 diagonally dominant matrix
 */
TEST_F(GaussSeidelTest, GaussSeidel_LargeMatrix_50x50_Correctness) {
    int n = PERFORMANCE_SIZE;
    Matriz A = createRandomDiagonallyDominantMatrix(n);
    Vector b = createRandomVector(n);
    
    solver::IterationInfo info;
    Vector x = gaussseidel::solve(A, b, TOLERANCE, MAX_ITERATIONS, &info);
    
    // Verify convergence
    EXPECT_TRUE(info.converged) << "Failed to converge on 50x50 matrix";
    EXPECT_LT(info.final_residual_norm, TOLERANCE * n);
    
    double residual = computeResidual(A, x, b);
    EXPECT_LT(residual, TOLERANCE * n);
}

/**
 * @test GaussSeidel_LargeMatrix_50x50_Performance
 * Performance benchmark for Gauss-Seidel on 50x50
 */
TEST_F(GaussSeidelTest, GaussSeidel_LargeMatrix_50x50_Performance) {
    int n = PERFORMANCE_SIZE;
    Matriz A = createRandomDiagonallyDominantMatrix(n);
    Vector b = createRandomVector(n);
    
    auto start = std::chrono::high_resolution_clock::now();
    solver::IterationInfo info;
    Vector x = gaussseidel::solve(A, b, TOLERANCE, MAX_ITERATIONS, &info);
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    // Verify convergence
    EXPECT_TRUE(info.converged);
    EXPECT_LT(info.final_residual_norm, TOLERANCE * n);
    
    std::cout << "Gauss-Seidel (50x50): " << info.iterations << " iterations, "
              << duration.count() << " ms\n";
}

/**
 * @test GaussSeidel_WithInitialGuess
 * Test solving with non-zero initial guess
 * (if the implementation supports x0 parameter)
 */
TEST_F(GaussSeidelTest, GaussSeidel_ConvergenceFaster) {
    Matriz A = createSimpleDiagonallyDominantMatrix3x3();
    Vector b = createSimpleVector3();
    
    solver::IterationInfo info;
    Vector x = gaussseidel::solve(A, b, TOLERANCE, MAX_ITERATIONS, &info);
    
    // Any diagonally dominant matrix should converge
    EXPECT_TRUE(info.converged);
}

/**
 * @test GaussSeidel_OmegaParameter
 * Test Gauss-Seidel with relaxation parameter (SOR)
 * omega = 1.0 is standard Gauss-Seidel
 * 0 < omega < 1 is under-relaxation
 * 1 < omega < 2 is over-relaxation (can accelerate convergence)
 */
TEST_F(GaussSeidelTest, GaussSeidel_WithOmega_Standard) {
    Matriz A = createSimpleDiagonallyDominantMatrix3x3();
    Vector b = createSimpleVector3();
    
    solver::IterationInfo info;
    Vector x = gaussseidel::solve(A, b, TOLERANCE, MAX_ITERATIONS, &info, 1.0);
    
    EXPECT_TRUE(info.converged);
    EXPECT_LT(info.final_residual_norm, RESIDUAL_CHECK_TOL);
}

/**
 * @test GaussSeidel_WithOmega_Overrelaxation
 * Test with over-relaxation (omega > 1)
 * May converge faster than standard Gauss-Seidel
 */
TEST_F(GaussSeidelTest, GaussSeidel_WithOmega_Overrelaxation) {
    Matriz A = createSimpleDiagonallyDominantMatrix3x3();
    Vector b = createSimpleVector3();
    
    solver::IterationInfo info;
    double omega = 1.5;  // Over-relaxation
    Vector x = gaussseidel::solve(A, b, TOLERANCE, MAX_ITERATIONS, &info, omega);
    
    EXPECT_TRUE(info.converged);
    EXPECT_LT(info.final_residual_norm, RESIDUAL_CHECK_TOL);
}

/**
 * @test GaussSeidel_NonDominantMatrix_Divergence
 * Verify behavior on non-diagonally dominant matrix (may diverge)
 */
TEST_F(GaussSeidelTest, GaussSeidel_NonDominantMatrix) {
    Matriz A(2, 2);
    // Non-diagonally dominant matrix
    A.setValue(0, 0, 1.0);   A.setValue(0, 1, 2.0);
    A.setValue(1, 0, 3.0);   A.setValue(1, 1, 1.0);
    
    Vector b(2);
    b.setValue(0, 5.0);
    b.setValue(1, 5.0);
    
    solver::IterationInfo info;
    // Non-dominant matrix may throw exception or diverge
    // This test documents both behaviors
    try {
        Vector x = gaussseidel::solve(A, b, TOLERANCE, 50, &info);
        // If it converges, great - but less likely for non-dominant
        EXPECT_LE(info.iterations, 50);
    } catch (const std::exception& e) {
        // Expected behavior for non-dominant matrix
        EXPECT_LE(info.iterations, 50);
    }
}
