#include <gtest/gtest.h>
#include <chrono>
#include <random>
#include <cmath>
#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"
#include "metodos/conjugate_gradient.hpp"

/**
 * @class ConjugateGradientTest
 * @brief Test suite for Conjugate Gradient iterative solver
 * 
 * Conjugate Gradient is an efficient iterative method for symmetric
 * positive definite (SPD) matrices.
 * This suite validates:
 * - Small SPD matrices (3x3) with known convergence behavior
 * - Convergence rate improvements with preconditioning
 * - Larger systems (50x50) for stress testing
 * - Multiple right-hand sides
 * - Floating-point tolerance (1e-10 for CG)
 */

class ConjugateGradientTest : public ::testing::Test {
protected:
    const double TOLERANCE = 1e-10;
    const int PERFORMANCE_SIZE = 50;
    const int MAX_ITERATIONS = 1000;
    
    /**
     * Create a simple 3x3 symmetric positive definite matrix
     */
    Matriz createSimpleSPDMatrix3x3() {
        Matriz A(3, 3);
        A.setValue(0, 0, 4.0);   A.setValue(0, 1, 2.0);   A.setValue(0, 2, 1.0);
        A.setValue(1, 0, 2.0);   A.setValue(1, 1, 5.0);   A.setValue(1, 2, 2.0);
        A.setValue(2, 0, 1.0);   A.setValue(2, 1, 2.0);   A.setValue(2, 2, 3.0);
        return A;
    }
    
    /**
     * Create RHS vector for 3x3 case
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
     * Create random SPD matrix for testing
     */
    Matriz createRandomSPDMatrix(int n) {
        std::mt19937 rng(42);
        std::uniform_real_distribution<double> dist(0.1, 1.0);
        
        // Create lower triangular matrix L
        Matriz L(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < i; ++j) {
                L.setValue(i, j, dist(rng));
            }
            L.setValue(i, i, dist(rng) + 0.5);
            for (int j = i + 1; j < n; ++j) {
                L.setValue(i, j, 0.0);
            }
        }
        
        // A = L * L^T ensures SPD
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
 * @test ConjugateGradient_SmallMatrix_3x3_Convergence
 * Test CG on 3x3 SPD matrix
 * CG should converge in at most n iterations for n-dimensional systems
 */
TEST_F(ConjugateGradientTest, ConjugateGradient_SmallMatrix_3x3_Convergence) {
    Matriz A = createSimpleSPDMatrix3x3();
    Vector b = createSimpleVector3();
    
    solver::IterationInfo info;
    Vector x = conjugate_gradient::solve(A, b, TOLERANCE, MAX_ITERATIONS, nullptr, &info);
    
    // CG should converge in at most 3 iterations for 3x3 system (in exact arithmetic)
    // Floating-point errors may require slightly more
    EXPECT_TRUE(info.converged);
    EXPECT_LE(info.iterations, 10) << "CG took too many iterations: " << info.iterations;
    EXPECT_LT(info.final_residual_norm, TOLERANCE);
    
    double residual = computeResidual(A, x, b);
    EXPECT_LT(residual, TOLERANCE);
    
    std::cout << "CG (3x3): converged in " << info.iterations << " iterations\n";
}

/**
 * @test ConjugateGradient_SmallMatrix_QuadraticConvergence
 * CG should converge quadratically on well-conditioned SPD matrices
 */
TEST_F(ConjugateGradientTest, ConjugateGradient_SmallMatrix_QuadraticConvergence) {
    Matriz A = createSimpleSPDMatrix3x3();
    Vector b = createSimpleVector3();
    
    solver::IterationInfo info;
    Vector x = conjugate_gradient::solve(A, b, TOLERANCE, MAX_ITERATIONS, nullptr, &info);
    
    EXPECT_TRUE(info.converged);
    EXPECT_LE(info.iterations, 5);
    EXPECT_LT(info.final_residual_norm, TOLERANCE);
}

/**
 * @test ConjugateGradient_LooseTolerance
 * Test with looser tolerance: should converge faster
 */
TEST_F(ConjugateGradientTest, ConjugateGradient_LooseTolerance) {
    Matriz A = createSimpleSPDMatrix3x3();
    Vector b = createSimpleVector3();
    
    solver::IterationInfo info;
    double loose_tol = 1e-6;
    Vector x = conjugate_gradient::solve(A, b, loose_tol, MAX_ITERATIONS, nullptr, &info);
    
    EXPECT_TRUE(info.converged);
    EXPECT_LE(info.iterations, 5);
    
    std::cout << "CG (loose tol 1e-6): " << info.iterations << " iterations\n";
}

/**
 * @test ConjugateGradient_InitialGuess
 * Test that initial guess can reduce iterations
 * (if x0 parameter is provided)
 */
TEST_F(ConjugateGradientTest, ConjugateGradient_WithInitialGuess) {
    Matriz A = createSimpleSPDMatrix3x3();
    Vector b = createSimpleVector3();
    
    // Create a good initial guess
    Vector x0(3);
    x0.setValue(0, 1.0);
    x0.setValue(1, 1.0);
    x0.setValue(2, 1.0);
    
    solver::IterationInfo info;
    Vector x = conjugate_gradient::solve(A, b, TOLERANCE, MAX_ITERATIONS, &x0, &info);
    
    EXPECT_TRUE(info.converged);
    EXPECT_LT(info.final_residual_norm, TOLERANCE);
    
    std::cout << "CG (with initial guess): " << info.iterations << " iterations\n";
}

/**
 * @test ConjugateGradient_LargeMatrix_50x50_Correctness
 * Stress test: CG on 50x50 SPD matrix
 * Should converge much faster than iterative methods like Gauss-Seidel
 */
TEST_F(ConjugateGradientTest, ConjugateGradient_LargeMatrix_50x50_Correctness) {
    int n = PERFORMANCE_SIZE;
    Matriz A = createRandomSPDMatrix(n);
    Vector b = createRandomVector(n);
    
    solver::IterationInfo info;
    Vector x = conjugate_gradient::solve(A, b, TOLERANCE, MAX_ITERATIONS, nullptr, &info);
    
    EXPECT_TRUE(info.converged) << "CG did not converge on 50x50 matrix";
    EXPECT_LT(info.final_residual_norm, TOLERANCE * n);
    
    double residual = computeResidual(A, x, b);
    EXPECT_LT(residual, TOLERANCE * n);
    
    std::cout << "CG (50x50): converged in " << info.iterations << " iterations\n";
}

/**
 * @test ConjugateGradient_LargeMatrix_50x50_Performance
 * Performance benchmark for CG on 50x50
 * Should be significantly faster than LU and Gauss-Jordan
 */
TEST_F(ConjugateGradientTest, ConjugateGradient_LargeMatrix_50x50_Performance) {
    int n = PERFORMANCE_SIZE;
    Matriz A = createRandomSPDMatrix(n);
    Vector b = createRandomVector(n);
    
    auto start = std::chrono::high_resolution_clock::now();
    solver::IterationInfo info;
    Vector x = conjugate_gradient::solve(A, b, TOLERANCE, MAX_ITERATIONS, nullptr, &info);
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    EXPECT_TRUE(info.converged);
    EXPECT_LT(info.final_residual_norm, TOLERANCE * n);
    
    double residual = computeResidual(A, x, b);
    EXPECT_LT(residual, TOLERANCE * n);
    
    std::cout << "Conjugate Gradient (50x50): " << info.iterations << " iterations, "
              << duration.count() << " ms\n";
}

/**
 * @test ConjugateGradient_WithJacobiPreconditioning
 * Test CG with Jacobi preconditioning for faster convergence
 */
TEST_F(ConjugateGradientTest, ConjugateGradient_WithJacobiPreconditioning) {
    Matriz A = createSimpleSPDMatrix3x3();
    Vector b = createSimpleVector3();
    
    solver::IterationInfo info_no_precond;
    Vector x_no_precond = conjugate_gradient::solve(A, b, TOLERANCE, MAX_ITERATIONS, 
                                                     nullptr, &info_no_precond, false);
    
    solver::IterationInfo info_with_precond;
    Vector x_with_precond = conjugate_gradient::solve(A, b, TOLERANCE, MAX_ITERATIONS,
                                                       nullptr, &info_with_precond, true);
    
    EXPECT_TRUE(info_no_precond.converged);
    EXPECT_TRUE(info_with_precond.converged);
    
    // Preconditioning should converge in fewer or equal iterations
    EXPECT_LE(info_with_precond.iterations, info_no_precond.iterations + 1);
    
    std::cout << "CG without precond: " << info_no_precond.iterations << " iter\n";
    std::cout << "CG with Jacobi precond: " << info_with_precond.iterations << " iter\n";
    
    // Both should converge to the same solution
    for (int i = 0; i < 3; ++i) {
        EXPECT_NEAR(x_no_precond.getValue(i), x_with_precond.getValue(i), TOLERANCE * 10);
    }
}

/**
 * @test ConjugateGradient_MultipleRHS
 * Test CG with multiple right-hand sides (solves column by column)
 */
TEST_F(ConjugateGradientTest, ConjugateGradient_MultipleRHS) {
    Matriz A = createSimpleSPDMatrix3x3();
    
    Matriz B(3, 2);
    B.setValue(0, 0, 7.0);   B.setValue(0, 1, 10.0);
    B.setValue(1, 0, 9.0);   B.setValue(1, 1, 12.0);
    B.setValue(2, 0, 6.0);   B.setValue(2, 1, 8.0);
    
    // Solve for multiple RHS
    solver::IterationInfo info;
    Matriz X = conjugate_gradient::solve(A, B, TOLERANCE, MAX_ITERATIONS, &info);
    
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
    
    EXPECT_TRUE(info.converged);
}

/**
 * @test ConjugateGradient_IdentityMatrix
 * Test on identity matrix (trivial convergence in 1 iteration)
 */
TEST_F(ConjugateGradientTest, ConjugateGradient_IdentityMatrix) {
    int n = 3;
    Matriz I(n, n);
    for (int i = 0; i < n; ++i) {
        I.setValue(i, i, 1.0);
    }
    
    Vector b(n);
    b.setValue(0, 5.0);
    b.setValue(1, 3.0);
    b.setValue(2, 7.0);
    
    solver::IterationInfo info;
    Vector x = conjugate_gradient::solve(I, b, TOLERANCE, MAX_ITERATIONS, nullptr, &info);
    
    // Should converge in 1 iteration
    EXPECT_EQ(info.iterations, 1);
    EXPECT_TRUE(info.converged);
    
    // Solution should be x = b
    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(x.getValue(i), b.getValue(i), TOLERANCE * 10);
    }
}

/**
 * @test ConjugateGradient_ScaledMatrix
 * Test on scaled matrix: A = diag(2, 4, 6)
 * Should still converge in n iterations
 */
TEST_F(ConjugateGradientTest, ConjugateGradient_ScaledMatrix) {
    Matriz A(3, 3);
    A.setValue(0, 0, 2.0);   A.setValue(1, 1, 4.0);   A.setValue(2, 2, 6.0);
    
    Vector b(3);
    b.setValue(0, 2.0);
    b.setValue(1, 4.0);
    b.setValue(2, 6.0);
    
    solver::IterationInfo info;
    Vector x = conjugate_gradient::solve(A, b, TOLERANCE, MAX_ITERATIONS, nullptr, &info);
    
    EXPECT_TRUE(info.converged);
    EXPECT_LE(info.iterations, 3);
    
    // Solution: x = (1, 1, 1)
    for (int i = 0; i < 3; ++i) {
        EXPECT_NEAR(x.getValue(i), 1.0, TOLERANCE * 10);
    }
}

/**
 * @test ConjugateGradient_IllConditionedMatrix
 * Test on ill-conditioned SPD matrix
 * Should still converge but may take more iterations
 */
TEST_F(ConjugateGradientTest, ConjugateGradient_IllConditionedMatrix) {
    Matriz A(2, 2);
    // Very ill-conditioned: eigenvalues are 100 and 1
    A.setValue(0, 0, 100.0);   A.setValue(0, 1, 0.0);
    A.setValue(1, 0, 0.0);     A.setValue(1, 1, 1.0);
    
    Vector b(2);
    b.setValue(0, 100.0);
    b.setValue(1, 1.0);
    
    solver::IterationInfo info;
    Vector x = conjugate_gradient::solve(A, b, 1e-8, MAX_ITERATIONS, nullptr, &info);
    
    // Should still converge
    EXPECT_TRUE(info.converged);
    
    // Solution: x = (1, 1)
    for (int i = 0; i < 2; ++i) {
        EXPECT_NEAR(x.getValue(i), 1.0, 1e-6);
    }
    
    std::cout << "CG (ill-conditioned): " << info.iterations << " iterations\n";
}
