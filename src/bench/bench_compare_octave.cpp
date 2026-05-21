#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <cmath>
#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"
#include "metodos/conjugate_gradient.hpp"
#include "metodos/gauss_seidel.hpp"
#include "metodos/lu.hpp"
#include "metodos/cholesky.hpp"

/**
 * Benchmark: Compare C++ implementations with Octave results
 * This program generates systems, solves them with our methods,
 * and compares with expected solutions.
 */

void generateRandomSPDMatrix(Matriz& A, int n) {
    // Generate a random symmetric positive definite matrix
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            double val = (double)rand() / RAND_MAX;
            A.setValue(i, j, val);
            A.setValue(j, i, val);
        }
        // Ensure diagonal dominance for SPD
        A.setValue(i, i, A.getValue(i, i) + n);
    }
}

void generateRandomVector(Vector& v) {
        for (int i = 0; i < v.getLength(); ++i) {
        v.setValue(i, (double)rand() / RAND_MAX);
    }
}

double computeResidualNorm(Matriz& A, Vector& x, Vector& b) {
    Vector Ax = multiplicar(A, x);
    Vector residual = Ax - b;
    return residual.euclidian_length();
}

double computeRelativeError(Vector& x_computed, Vector& x_expected) {
    Vector error = x_computed - x_expected;
    double error_norm = error.euclidian_length();
    double expected_norm = x_expected.euclidian_length();
    
    if (expected_norm < 1e-14) return error_norm;
    return error_norm / expected_norm;
}

int main() {
    std::cout << "=== Linear Systems Solver Benchmark ===\n" << std::endl;
    
    int seed = 42;
    srand(seed);
    
    std::vector<int> matrix_sizes = {10, 50, 100, 200};
    
    for (int n : matrix_sizes) {
        std::cout << "\n--- Matrix Size: " << n << "x" << n << " ---" << std::endl;
        
        Matriz A(n, n);
        generateRandomSPDMatrix(A, n);
        
        Vector b(n);
        generateRandomVector(b);
        
        // Create expected solution (x = A^-1 * b conceptually)
        Vector x_expected(n);
        for (int i = 0; i < n; ++i) {
            x_expected.setValue(i, (double)rand() / RAND_MAX);
        }
        Vector b_expected = multiplicar(A, x_expected);
        
        // === Conjugate Gradient ===
        try {
            auto start = std::chrono::high_resolution_clock::now();
            solver::IterationInfo info_cg;
            Vector x_cg = conjugate_gradient::solve(A, b_expected, 1e-10, n, nullptr, &info_cg);
            auto end = std::chrono::high_resolution_clock::now();
            
            double residual_cg = computeResidualNorm(A, x_cg, b_expected);
            double relative_error = computeRelativeError(x_cg, x_expected);
            double time_ms = std::chrono::duration<double, std::milli>(end - start).count();
            
            std::cout << "\nConjugate Gradient:" << std::endl;
            std::cout << "  Iterations: " << info_cg.iterations << std::endl;
            std::cout << "  Residual Norm: " << residual_cg << std::endl;
            std::cout << "  Relative Error: " << relative_error << std::endl;
            std::cout << "  Time (ms): " << time_ms << std::endl;
            std::cout << "  Converged: " << (info_cg.converged ? "YES" : "NO") << std::endl;
        } catch (const std::exception& e) {
            std::cout << "\nConjugate Gradient: FAILED - " << e.what() << std::endl;
        }
        
        // === LU Decomposition ===
        try {
            auto start = std::chrono::high_resolution_clock::now();
            Vector x_lu = lu::solve(A, b_expected);
            auto end = std::chrono::high_resolution_clock::now();
            
            double residual_lu = computeResidualNorm(A, x_lu, b_expected);
            double relative_error = computeRelativeError(x_lu, x_expected);
            double time_ms = std::chrono::duration<double, std::milli>(end - start).count();
            
            std::cout << "\nLU Decomposition:" << std::endl;
            std::cout << "  Residual Norm: " << residual_lu << std::endl;
            std::cout << "  Relative Error: " << relative_error << std::endl;
            std::cout << "  Time (ms): " << time_ms << std::endl;
        } catch (const std::exception& e) {
            std::cout << "\nLU Decomposition: FAILED - " << e.what() << std::endl;
        }
        
        // === Cholesky ===
        try {
            auto start = std::chrono::high_resolution_clock::now();
            Vector x_cholesky = cholesky::solve(A, b_expected);
            auto end = std::chrono::high_resolution_clock::now();
            
            double residual_cholesky = computeResidualNorm(A, x_cholesky, b_expected);
            double relative_error = computeRelativeError(x_cholesky, x_expected);
            double time_ms = std::chrono::duration<double, std::milli>(end - start).count();
            
            std::cout << "\nCholesky Decomposition:" << std::endl;
            std::cout << "  Residual Norm: " << residual_cholesky << std::endl;
            std::cout << "  Relative Error: " << relative_error << std::endl;
            std::cout << "  Time (ms): " << time_ms << std::endl;
        } catch (const std::exception& e) {
            std::cout << "\nCholesky Decomposition: FAILED - " << e.what() << std::endl;
        }
        
        // === Gauss-Seidel ===
        try {
            auto start = std::chrono::high_resolution_clock::now();
            solver::IterationInfo info_gs;
            Vector x_gs = gaussseidel::solve(A, b_expected, 1e-10, 1000, &info_gs);
            auto end = std::chrono::high_resolution_clock::now();
            
            double residual_gs = computeResidualNorm(A, x_gs, b_expected);
            double relative_error = computeRelativeError(x_gs, x_expected);
            double time_ms = std::chrono::duration<double, std::milli>(end - start).count();
            
            std::cout << "\nGauss-Seidel:" << std::endl;
            std::cout << "  Iterations: " << info_gs.iterations << std::endl;
            std::cout << "  Residual Norm: " << residual_gs << std::endl;
            std::cout << "  Relative Error: " << relative_error << std::endl;
            std::cout << "  Time (ms): " << time_ms << std::endl;
            std::cout << "  Converged: " << (info_gs.converged ? "YES" : "NO") << std::endl;
        } catch (const std::exception& e) {
            std::cout << "\nGauss-Seidel: FAILED - " << e.what() << std::endl;
        }
    }
    
    std::cout << "\n=== Benchmark Complete ===" << std::endl;
    return 0;
}
