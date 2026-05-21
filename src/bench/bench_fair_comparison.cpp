#include <iostream>
#include <chrono>
#include <vector>
#include <iomanip>
#include <cmath>
#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"
#include "metodos/conjugate_gradient.hpp"
#include "metodos/gauss_seidel.hpp"
#include "metodos/lu.hpp"
#include "metodos/cholesky.hpp"

/**
 * Fair Comparison Benchmark
 * Compares computational performance of different solvers with same stopping criteria
 */

class BenchmarkTimer {
private:
    std::chrono::high_resolution_clock::time_point start_time;
    
public:
    void start() {
        start_time = std::chrono::high_resolution_clock::now();
    }
    
    double elapsed_ms() {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double, std::milli>(end - start_time).count();
    }
};

void generateSPDMatrix(Matriz& A) {
    int n = A.getRows();
    for (int i = 0; i < n; ++i) {
        A.setValue(i, i, 2.0 * n);
        for (int j = i + 1; j < n; ++j) {
            double val = 0.1;
            A.setValue(i, j, val);
            A.setValue(j, i, val);
        }
    }
}

Vector createRandomRHS(int n) {
    Vector b(n);
    for (int i = 0; i < n; ++i) {
        b.setValue(i, 1.0 + 0.1 * (i % 10));
    }
    return b;
}

double computeRelativeResidual(Matriz& A, Vector& x, Vector& b) {
    Vector Ax = multiplicar(A, x);
    Vector residual = Ax - b;
    double residual_norm = residual.euclidian_length();
    double b_norm = b.euclidian_length();
    return residual_norm / b_norm;
}

struct BenchmarkResult {
    std::string method_name;
    int matrix_size;
    double time_ms;
    double residual;
    int iterations;
    bool converged;
    std::string status;
};

int main() {
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "=== Fair Comparison Benchmark ===" << std::endl;
    std::cout << "All solvers tested with consistent tolerances" << std::endl;
    
    std::vector<int> sizes = {50, 100, 200, 500};
    double tolerance = 1e-8;
    int max_iterations = 1000;
    
    std::vector<BenchmarkResult> results;
    
    for (int n : sizes) {
        std::cout << "\n" << std::string(70, '=') << std::endl;
        std::cout << "Matrix Size: " << n << "x" << n << std::endl;
        std::cout << "Tolerance: " << tolerance << ", Max Iterations: " << max_iterations << std::endl;
        std::cout << std::string(70, '=') << std::endl;
        
        Matriz A(n, n);
        generateSPDMatrix(A);
        
        Vector b = createRandomRHS(n);
        
        // === Direct: LU Decomposition ===
        {
            BenchmarkResult res;
            res.method_name = "LU Decomposition";
            res.matrix_size = n;
            res.iterations = 1;
            res.status = "OK";
            
            try {
                BenchmarkTimer timer;
                timer.start();
                Vector x = lu::solve(A, b);
                res.time_ms = timer.elapsed_ms();
                res.residual = computeRelativeResidual(A, x, b);
                res.converged = (res.residual < tolerance);
            } catch (const std::exception& e) {
                res.status = "FAILED";
                res.time_ms = -1.0;
                res.residual = -1.0;
            }
            
            std::cout << "\nLU Decomposition:" << std::endl;
            std::cout << "  Status: " << res.status << std::endl;
            if (res.status == "OK") {
                std::cout << "  Time (ms): " << res.time_ms << std::endl;
                std::cout << "  Relative Residual: " << res.residual << std::endl;
                std::cout << "  Converged: " << (res.converged ? "YES" : "NO") << std::endl;
            }
            results.push_back(res);
        }
        
        // === Direct: Cholesky Decomposition ===
        {
            BenchmarkResult res;
            res.method_name = "Cholesky";
            res.matrix_size = n;
            res.iterations = 1;
            res.status = "OK";
            
            try {
                BenchmarkTimer timer;
                timer.start();
                Vector x = cholesky::solve(A, b);
                res.time_ms = timer.elapsed_ms();
                res.residual = computeRelativeResidual(A, x, b);
                res.converged = (res.residual < tolerance);
            } catch (const std::exception& e) {
                res.status = "FAILED";
                res.time_ms = -1.0;
                res.residual = -1.0;
            }
            
            std::cout << "\nCholesky Decomposition:" << std::endl;
            std::cout << "  Status: " << res.status << std::endl;
            if (res.status == "OK") {
                std::cout << "  Time (ms): " << res.time_ms << std::endl;
                std::cout << "  Relative Residual: " << res.residual << std::endl;
                std::cout << "  Converged: " << (res.converged ? "YES" : "NO") << std::endl;
            }
            results.push_back(res);
        }
        
        // === Iterative: Conjugate Gradient ===
        {
            BenchmarkResult res;
            res.method_name = "Conjugate Gradient";
            res.matrix_size = n;
            res.status = "OK";
            
            try {
                BenchmarkTimer timer;
                timer.start();
                solver::IterationInfo info;
                Vector x = conjugate_gradient::solve(A, b, tolerance, max_iterations, nullptr, &info);
                res.time_ms = timer.elapsed_ms();
                res.residual = computeRelativeResidual(A, x, b);
                res.iterations = info.iterations;
                res.converged = info.converged;
            } catch (const std::exception& e) {
                res.status = "FAILED";
                res.time_ms = -1.0;
                res.residual = -1.0;
            }
            
            std::cout << "\nConjugate Gradient:" << std::endl;
            std::cout << "  Status: " << res.status << std::endl;
            if (res.status == "OK") {
                std::cout << "  Iterations: " << res.iterations << std::endl;
                std::cout << "  Time (ms): " << res.time_ms << std::endl;
                std::cout << "  Relative Residual: " << res.residual << std::endl;
                std::cout << "  Converged: " << (res.converged ? "YES" : "NO") << std::endl;
            }
            results.push_back(res);
        }
        
        // === Iterative: Gauss-Seidel ===
        {
            BenchmarkResult res;
            res.method_name = "Gauss-Seidel";
            res.matrix_size = n;
            res.status = "OK";
            
            try {
                BenchmarkTimer timer;
                timer.start();
                solver::IterationInfo info;
                Vector x = gaussseidel::solve(A, b, tolerance, max_iterations, &info, 1.0);
                res.time_ms = timer.elapsed_ms();
                res.residual = computeRelativeResidual(A, x, b);
                res.iterations = info.iterations;
                res.converged = info.converged;
            } catch (const std::exception& e) {
                res.status = "FAILED";
                res.time_ms = -1.0;
                res.residual = -1.0;
            }
            
            std::cout << "\nGauss-Seidel:" << std::endl;
            std::cout << "  Status: " << res.status << std::endl;
            if (res.status == "OK") {
                std::cout << "  Iterations: " << res.iterations << std::endl;
                std::cout << "  Time (ms): " << res.time_ms << std::endl;
                std::cout << "  Relative Residual: " << res.residual << std::endl;
                std::cout << "  Converged: " << (res.converged ? "YES" : "NO") << std::endl;
            }
            results.push_back(res);
        }
    }
    
    // === Summary Table ===
    std::cout << "\n" << std::string(100, '=') << std::endl;
    std::cout << "SUMMARY TABLE" << std::endl;
    std::cout << std::string(100, '=') << std::endl;
    std::cout << std::setw(20) << "Method" 
              << std::setw(15) << "Matrix Size"
              << std::setw(15) << "Time (ms)"
              << std::setw(15) << "Iterations"
              << std::setw(20) << "Rel. Residual"
              << std::setw(10) << "Status" << std::endl;
    std::cout << std::string(100, '-') << std::endl;
    
    for (const auto& res : results) {
        std::cout << std::setw(20) << res.method_name
                  << std::setw(15) << res.matrix_size
                  << std::setw(15) << (res.time_ms >= 0 ? std::to_string(res.time_ms) : "N/A")
                  << std::setw(15) << res.iterations
                  << std::setw(20) << (res.residual >= 0 ? std::to_string(res.residual) : "N/A")
                  << std::setw(10) << res.status << std::endl;
    }
    
    std::cout << "\nBenchmark Complete" << std::endl;
    return 0;
}
