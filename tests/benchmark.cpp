#include <chrono>
#include <cstdio>
#include <random>
#include <vector>

#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"
#include "metodos/cholesky.hpp"
#include "metodos/conjugate_gradient.hpp"
#include "metodos/gauss_seidel.hpp"
#include "metodos/lu.hpp"

#ifdef USE_LAPACK
#include <lapacke.h>
#endif

using Clock = std::chrono::high_resolution_clock;
using Ms    = std::chrono::duration<double, std::milli>;

#ifndef VARIANTE
#define VARIANTE "padrão"
#endif

static Matriz make_spd(int n, std::mt19937 &rng) {
    std::uniform_real_distribution<double> dist(0.1, 1.0);
    Matriz A(n, n);
    std::vector<double> row_sums(n, 0.0);
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double v = dist(rng);
            A.setValue(i, j, v);
            A.setValue(j, i, v);
            row_sums[i] += v;
            row_sums[j] += v;
        }
    }
    for (int i = 0; i < n; i++)
        A.setValue(i, i, row_sums[i] + 1.0);
    return A;
}

static Vector make_b(int n, std::mt19937 &rng) {
    std::uniform_real_distribution<double> dist(-10.0, 10.0);
    Vector b(n);
    for (int i = 0; i < n; i++)
        b.setValue(i, dist(rng));
    return b;
}

static double residual(const Matriz &A, const Vector &b, const Vector &x) {
    return (multiplicar(A, x) - b).euclidian_length();
}

struct Solver {
    const char *name;
    Vector (*fn)(const Matriz &, const Vector &);
};

static Vector run_lu(const Matriz &A, const Vector &b) {
    return lu::solve(A, b);
}
static Vector run_cholesky(const Matriz &A, const Vector &b) {
    return cholesky::solve(A, b);
}
static Vector run_cg(const Matriz &A, const Vector &b) {
    return conjugate_gradient::solve(A, b);
}
static Vector run_gs(const Matriz &A, const Vector &b) {
    return gaussseidel::solve(A, b);
}

#ifdef USE_LAPACK
// Converte Matriz row-major → vetor column-major (layout que o LAPACK espera)
static std::vector<double> to_colmajor(const Matriz &A) {
    int n = A.getRows(), m = A.getColumns();
    std::vector<double> buf(n * m);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            buf[j * n + i] = A.getValue(i, j);
    return buf;
}

static Vector run_lapack_lu(const Matriz &A, const Vector &b) {
    int n = A.getRows();
    std::vector<double> a = to_colmajor(A);
    std::vector<double> rhs(n);
    for (int i = 0; i < n; i++) rhs[i] = b.getValue(i);
    std::vector<lapack_int> ipiv(n);
    LAPACKE_dgesv(LAPACK_COL_MAJOR, n, 1, a.data(), n, ipiv.data(), rhs.data(), n);
    Vector x(n);
    for (int i = 0; i < n; i++) x.setValue(i, rhs[i]);
    return x;
}

static Vector run_lapack_cholesky(const Matriz &A, const Vector &b) {
    int n = A.getRows();
    std::vector<double> a = to_colmajor(A);
    std::vector<double> rhs(n);
    for (int i = 0; i < n; i++) rhs[i] = b.getValue(i);
    LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', n, a.data(), n);
    LAPACKE_dpotrs(LAPACK_COL_MAJOR, 'L', n, 1, a.data(), n, rhs.data(), n);
    Vector x(n);
    for (int i = 0; i < n; i++) x.setValue(i, rhs[i]);
    return x;
}
#endif

int main() {
    std::mt19937 rng(42);
    const std::vector<int> sizes = {64, 128, 256, 512, 1000, 2000, 5000};
    const Solver solvers[] = {
        {"lu",              run_lu},
        // {"cholesky",        run_cholesky},
        // {"cg",              run_cg},
        // {"gauss_seidel",    run_gs},
#ifdef USE_LAPACK
        {"lapack_lu",       run_lapack_lu},
        // {"lapack_cholesky", run_lapack_cholesky},
#endif
    };

    printf("n,metodo,variante,tempo_ms,residuo\n");

    for (int n : sizes) {
        Matriz A = make_spd(n, rng);
        Vector b = make_b(n, rng);

        for (const auto &s : solvers) {
            auto t0 = Clock::now();
            Vector x = s.fn(A, b);
            double ms = Ms(Clock::now() - t0).count();
            printf("%d,%s,%s,%.3f,%.6e\n", n, s.name, VARIANTE, ms, residual(A, b, x));
        }
    }

    return 0;
}
