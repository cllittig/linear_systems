#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"
#include "metodos/conjugate_gradient.hpp"
#include "metodos/cholesky.hpp"
#include "metodos/lu.hpp"
#include "metodos/gauss_seidel.hpp"
#include "metodos/gauss_jordan.hpp"

using namespace std;

static void write_matrix_to_file(const Matriz &A, const string &path) {
    ofstream f(path);
    for (int i = 0; i < A.getRows(); ++i) {
        for (int j = 0; j < A.getColumns(); ++j) {
            f << A.getValue(i, j);
            if (j + 1 < A.getColumns()) f << " ";
        }
        f << "\n";
    }
}

static void write_vector_to_file(const Vector &v, const string &path) {
    ofstream f(path);
    for (int i = 0; i < v.getLength(); ++i) {
        f << v.getValue(i);
        if (i + 1 < v.getLength()) f << "\n";
    }
    f << "\n";
}

static Vector read_vector_from_file(const string &path) {
    ifstream f(path);
    vector<double> data;
    double x;
    while (f >> x) data.push_back(x);
    Vector v((int)data.size());
    for (int i = 0; i < (int)data.size(); ++i) v.setValue(i, data[i]);
    return v;
}

static Matriz read_matrix_from_file(const string &path) {
    ifstream f(path);
    vector<vector<double>> rows;
    string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        stringstream ss(line);
        vector<double> row;
        double val;
        while (ss >> val) row.push_back(val);
        if (!row.empty()) rows.push_back(row);
    }
    if (rows.empty()) return Matriz();
    int n = rows.size();
    int m = rows[0].size();
    Matriz M(n, m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) M.setValue(i, j, rows[i][j]);
    }
    return M;
}

static void generate_octave_script_fair() {
    ofstream script("octave_fair_benchmark.m");
    script << "% Fair benchmark: single session, tic/toc timing\n";
    script << "format long e;\n\n";

    // Test 1: CG 3x3
    script << "% Test 1: CG 3x3\n";
    script << "A1 = [4 1 0; 1 3 1; 0 1 2];\n";
    script << "b1 = [1; 2; 3];\n";
    script << "tic; x1 = A1 \\ b1; t1 = toc; % Using backslash for simple case\n";
    script << "dlmwrite('octave_sol_cg3.txt', x1, 'delimiter', ' ');\n";
    script << "dlmwrite('octave_time_cg3.txt', t1, 'delimiter', ' ');\n\n";

    // Test 2: CG 5x5
    script << "% Test 2: CG 5x5 (tridiagonal)\n";
    script << "A2 = diag(4*ones(5,1)) + diag(ones(4,1),1) + diag(ones(4,1),-1);\n";
    script << "b2 = [1; 2; 3; 4; 5];\n";
    script << "tic; x2 = A2 \\ b2; t2 = toc;\n";
    script << "dlmwrite('octave_sol_cg5.txt', x2, 'delimiter', ' ');\n";
    script << "dlmwrite('octave_time_cg5.txt', t2, 'delimiter', ' ');\n\n";

    // Test 3: Cholesky 3x3 (symmetric positive definite)
    script << "% Test 3: Cholesky 3x3\n";
    script << "A3 = [4 1 0; 1 3 1; 0 1 2];\n";
    script << "b3 = [1; 2; 3];\n";
    script << "tic;\n";
    script << "  R = chol(A3);\n";
    script << "  y = R' \\ b3;\n";
    script << "  x3 = R \\ y;\n";
    script << "t3 = toc;\n";
    script << "dlmwrite('octave_sol_chol3.txt', x3, 'delimiter', ' ');\n";
    script << "dlmwrite('octave_time_chol3.txt', t3, 'delimiter', ' ');\n\n";

    // Test 4: LU 3x3
    script << "% Test 4: LU 3x3\n";
    script << "A4 = [2 1 -1; -3 4 2; 1 -1 3];\n";
    script << "b4 = [1; 7; 3];\n";
    script << "tic;\n";
    script << "  [L,U,P] = lu(A4);\n";
    script << "  x4 = U \\ (L \\ (P*b4));\n";
    script << "t4 = toc;\n";
    script << "dlmwrite('octave_sol_lu3.txt', x4, 'delimiter', ' ');\n";
    script << "dlmwrite('octave_time_lu3.txt', t4, 'delimiter', ' ');\n\n";

    // Test 5: Gauss-Jordan 3x3
    script << "% Test 5: Gauss-Jordan 3x3 (via rref)\n";
    script << "A5 = [3 2 -1; 2 -2 4; -1 0.5 -1];\n";
    script << "b5 = [1; -2; 0];\n";
    script << "tic;\n";
    script << "  M = [A5 b5];\n";
    script << "  R = rref(M);\n";
    script << "  x5 = R(:, end);\n";
    script << "t5 = toc;\n";
    script << "dlmwrite('octave_sol_gj3.txt', x5, 'delimiter', ' ');\n";
    script << "dlmwrite('octave_time_gj3.txt', t5, 'delimiter', ' ');\n\n";

    // Test 6: Gauss-Seidel 5x5
    script << "% Test 6: Gauss-Seidel 5x5\n";
    script << "A6 = diag(4*ones(5,1)) + diag(ones(4,1),1) + diag(ones(4,1),-1);\n";
    script << "b6 = [1; 2; 3; 4; 5];\n";
    script << "function x = gs_iter(A,b,tol,maxIter)\n";
    script << "  n = length(b); x = zeros(n,1);\n";
    script << "  for iter=1:maxIter\n";
    script << "    x_old = x;\n";
    script << "    for i=1:n\n";
    script << "      s = A(i,1:i-1)*x(1:i-1) + A(i,i+1:n)*x_old(i+1:n);\n";
    script << "      x(i) = (b(i) - s)/A(i,i);\n";
    script << "    end\n";
    script << "    if norm(x - x_old, Inf) < tol, break; end\n";
    script << "  end\n";
    script << "endfunction\n";
    script << "tic; x6 = gs_iter(A6, b6, 1e-8, 10000); t6 = toc;\n";
    script << "dlmwrite('octave_sol_gs5.txt', x6, 'delimiter', ' ');\n";
    script << "dlmwrite('octave_time_gs5.txt', t6, 'delimiter', ' ');\n\n";

    // Test 7: Multi-RHS 4x2
    script << "% Test 7: Multi-RHS 4x2\n";
    script << "A7 = diag(6*ones(4,1)) + diag(ones(3,1),1) + diag(ones(3,1),-1);\n";
    script << "B7 = [1 2; 2 2.5; 3 3; 4 4];\n";
    script << "tic; X7 = A7 \\ B7; t7 = toc;\n";
    script << "dlmwrite('octave_sol_multi.txt', X7, 'delimiter', ' ');\n";
    script << "dlmwrite('octave_time_multi.txt', t7, 'delimiter', ' ');\n\n";

    // Test 8: Hilbert 6x6
    script << "% Test 8: Hilbert 6x6 (ill-conditioned)\n";
    script << "n8 = 6;\n";
    script << "H = zeros(n8, n8);\n";
    script << "for i=1:n8\n";
    script << "  for j=1:n8\n";
    script << "    H(i,j) = 1/(i+j-1);\n";
    script << "  end\n";
    script << "end\n";
    script << "b8 = ones(n8, 1);\n";
    script << "tic; x8 = H \\ b8; t8 = toc;\n";
    script << "dlmwrite('octave_sol_hilbert.txt', x8, 'delimiter', ' ');\n";
    script << "dlmwrite('octave_time_hilbert.txt', t8, 'delimiter', ' ');\n\n";

    script << "exit;\n";
    script.close();
}

int main() {
    // Generate Octave script that tests within a single session
    generate_octave_script_fair();

    // Execute Octave once
    cout << "Running single Octave session for all tests...\n";
    auto oct_start = chrono::high_resolution_clock::now();
    system("octave --quiet octave_fair_benchmark.m");
    auto oct_end = chrono::high_resolution_clock::now();
    double octave_total_ms = chrono::duration<double, milli>(oct_end - oct_start).count();

    // Read Octave's own timings from files
    auto read_octave_time = [](const string &path) -> double {
        ifstream f(path);
        double t = 0;
        f >> t;
        return t * 1000; // convert to milliseconds
    };

    struct Row { string test; string method; double diff; double t_oct; double t_lib; };
    vector<Row> rows;

    // Now run C++ library tests
    const int runs = 1; // single run to match Octave

    // 1) CG 3x3
    {
        Matriz A(3,3);
        A.setValue(0,0,4); A.setValue(0,1,1); A.setValue(0,2,0);
        A.setValue(1,0,1); A.setValue(1,1,3); A.setValue(1,2,1);
        A.setValue(2,0,0); A.setValue(2,1,1); A.setValue(2,2,2);
        Vector b(3); b.setValue(0,1); b.setValue(1,2); b.setValue(2,3);

        Vector x_oct = read_vector_from_file("octave_sol_cg3.txt");
        double t_oct = read_octave_time("octave_time_cg3.txt");

        auto l0 = chrono::high_resolution_clock::now();
        Vector x_lib = conjugate_gradient::solve(A, b, 1e-12, 1000);
        auto l1 = chrono::high_resolution_clock::now();
        double t_lib = chrono::duration<double, milli>(l1 - l0).count();

        double diff = 0;
        for (int i = 0; i < 3; ++i) diff = max(diff, fabs(x_oct.getValue(i) - x_lib.getValue(i)));
        rows.push_back({"CG_3x3", "conjugate_gradient", diff, t_oct, t_lib});
    }

    // 2) CG 5x5
    {
        const int n = 5;
        Matriz A(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                double val = (i == j) ? 4.0 : ((abs(i - j) == 1) ? 1.0 : 0.0);
                A.setValue(i, j, val);
            }
        }
        Vector b(n);
        for (int i = 0; i < n; ++i) b.setValue(i, 1.0 + i);

        Vector x_oct = read_vector_from_file("octave_sol_cg5.txt");
        double t_oct = read_octave_time("octave_time_cg5.txt");

        auto l0 = chrono::high_resolution_clock::now();
        Vector x_lib = conjugate_gradient::solve(A, b, 1e-12, 1000);
        auto l1 = chrono::high_resolution_clock::now();
        double t_lib = chrono::duration<double, milli>(l1 - l0).count();

        double diff = 0;
        for (int i = 0; i < n; ++i) diff = max(diff, fabs(x_oct.getValue(i) - x_lib.getValue(i)));
        rows.push_back({"CG_5x5", "conjugate_gradient", diff, t_oct, t_lib});
    }

    // 3) Cholesky 3x3
    {
        Matriz A(3, 3);
        A.setValue(0, 0, 4); A.setValue(0, 1, 1); A.setValue(0, 2, 0);
        A.setValue(1, 0, 1); A.setValue(1, 1, 3); A.setValue(1, 2, 1);
        A.setValue(2, 0, 0); A.setValue(2, 1, 1); A.setValue(2, 2, 2);
        Vector b(3);
        b.setValue(0, 1); b.setValue(1, 2); b.setValue(2, 3);

        Vector x_oct = read_vector_from_file("octave_sol_chol3.txt");
        double t_oct = read_octave_time("octave_time_chol3.txt");

        auto l0 = chrono::high_resolution_clock::now();
        Vector x_lib = cholesky::solve(A, b);
        auto l1 = chrono::high_resolution_clock::now();
        double t_lib = chrono::duration<double, milli>(l1 - l0).count();

        double diff = 0;
        for (int i = 0; i < 3; ++i) diff = max(diff, fabs(x_oct.getValue(i) - x_lib.getValue(i)));
        rows.push_back({"Cholesky_3x3", "cholesky", diff, t_oct, t_lib});
    }

    // 4) LU 3x3
    {
        Matriz A(3, 3);
        A.setValue(0, 0, 2); A.setValue(0, 1, 1); A.setValue(0, 2, -1);
        A.setValue(1, 0, -3); A.setValue(1, 1, 4); A.setValue(1, 2, 2);
        A.setValue(2, 0, 1); A.setValue(2, 1, -1); A.setValue(2, 2, 3);
        Vector b(3);
        b.setValue(0, 1); b.setValue(1, 7); b.setValue(2, 3);

        Vector x_oct = read_vector_from_file("octave_sol_lu3.txt");
        double t_oct = read_octave_time("octave_time_lu3.txt");

        auto l0 = chrono::high_resolution_clock::now();
        Vector x_lib = lu::solve(A, b);
        auto l1 = chrono::high_resolution_clock::now();
        double t_lib = chrono::duration<double, milli>(l1 - l0).count();

        double diff = 0;
        for (int i = 0; i < 3; ++i) diff = max(diff, fabs(x_oct.getValue(i) - x_lib.getValue(i)));
        rows.push_back({"LU_3x3", "lu", diff, t_oct, t_lib});
    }

    // 5) Gauss-Jordan 3x3
    {
        Matriz A(3, 3);
        A.setValue(0, 0, 3); A.setValue(0, 1, 2); A.setValue(0, 2, -1);
        A.setValue(1, 0, 2); A.setValue(1, 1, -2); A.setValue(1, 2, 4);
        A.setValue(2, 0, -1); A.setValue(2, 1, 0.5); A.setValue(2, 2, -1);
        Vector b(3);
        b.setValue(0, 1); b.setValue(1, -2); b.setValue(2, 0);

        Vector x_oct = read_vector_from_file("octave_sol_gj3.txt");
        double t_oct = read_octave_time("octave_time_gj3.txt");

        auto l0 = chrono::high_resolution_clock::now();
        Vector x_lib = gaussjordan::solve(A, b);
        auto l1 = chrono::high_resolution_clock::now();
        double t_lib = chrono::duration<double, milli>(l1 - l0).count();

        double diff = 0;
        for (int i = 0; i < 3; ++i) diff = max(diff, fabs(x_oct.getValue(i) - x_lib.getValue(i)));
        rows.push_back({"GaussJordan_3x3", "gaussjordan", diff, t_oct, t_lib});
    }

    // 6) Gauss-Seidel 5x5
    {
        const int n = 5;
        Matriz A(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                double val = (i == j) ? 4.0 : ((abs(i - j) == 1) ? 1.0 : 0.0);
                A.setValue(i, j, val);
            }
        }
        Vector b(n);
        for (int i = 0; i < n; ++i) b.setValue(i, 1.0 + i);

        Vector x_oct = read_vector_from_file("octave_sol_gs5.txt");
        double t_oct = read_octave_time("octave_time_gs5.txt");

        auto l0 = chrono::high_resolution_clock::now();
        Vector x_lib = gaussseidel::solve(A, b);
        auto l1 = chrono::high_resolution_clock::now();
        double t_lib = chrono::duration<double, milli>(l1 - l0).count();

        double diff = 0;
        for (int i = 0; i < n; ++i) diff = max(diff, fabs(x_oct.getValue(i) - x_lib.getValue(i)));
        rows.push_back({"GaussSeidel_5x5", "gaussseidel", diff, t_oct, t_lib});
    }

    // 7) Multi-RHS 4x2
    {
        const int n = 4;
        Matriz A(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                double val = (i == j) ? 6.0 : ((abs(i - j) == 1) ? 1.0 : 0.0);
                A.setValue(i, j, val);
            }
        }
        Matriz B(n, 2);
        for (int i = 0; i < n; ++i) {
            B.setValue(i, 0, 1.0 + i);
            B.setValue(i, 1, 2.0 + 0.5 * i);
        }

        Matriz X_oct = read_matrix_from_file("octave_sol_multi.txt");
        double t_oct = read_octave_time("octave_time_multi.txt");

        auto l0 = chrono::high_resolution_clock::now();
        Matriz X_lib = cholesky::solve(A, B);
        auto l1 = chrono::high_resolution_clock::now();
        double t_lib = chrono::duration<double, milli>(l1 - l0).count();

        double diff = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < 2; ++j) {
                diff = max(diff, fabs(X_oct.getValue(i, j) - X_lib.getValue(i, j)));
            }
        }
        rows.push_back({"MultiRHS_4x2", "cholesky", diff, t_oct, t_lib});
    }

    // 8) Hilbert 6x6
    {
        int n = 6;
        Matriz H(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                H.setValue(i, j, 1.0 / (i + j + 1));
            }
        }
        Vector b(n);
        for (int i = 0; i < n; ++i) b.setValue(i, 1.0);

        Vector x_oct = read_vector_from_file("octave_sol_hilbert.txt");
        double t_oct = read_octave_time("octave_time_hilbert.txt");

        auto l0 = chrono::high_resolution_clock::now();
        solver::IterationInfo info;
        Vector x_lib = conjugate_gradient::solve(H, b, 1e-10, 5000, nullptr, &info, true);
        auto l1 = chrono::high_resolution_clock::now();
        double t_lib = chrono::duration<double, milli>(l1 - l0).count();

        double diff = 0;
        for (int i = 0; i < n; ++i) diff = max(diff, fabs(x_oct.getValue(i) - x_lib.getValue(i)));
        rows.push_back({"Hilbert_6x6", "conjugate_gradient_precond", diff, t_oct, t_lib});
    }

    // Print results
    cout << "\n==========================================================\n";
    cout << "FAIR COMPARISON: Single Octave Session vs C++ Library\n";
    cout << "==========================================================\n\n";
    cout << "Test,Method,MaxAbsDiff,Octave_ms,Lib_ms,Speedup\n";
    for (auto &r : rows) {
        double speedup = r.t_oct / r.t_lib;
        cout << r.test << "," << r.method << "," << r.diff << ","
             << r.t_oct << "," << r.t_lib << "," << speedup << "x\n";
    }

    double total_oct_ms = 0, total_lib_ms = 0;
    for (auto &r : rows) {
        total_oct_ms += r.t_oct;
        total_lib_ms += r.t_lib;
    }
    double total_speedup = total_oct_ms / total_lib_ms;

    cout << "\n==========================================================\n";
    cout << "SUMMARY\n";
    cout << "Total Octave (tic/toc): " << total_oct_ms << " ms\n";
    cout << "Total Library: " << total_lib_ms << " ms\n";
    cout << "Overall Speedup: " << total_speedup << "x\n";
    cout << "Octave Session Startup + I/O: " << octave_total_ms - total_oct_ms << " ms\n";
    cout << "==========================================================\n";

    return 0;
}
