#include <gtest/gtest.h>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <string>
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

static void run_octave_solve(const string &matfile, const string &vecfile, const string &outfile, const string &method = "backslash") {
    // Build an inline Octave command string and invoke with --eval to avoid temporary files
    std::ostringstream eval;
    eval << "A = dlmread('" << matfile << "');";
    eval << "b = dlmread('" << vecfile << "');";

    if (method == "chol") {
        eval << "R = chol(A);";
        eval << "y = R' \\ b;";
        eval << "x = R \\ y;";
    } else if (method == "lu") {
        eval << "[L,U,P] = lu(A);";
        eval << "x = U \\ (L \\ (P*b));";
    } else if (method == "gaussseidel") {
        // Inline Gauss-Seidel iteration
        eval << "n = length(b); x = zeros(n,1); tol = 1e-8; maxIter=10000;";
        eval << "for iter=1:maxIter; x_old = x;";
        eval << " for i=1:n; s = 0; if i>1, s = s + A(i,1:i-1)*x(1:i-1); end;";
        eval << " if i<n, s = s + A(i,i+1:n)*x_old(i+1:n); end;";
        eval << " x(i) = (b(i) - s)/A(i,i); end;";
        eval << " if norm(x - x_old, Inf) < tol, break; end; end;";
    } else if (method == "gaussjordan") {
        eval << "M = [A b]; R = rref(M); x = R(:, end);";
    } else {
        eval << "x = A \\ b;";
    }

    eval << "dlmwrite('" << outfile << "', x, 'delimiter', ' ');";

    // Wrap eval string in double quotes for shell, filenames use single quotes inside
    std::string cmd = "octave --quiet --eval \"" + eval.str() + "\" > /dev/null 2>&1";
    int rc = system(cmd.c_str());
    (void)rc;
}

static bool vectors_close(const Vector &a, const Vector &b, double tol) {
    if (a.getLength() != b.getLength()) return false;
    for (int i = 0; i < a.getLength(); ++i) {
        if (fabs(a.getValue(i) - b.getValue(i)) > tol) return false;
    }
    return true;
}

// helper removed: use unified solve(A,b) returning Vector where available

TEST(CompareWithOctave, ConjugateGradientSmallSystems) {
    // 3x3 symmetric positive definite
    Matriz A(3,3);
    A.setValue(0,0,4); A.setValue(0,1,1); A.setValue(0,2,0);
    A.setValue(1,0,1); A.setValue(1,1,3); A.setValue(1,2,1);
    A.setValue(2,0,0); A.setValue(2,1,1); A.setValue(2,2,2);

    Vector b(3);
    b.setValue(0,1); b.setValue(1,2); b.setValue(2,3);

    // Write files
    const string matfile = "matrix_cmp.txt";
    const string vecfile = "vector_cmp.txt";
    const string outfile = "solution_cmp.txt";
    write_matrix_to_file(A, matfile);
    write_vector_to_file(b, vecfile);

    // Run Octave
    run_octave_solve(matfile, vecfile, outfile);

    // Read octave solution
    Vector x_oct = read_vector_from_file(outfile);

    // Solve with library (conjugate gradient)
    Vector x_lib = conjugate_gradient::solve(A, b, 1e-12, 1000);

    ASSERT_TRUE(vectors_close(x_lib, x_oct, 1e-8));
}

TEST(CompareWithOctave, ConjugateGradientFiveByFive) {
    const int n = 5;
    Matriz A(n,n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double val = (i==j) ? 4.0 : ((abs(i-j)==1) ? 1.0 : 0.0);
            A.setValue(i,j,val);
        }
    }
    Vector b(n);
    for (int i = 0; i < n; ++i) b.setValue(i, 1.0 + i);

    const string matfile = "matrix_cmp2.txt";
    const string vecfile = "vector_cmp2.txt";
    const string outfile = "solution_cmp2.txt";
    write_matrix_to_file(A, matfile);
    write_vector_to_file(b, vecfile);
    run_octave_solve(matfile, vecfile, outfile);
    Vector x_oct = read_vector_from_file(outfile);
    Vector x_lib = conjugate_gradient::solve(A, b, 1e-12, 1000);
    ASSERT_TRUE(vectors_close(x_lib, x_oct, 1e-8));
}

TEST(CompareWithOctave, CholeskyMethod) {
    Matriz A(3,3);
    A.setValue(0,0,4); A.setValue(0,1,1); A.setValue(0,2,0);
    A.setValue(1,0,1); A.setValue(1,1,3); A.setValue(1,2,1);
    A.setValue(2,0,0); A.setValue(2,1,1); A.setValue(2,2,2);
    Vector b(3); b.setValue(0,1); b.setValue(1,2); b.setValue(2,3);

    const string matfile = "matrix_chol.txt";
    const string vecfile = "vector_chol.txt";
    const string outfile = "solution_chol.txt";
    write_matrix_to_file(A, matfile);
    write_vector_to_file(b, vecfile);
    run_octave_solve(matfile, vecfile, outfile, "chol");
    Vector x_oct = read_vector_from_file(outfile);
    Vector x_lib = cholesky::solve(A, b);
    ASSERT_TRUE(vectors_close(x_lib, x_oct, 1e-8));
}

TEST(CompareWithOctave, LUmethod) {
    Matriz A(3,3);
    A.setValue(0,0,2); A.setValue(0,1,1); A.setValue(0,2,-1);
    A.setValue(1,0,-3); A.setValue(1,1,4); A.setValue(1,2,2);
    A.setValue(2,0,1); A.setValue(2,1,-1); A.setValue(2,2,3);
    Vector b(3); b.setValue(0,1); b.setValue(1,7); b.setValue(2,3);

    const string matfile = "matrix_lu.txt";
    const string vecfile = "vector_lu.txt";
    const string outfile = "solution_lu.txt";
    write_matrix_to_file(A, matfile);
    write_vector_to_file(b, vecfile);
    run_octave_solve(matfile, vecfile, outfile, "lu");
    Vector x_oct = read_vector_from_file(outfile);

    Vector x_lib = lu::solve(A, b);
    ASSERT_TRUE(vectors_close(x_lib, x_oct, 1e-8));
}

TEST(CompareWithOctave, GaussJordanMethod) {
    Matriz A(3,3);
    A.setValue(0,0,3); A.setValue(0,1,2); A.setValue(0,2,-1);
    A.setValue(1,0,2); A.setValue(1,1,-2); A.setValue(1,2,4);
    A.setValue(2,0,-1); A.setValue(2,1,0.5); A.setValue(2,2,-1);
    Vector b(3); b.setValue(0,1); b.setValue(1,-2); b.setValue(2,0);

    const string matfile = "matrix_gj.txt";
    const string vecfile = "vector_gj.txt";
    const string outfile = "solution_gj.txt";
    write_matrix_to_file(A, matfile);
    write_vector_to_file(b, vecfile);
    run_octave_solve(matfile, vecfile, outfile, "gaussjordan");
    Vector x_oct = read_vector_from_file(outfile);

    Vector x_lib = gaussjordan::solve(A, b);
    ASSERT_TRUE(vectors_close(x_lib, x_oct, 1e-8));
}

TEST(CompareWithOctave, GaussSeidelMethod) {
    const int n = 5;
    Matriz A(n,n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double val = (i==j) ? 4.0 : ((abs(i-j)==1) ? 1.0 : 0.0);
            A.setValue(i,j,val);
        }
    }
    Vector b(n);
    for (int i = 0; i < n; ++i) b.setValue(i, 1.0 + i);

    const string matfile = "matrix_gs.txt";
    const string vecfile = "vector_gs.txt";
    const string outfile = "solution_gs.txt";
    write_matrix_to_file(A, matfile);
    write_vector_to_file(b, vecfile);
    run_octave_solve(matfile, vecfile, outfile, "gaussseidel");
    Vector x_oct = read_vector_from_file(outfile);

    Vector x_lib = gaussseidel::solve(A, b);
    ASSERT_TRUE(vectors_close(x_lib, x_oct, 1e-6));
}

TEST(CompareWithOctave, MultiRHS_CholeskyAndCG) {
        const int n = 4;
        Matriz A(n,n);
        for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                        double val = (i==j) ? 6.0 : ((abs(i-j)==1) ? 1.0 : 0.0);
                        A.setValue(i,j,val);
                }
        }
        Matriz B(n,2);
        for (int i = 0; i < n; ++i) { B.setValue(i,0, 1.0 + i); B.setValue(i,1, 2.0 + 0.5*i); }

        const string matfile = "matrix_mult.txt";
        const string vecfile = "matrixB_mult.txt";
        const string outfile = "solution_mult.txt";
        write_matrix_to_file(A, matfile);
        write_matrix_to_file(B, vecfile);
        run_octave_solve(matfile, vecfile, outfile);
        Matriz X_oct = read_matrix_from_file(outfile);

        Matriz X_chol = cholesky::solve(A, B);
        Matriz X_cg = conjugate_gradient::solve(A, B, 1e-12, 1000);

        // compare elementwise
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < 2; ++j) {
                ASSERT_NEAR(X_chol.getValue(i,j), X_oct.getValue(i,j), 1e-8);
                ASSERT_NEAR(X_cg.getValue(i,j), X_oct.getValue(i,j), 1e-6);
            }
        }
}

TEST(CompareWithOctave, SingularMatrixThrows) {
        Matriz A(3,3);
        A.setValue(0,0,1); A.setValue(0,1,2); A.setValue(0,2,3);
        A.setValue(1,0,2); A.setValue(1,1,4); A.setValue(1,2,6);
        A.setValue(2,0,1); A.setValue(2,1,0); A.setValue(2,2,1);
        Vector b(3); b.setValue(0,1); b.setValue(1,2); b.setValue(2,3);

        Matriz bmat(3,1); for (int i=0;i<3;++i) bmat.setValue(i,0,b.getValue(i));
        EXPECT_THROW(lu::lu(A, bmat), std::runtime_error);
        EXPECT_THROW(gaussjordan::gauss_jordan(A, bmat), std::runtime_error);
}

TEST(CompareWithOctave, IllConditionedHilbert) {
        int n = 6;
        Matriz H(n,n);
        for (int i=0;i<n;++i) for (int j=0;j<n;++j) H.setValue(i,j, 1.0/(i+j+1));
        Vector b(n); for (int i=0;i<n;++i) b.setValue(i, 1.0);

        const string matfile = "matrix_hilbert.txt";
        const string vecfile = "vector_hilbert.txt";
        const string outfile = "solution_hilbert.txt";
        write_matrix_to_file(H, matfile);
        write_vector_to_file(b, vecfile);
        run_octave_solve(matfile, vecfile, outfile);
        Vector x_oct = read_vector_from_file(outfile);

        // CG with preconditioning may help
        solver::IterationInfo info;
        Vector x_cg;
        try {
            x_cg = conjugate_gradient::solve(H, b, 1e-10, 5000, nullptr, &info, true);
        } catch (...) {
            FAIL() << "CG failed on Hilbert";
        }
        // For ill-conditioned systems, compare by residual norm instead of direct
        // component-wise equality.
        Vector r = b;
        Vector Hx = multiplicar(H, x_cg);
        for (int i = 0; i < n; ++i) r.setValue(i, r.getValue(i) - Hx.getValue(i));
        double resnorm = std::sqrt(r.linear_product(r));
        ASSERT_LT(resnorm, 1e-8);
}
