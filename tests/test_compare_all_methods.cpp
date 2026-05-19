#include <gtest/gtest.h>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>

#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"
#include "metodos/cholesky.hpp"
#include "metodos/gauss_jordan.hpp"
#include "metodos/gauss_seidel.hpp"
#include "metodos/lu.hpp"

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

static void write_vector_to_file_col(const Vector &v, const string &path) {
    ofstream f(path);
    for (int i = 0; i < v.getLength(); ++i) {
        f << v.getValue(i) << "\n";
    }
}

static void write_vector_to_file_row_from_matriz(const Matriz &b, const string &path) {
    ofstream f(path);
    // write as row vector
    int cols = b.getColumns();
    int rows = b.getRows();
    if (rows == 1) {
        for (int j = 0; j < cols; ++j) {
            f << b.getValue(0, j);
            if (j+1<cols) f << " ";
        }
        f << "\n";
    } else if (cols == 1) {
        for (int i = 0; i < rows; ++i) {
            f << b.getValue(i, 0);
            if (i+1<rows) f << "\n";
        }
        f << "\n";
    }
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

static void run_octave_solve(const string &matfile, const string &vecfile, const string &outfile) {
    const string script = "run_octave_tmp_all.m";
    ofstream s(script);
    s << "A = dlmread('" << matfile << "');\n";
    s << "b = dlmread('" << vecfile << "');\n";
    s << "x = A \\ b;\n";
    s << "dlmwrite('" << outfile << "', x, 'delimiter', ' ');\n";
    s << "exit;\n";
    s.close();
    string cmd = "octave --quiet " + script + " > /dev/null 2>&1";
    system(cmd.c_str());
}

static bool vectors_close(const Vector &a, const Vector &b, double tol) {
    if (a.getLength() != b.getLength()) return false;
    for (int i = 0; i < a.getLength(); ++i) {
        if (fabs(a.getValue(i) - b.getValue(i)) > tol) return false;
    }
    return true;
}

static Vector matriz_to_vector_col(const Matriz &m) {
    // if m is n x 1 => column; if 1 x n => convert to column
    if (m.getColumns() == 1) {
        Vector v(m.getRows());
        for (int i = 0; i < m.getRows(); ++i) v.setValue(i, m.getValue(i, 0));
        return v;
    } else if (m.getRows() == 1) {
        Vector v(m.getColumns());
        for (int j = 0; j < m.getColumns(); ++j) v.setValue(j, m.getValue(0, j));
        return v;
    }
    // fallback: flatten first column
    Vector v(m.getRows());
    for (int i = 0; i < m.getRows(); ++i) v.setValue(i, m.getValue(i, 0));
    return v;
}

TEST(CompareAllMethods, Cholesky) {
    Matriz A(3,3);
    A.setValue(0,0,4); A.setValue(0,1,1); A.setValue(0,2,1);
    A.setValue(1,0,1); A.setValue(1,1,3); A.setValue(1,2,1);
    A.setValue(2,0,1); A.setValue(2,1,1); A.setValue(2,2,2);

    Vector b(3); b.setValue(0,1); b.setValue(1,2); b.setValue(2,3);

    write_matrix_to_file(A, "mat_chol.txt");
    write_vector_to_file_col(b, "vec_chol.txt");
    run_octave_solve("mat_chol.txt","vec_chol.txt","sol_chol.txt");
    Vector x_oct = read_vector_from_file("sol_chol.txt");

    Vector x_lib = cholesky::solve(A, b);
    ASSERT_TRUE(vectors_close(x_lib, x_oct, 1e-8));
}

TEST(CompareAllMethods, GaussJordan) {
    Matriz A(3,3);
    A.setValue(0,0,3); A.setValue(0,1,2); A.setValue(0,2,-1);
    A.setValue(1,0,2); A.setValue(1,1,-2); A.setValue(1,2,4);
    A.setValue(2,0,-1); A.setValue(2,1,0.5); A.setValue(2,2,-1);

    Matriz b(3,1);
    b.setValue(0,0,1); b.setValue(1,0,-2); b.setValue(2,0,0);

    write_matrix_to_file(A, "mat_gj.txt");
    write_matrix_to_file(b, "vec_gj.txt");
    run_octave_solve("mat_gj.txt","vec_gj.txt","sol_gj.txt");
    Vector x_oct = read_vector_from_file("sol_gj.txt");

    Matriz xm = gaussjordan::gauss_jordan(A, b);
    Vector x_lib = matriz_to_vector_col(xm);
    ASSERT_TRUE(vectors_close(x_lib, x_oct, 1e-8));
}

TEST(CompareAllMethods, GaussSeidel) {
    Matriz A(3,3);
    A.setValue(0,0,4); A.setValue(0,1,-1); A.setValue(0,2,0);
    A.setValue(1,0,-1); A.setValue(1,1,4); A.setValue(1,2,-1);
    A.setValue(2,0,0); A.setValue(2,1,-1); A.setValue(2,2,4);

    Matriz b(3,1);
    b.setValue(0,0,1); b.setValue(1,0,2); b.setValue(2,0,3);

    write_matrix_to_file(A, "mat_gs.txt");
    write_matrix_to_file(b, "vec_gs.txt");
    run_octave_solve("mat_gs.txt","vec_gs.txt","sol_gs.txt");
    Vector x_oct = read_vector_from_file("sol_gs.txt");

    Matriz xm = gaussseidel::gauss_seidel(A, b);
    Vector x_lib = matriz_to_vector_col(xm);
    ASSERT_TRUE(vectors_close(x_lib, x_oct, 1e-6));
}

TEST(CompareAllMethods, LU) {
    Matriz A(3,3);
    A.setValue(0,0,2); A.setValue(0,1,1); A.setValue(0,2,1);
    A.setValue(1,0,4); A.setValue(1,1,-6); A.setValue(1,2,0);
    A.setValue(2,0,-2); A.setValue(2,1,7); A.setValue(2,2,2);

    // lu implementation expects b as 1 x n in this codebase
    Matriz b(1,3);
    b.setValue(0,0,5); b.setValue(0,1,-2); b.setValue(0,2,9);

    // For Octave, write b as column
    Vector bcol(3); bcol.setValue(0,5); bcol.setValue(1,-2); bcol.setValue(2,9);
    write_matrix_to_file(A, "mat_lu.txt");
    write_vector_to_file_col(bcol, "vec_lu.txt");
    run_octave_solve("mat_lu.txt","vec_lu.txt","sol_lu.txt");
    Vector x_oct = read_vector_from_file("sol_lu.txt");

    Matriz xm = lu::lu(A, b);
    Vector x_lib = matriz_to_vector_col(xm);
    ASSERT_TRUE(vectors_close(x_lib, x_oct, 1e-6));
}
