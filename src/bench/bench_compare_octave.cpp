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

static void run_octave_solve(const string &matfile, const string &vecfile, const string &outfile, const string &method = "backslash") {
    const string script = "run_octave_tmp.m";
    ofstream s(script);
    s << "A = dlmread('" << matfile << "');\n";
    s << "b = dlmread('" << vecfile << "');\n";

    if (method == "chol") {
        s << "R = chol(A);\n";
        s << "y = R' \\ b;\n";
        s << "x = R \\ y;\n";
    } else if (method == "lu") {
        s << "[L,U,P] = lu(A);\n";
        s << "x = U \\ (L \\ (P*b));\n";
    } else if (method == "gaussseidel") {
        s << "function x = gs(A,b,tol,maxIter)\n";
        s << "n = length(b); x = zeros(n,1);\n";
        s << "for iter=1:maxIter\n";
        s << "  x_old = x;\n";
        s << "  for i=1:n\n";
        s << "    s = A(i,1:i-1)*x(1:i-1) + A(i,i+1:n)*x_old(i+1:n);\n";
        s << "    x(i) = (b(i) - s)/A(i,i);\n";
        s << "  end\n";
        s << "  if norm(x - x_old, Inf) < tol, break; end\n";
        s << "end\n";
        s << "endfunction\n";
        s << "x = gs(A,b,1e-8,10000);\n";
    } else if (method == "gaussjordan") {
        s << "M = [A b];\n";
        s << "R = rref(M);\n";
        s << "x = R(:, end);\n";
    } else {
        s << "x = A \\ b;\n";
    }

    s << "dlmwrite('" << outfile << "', x, 'delimiter', ' ');\n";
    s << "exit;\n";
    s.close();

    string cmd = "octave --quiet " + script + " > /dev/null 2>&1";
    int rc = system(cmd.c_str());
    (void)rc;
}

int main() {
    struct Row { string test; string method; double diff; double t_oct; double t_lib; };
    vector<Row> rows;

    const int runs = 5;

    // 1) CG 3x3
    {
        Matriz A(3,3);
        A.setValue(0,0,4); A.setValue(0,1,1); A.setValue(0,2,0);
        A.setValue(1,0,1); A.setValue(1,1,3); A.setValue(1,2,1);
        A.setValue(2,0,0); A.setValue(2,1,1); A.setValue(2,2,2);
        Vector b(3); b.setValue(0,1); b.setValue(1,2); b.setValue(2,3);
        write_matrix_to_file(A, "matrix_cmp.txt");
        write_vector_to_file(b, "vector_cmp.txt");
        double toct=0, tlib=0; Vector x_oct, x_lib;
        for (int i=0;i<runs;++i) {
            auto t0 = chrono::high_resolution_clock::now();
            run_octave_solve("matrix_cmp.txt","vector_cmp.txt","solution_cmp.txt");
            auto t1 = chrono::high_resolution_clock::now();
            toct += chrono::duration<double, milli>(t1-t0).count();
            x_oct = read_vector_from_file("solution_cmp.txt");

            auto l0 = chrono::high_resolution_clock::now();
            x_lib = conjugate_gradient::solve(A,b,1e-12,1000);
            auto l1 = chrono::high_resolution_clock::now();
            tlib += chrono::duration<double, milli>(l1-l0).count();
        }
        // diff
        double diff = 0; for (int i=0;i<3;++i) diff = max(diff, fabs(x_oct.getValue(i)-x_lib.getValue(i)));
        rows.push_back({"CG_3","conjugate_gradient", diff, toct/runs, tlib/runs});
    }

    // 2) CG 5x5
    {
        const int n=5;
        Matriz A(n,n);
        for (int i=0;i<n;++i) for (int j=0;j<n;++j) {
            double val = (i==j) ? 4.0 : ((abs(i-j)==1) ? 1.0 : 0.0);
            A.setValue(i,j,val);
        }
        Vector b(n); for (int i=0;i<n;++i) b.setValue(i,1.0+i);
        write_matrix_to_file(A, "matrix_cmp2.txt"); write_vector_to_file(b, "vector_cmp2.txt");
        double toct=0,tlib=0; Vector x_oct,x_lib;
        for (int i=0;i<runs;++i){
            auto t0 = chrono::high_resolution_clock::now();
            run_octave_solve("matrix_cmp2.txt","vector_cmp2.txt","solution_cmp2.txt");
            auto t1 = chrono::high_resolution_clock::now(); toct+=chrono::duration<double,milli>(t1-t0).count();
            x_oct = read_vector_from_file("solution_cmp2.txt");
            auto l0 = chrono::high_resolution_clock::now();
            x_lib = conjugate_gradient::solve(A,b,1e-12,1000);
            auto l1 = chrono::high_resolution_clock::now(); tlib+=chrono::duration<double,milli>(l1-l0).count();
        }
        double diff=0; for (int i=0;i<n;++i) diff=max(diff,fabs(x_oct.getValue(i)-x_lib.getValue(i)));
        rows.push_back({"CG_5","conjugate_gradient", diff, toct/runs, tlib/runs});
    }

    // 3) Cholesky 3x3
    {
        Matriz A(3,3);
        A.setValue(0,0,4); A.setValue(0,1,1); A.setValue(0,2,0);
        A.setValue(1,0,1); A.setValue(1,1,3); A.setValue(1,2,1);
        A.setValue(2,0,0); A.setValue(2,1,1); A.setValue(2,2,2);
        Vector b(3); b.setValue(0,1); b.setValue(1,2); b.setValue(2,3);
        write_matrix_to_file(A, "matrix_chol.txt"); write_vector_to_file(b, "vector_chol.txt");
        double toct=0,tlib=0; Vector x_oct,x_lib;
        for (int i=0;i<runs;++i){
            auto t0=chrono::high_resolution_clock::now(); run_octave_solve("matrix_chol.txt","vector_chol.txt","solution_chol.txt","chol"); auto t1=chrono::high_resolution_clock::now(); toct+=chrono::duration<double,milli>(t1-t0).count();
            x_oct = read_vector_from_file("solution_chol.txt");
            auto l0=chrono::high_resolution_clock::now(); x_lib = cholesky::solve(A,b); auto l1=chrono::high_resolution_clock::now(); tlib+=chrono::duration<double,milli>(l1-l0).count();
        }
        double diff=0; for (int i=0;i<3;++i) diff=max(diff,fabs(x_oct.getValue(i)-x_lib.getValue(i)));
        rows.push_back({"Chol_3","cholesky", diff, toct/runs, tlib/runs});
    }

    // 4) LU 3x3
    {
        Matriz A(3,3);
        A.setValue(0,0,2); A.setValue(0,1,1); A.setValue(0,2,-1);
        A.setValue(1,0,-3); A.setValue(1,1,4); A.setValue(1,2,2);
        A.setValue(2,0,1); A.setValue(2,1,-1); A.setValue(2,2,3);
        Vector b(3); b.setValue(0,1); b.setValue(1,7); b.setValue(2,3);
        write_matrix_to_file(A,"matrix_lu.txt"); write_vector_to_file(b,"vector_lu.txt");
        double toct=0,tlib=0; Vector x_oct,x_lib;
        for (int i=0;i<runs;++i){ auto t0=chrono::high_resolution_clock::now(); run_octave_solve("matrix_lu.txt","vector_lu.txt","solution_lu.txt","lu"); auto t1=chrono::high_resolution_clock::now(); toct+=chrono::duration<double,milli>(t1-t0).count(); x_oct=read_vector_from_file("solution_lu.txt"); auto l0=chrono::high_resolution_clock::now(); x_lib = lu::solve(A,b); auto l1=chrono::high_resolution_clock::now(); tlib+=chrono::duration<double,milli>(l1-l0).count(); }
        double diff=0; for (int i=0;i<3;++i) diff=max(diff,fabs(x_oct.getValue(i)-x_lib.getValue(i)));
        rows.push_back({"LU_3","lu", diff, toct/runs, tlib/runs});
    }

    // 5) Gauss-Jordan 3x3
    {
        Matriz A(3,3);
        A.setValue(0,0,3); A.setValue(0,1,2); A.setValue(0,2,-1);
        A.setValue(1,0,2); A.setValue(1,1,-2); A.setValue(1,2,4);
        A.setValue(2,0,-1); A.setValue(2,1,0.5); A.setValue(2,2,-1);
        Vector b(3); b.setValue(0,1); b.setValue(1,-2); b.setValue(2,0);
        write_matrix_to_file(A,"matrix_gj.txt"); write_vector_to_file(b,"vector_gj.txt");
        double toct=0,tlib=0; Vector x_oct,x_lib;
        for (int i=0;i<runs;++i){ auto t0=chrono::high_resolution_clock::now(); run_octave_solve("matrix_gj.txt","vector_gj.txt","solution_gj.txt","gaussjordan"); auto t1=chrono::high_resolution_clock::now(); toct+=chrono::duration<double,milli>(t1-t0).count(); x_oct=read_vector_from_file("solution_gj.txt"); auto l0=chrono::high_resolution_clock::now(); x_lib = gaussjordan::solve(A,b); auto l1=chrono::high_resolution_clock::now(); tlib+=chrono::duration<double,milli>(l1-l0).count(); }
        double diff=0; for (int i=0;i<3;++i) diff=max(diff,fabs(x_oct.getValue(i)-x_lib.getValue(i)));
        rows.push_back({"GJ_3","gaussjordan", diff, toct/runs, tlib/runs});
    }

    // 6) Gauss-Seidel 5x5
    {
        const int n=5; Matriz A(n,n); for (int i=0;i<n;++i) for (int j=0;j<n;++j) { double val = (i==j)?4.0:((abs(i-j)==1)?1.0:0.0); A.setValue(i,j,val); }
        Vector b(n); for (int i=0;i<n;++i) b.setValue(i,1.0+i);
        write_matrix_to_file(A,"matrix_gs.txt"); write_vector_to_file(b,"vector_gs.txt");
        double toct=0,tlib=0; Vector x_oct,x_lib;
        for (int i=0;i<runs;++i){ auto t0=chrono::high_resolution_clock::now(); run_octave_solve("matrix_gs.txt","vector_gs.txt","solution_gs.txt","gaussseidel"); auto t1=chrono::high_resolution_clock::now(); toct+=chrono::duration<double,milli>(t1-t0).count(); x_oct=read_vector_from_file("solution_gs.txt"); auto l0=chrono::high_resolution_clock::now(); x_lib = gaussseidel::solve(A,b); auto l1=chrono::high_resolution_clock::now(); tlib+=chrono::duration<double,milli>(l1-l0).count(); }
        double diff=0; for (int i=0;i<n;++i) diff=max(diff,fabs(x_oct.getValue(i)-x_lib.getValue(i)));
        rows.push_back({"GS_5","gaussseidel", diff, toct/runs, tlib/runs});
    }

    // 7) MultiRHS 4x2
    {
        const int n=4; Matriz A(n,n); for (int i=0;i<n;++i) for (int j=0;j<n;++j) { double val=(i==j)?6.0:((abs(i-j)==1)?1.0:0.0); A.setValue(i,j,val); }
        Matriz B(n,2); for (int i=0;i<n;++i) { B.setValue(i,0,1.0+i); B.setValue(i,1,2.0+0.5*i); }
        write_matrix_to_file(A,"matrix_mult.txt"); write_matrix_to_file(B,"matrixB_mult.txt");
        double toct=0,tlib_ch=0,tlib_cg=0; Matriz X_oct,X_ch,X_cg;
        for (int i=0;i<runs;++i){ auto t0=chrono::high_resolution_clock::now(); run_octave_solve("matrix_mult.txt","matrixB_mult.txt","solution_mult.txt"); auto t1=chrono::high_resolution_clock::now(); toct+=chrono::duration<double,milli>(t1-t0).count(); X_oct = read_matrix_from_file("solution_mult.txt"); auto l0=chrono::high_resolution_clock::now(); X_ch = cholesky::solve(A,B); auto l1=chrono::high_resolution_clock::now(); tlib_ch+=chrono::duration<double,milli>(l1-l0).count(); auto c0=chrono::high_resolution_clock::now(); X_cg = conjugate_gradient::solve(A,B,1e-12,1000); auto c1=chrono::high_resolution_clock::now(); tlib_cg+=chrono::duration<double,milli>(c1-c0).count(); }
        double diff_ch=0, diff_cg=0; for (int i=0;i<n;++i) for (int j=0;j<2;++j) { diff_ch=max(diff_ch,fabs(X_oct.getValue(i,j)-X_ch.getValue(i,j))); diff_cg=max(diff_cg,fabs(X_oct.getValue(i,j)-X_cg.getValue(i,j))); }
        rows.push_back({"Multi_4x2_chol","cholesky", diff_ch, toct/runs, tlib_ch/runs});
        rows.push_back({"Multi_4x2_cg","conjugate_gradient", diff_cg, toct/runs, tlib_cg/runs});
    }

    // 8) Hilbert 6
    {
        int n=6; Matriz H(n,n); for (int i=0;i<n;++i) for (int j=0;j<n;++j) H.setValue(i,j,1.0/(i+j+1)); Vector b(n); for (int i=0;i<n;++i) b.setValue(i,1.0);
        write_matrix_to_file(H,"matrix_hilbert.txt"); write_vector_to_file(b,"vector_hilbert.txt");
        double toct=0,tlib=0; Vector x_oct,x_lib; for (int i=0;i<runs;++i){ auto t0=chrono::high_resolution_clock::now(); run_octave_solve("matrix_hilbert.txt","vector_hilbert.txt","solution_hilbert.txt"); auto t1=chrono::high_resolution_clock::now(); toct+=chrono::duration<double,milli>(t1-t0).count(); x_oct=read_vector_from_file("solution_hilbert.txt"); auto l0=chrono::high_resolution_clock::now(); solver::IterationInfo info; x_lib = conjugate_gradient::solve(H,b,1e-10,5000,nullptr,&info,true); auto l1=chrono::high_resolution_clock::now(); tlib+=chrono::duration<double,milli>(l1-l0).count(); }
        double diff=0; for (int i=0;i<n;++i) diff=max(diff,fabs(x_oct.getValue(i)-x_lib.getValue(i)));
        rows.push_back({"Hilbert_6","conjugate_gradient_precond", diff, toct/runs, tlib/runs});
    }

    // print table
    cout << "Test,Method,MaxAbsDiff,Octave_ms,Lib_ms\n";
    for (auto &r: rows) {
        cout << r.test << "," << r.method << "," << r.diff << "," << r.t_oct << "," << r.t_lib << "\n";
    }

    return 0;
}
