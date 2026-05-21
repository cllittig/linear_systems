#include <iostream>
#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"
#include "metodos/lu.hpp"
#include "metodos/cholesky.hpp"
#include "metodos/conjugate_gradient.hpp"
#include "metodos/gauss_seidel.hpp"
#include "metodos/gauss_jordan.hpp"
#include "metodos/solver_utils.hpp"

using namespace std;

static void printVector(const Vector &v, const string &label) {
  cout << label << ": [";
  for (int i = 0; i < v.getLength(); ++i) {
    cout << v.getValue(i);
    if (i + 1 < v.getLength()) cout << ", ";
  }
  cout << "]\n";
}

int main() {
  try {
    // Example 3x3 SPD matrix and RHS
    Matriz A(3, 3);
    A.setValue(0,0,4.0); A.setValue(0,1,1.0); A.setValue(0,2,1.0);
    A.setValue(1,0,1.0); A.setValue(1,1,3.0); A.setValue(1,2,0.0);
    A.setValue(2,0,1.0); A.setValue(2,1,0.0); A.setValue(2,2,2.0);

    Vector b(3);
    b.setValue(0, 1.0);
    b.setValue(1, 2.0);
    b.setValue(2, 3.0);

    cout << "Matrix A:\n";
    A.print();
    printVector(b, "Vector b");

    // LU
    try {
      Vector x_lu = lu::solve(A, b);
      printVector(x_lu, "LU solution");
    } catch (const exception &e) {
      cout << "LU failed: " << e.what() << "\n";
    }

    // Cholesky
    try {
      Vector x_chol = cholesky::solve(A, b);
      printVector(x_chol, "Cholesky solution");
    } catch (const exception &e) {
      cout << "Cholesky failed: " << e.what() << "\n";
    }

    // Conjugate Gradient
    try {
      solver::IterationInfo info_cg;
      Vector x_cg = conjugate_gradient::solve(A, b, 1e-8, 1000, nullptr, &info_cg, false);
      printVector(x_cg, "Conjugate Gradient solution");
      cout << "  CG iterations: " << info_cg.iterations << ", residual: " << info_cg.final_residual_norm << "\n";
    } catch (const exception &e) {
      cout << "Conjugate Gradient failed: " << e.what() << "\n";
    }

    // Gauss-Seidel
    try {
      solver::IterationInfo info_gs;
      Vector x_gs = gaussseidel::solve(A, b, 1e-8, 1000, &info_gs, 1.0);
      printVector(x_gs, "Gauss-Seidel solution");
      cout << "  GS iterations: " << info_gs.iterations << ", residual: " << info_gs.final_residual_norm << "\n";
    } catch (const exception &e) {
      cout << "Gauss-Seidel failed: " << e.what() << "\n";
    }

    // Gauss-Jordan
    try {
      Vector x_gj = gaussjordan::solve(A, b);
      printVector(x_gj, "Gauss-Jordan solution");
    } catch (const exception &e) {
      cout << "Gauss-Jordan failed: " << e.what() << "\n";
    }

  } catch (const exception &e) {
    cerr << "Unexpected error: " << e.what() << "\n";
    return 1;
  }

  return 0;
}
