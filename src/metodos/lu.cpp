
#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"
#include <stdexcept>
#include <vector>
#include <cmath>
#include <algorithm>

namespace lu {

Matriz lu(const Matriz &A, const Matriz &b) {
  if (A.getColumns() != A.getRows()) {
    throw std::runtime_error("matriz precisa ser quadrada");
  }
<<<<<<< HEAD
  int n = A.getRows();

  // determine b shape
  int bm_rows = b.getRows();
  int bm_cols = b.getColumns();
  int m = 0;
  bool b_is_column = false;
  if (bm_rows == n) {
    m = bm_cols;
    b_is_column = true;
  } else if (bm_rows == 1 && bm_cols == n) {
    m = 1;
    b_is_column = false; // row vector provided, treat as one RHS with entries across columns
  } else {
    throw std::runtime_error("Dimensões incompatíveis entre A e b.");
  }

  int rhs = (bm_rows == n) ? bm_cols : 1;
  int cols_aug = n + rhs;

  // build augmented matrix (n x (n+rhs)) as double vector
  std::vector<std::vector<double>> aug(n, std::vector<double>(cols_aug, 0.0));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) aug[i][j] = A.getValue(i, j);
    if (bm_rows == n) {
      for (int j = 0; j < bm_cols; ++j) aug[i][n + j] = b.getValue(i, j);
    } else {
      // b is 1 x n row, interpret as column vector
      aug[i][n + 0] = b.getValue(0, i);
=======
  auto n = A.getColumns();
  Matriz lu(n, n);
  double sum = 0.0;
  for (auto i = 0; i < n; i++) {
    for (auto j = i; j < n; j++) {
      sum = 0;
      for (auto k = 0; k < i; k++)
        sum += lu.getValue(i, k) * lu.getValue(k, j);
      lu.setValue(i, j, A.getValue(i, j) - sum);
>>>>>>> a0306ec7ad214946c9df52e0635fa4c496afe0d2
    }
  }

  // Gaussian elimination with partial pivoting
  for (int k = 0; k < n; ++k) {
    // find pivot
    int piv = k;
    double maxv = std::fabs(aug[k][k]);
    for (int i = k + 1; i < n; ++i) {
      if (std::fabs(aug[i][k]) > maxv) {
        maxv = fabs(aug[i][k]);
        piv = i;
      }
<<<<<<< HEAD
    }
    if (fabs(maxv) < 1e-12) throw std::runtime_error("Matriz singular ou quase singular.");
    if (piv != k) std::swap(aug[piv], aug[k]);

    for (int i = k + 1; i < n; ++i) {
      double factor = aug[i][k] / aug[k][k];
      aug[i][k] = 0.0;
      for (int j = k + 1; j < cols_aug; ++j) aug[i][j] -= factor * aug[k][j];
=======
      lu.setValue(j, i, (A.getValue(j, i) - sum) / lu.getValue(i, i));
>>>>>>> a0306ec7ad214946c9df52e0635fa4c496afe0d2
    }
  }

  // back substitution for each RHS
  Matriz X(n, rhs);
  for (int col = 0; col < rhs; ++col) {
    std::vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
      double s = aug[i][n + col];
      for (int j = i + 1; j < n; ++j) s -= aug[i][j] * x[j];
      x[i] = s / aug[i][i];
    }
    for (int i = 0; i < n; ++i) X.setValue(i, col, x[i]);
  }
<<<<<<< HEAD

  // if original b was 1 x n row (bm_rows==1 && bm_cols==n), return 1 x n row as before
  if (bm_rows == 1 && bm_cols == n) {
    Matriz out(1, n);
    for (int j = 0; j < n; ++j) out.setValue(0, j, X.getValue(j, 0));
    return out;
=======
  Matriz x(1, n);
  for (auto i = n - 1; i >= 0; i--) {
    sum = 0;
    for (auto k = i + 1; k < n; k++)
      sum += lu.getValue(i, k) * x.getValue(0, k);
    x.setValue(0, i, (1 / lu.getValue(i, i)) * (y.getValue(0, i) - sum));
>>>>>>> a0306ec7ad214946c9df52e0635fa4c496afe0d2
  }

  return X;
}

} // namespace lu

namespace lu {
Vector solve(const Matriz &A, const Vector &b) {
  if (A.getRows() != b.getLength()) {
    throw std::runtime_error("Dimensoes incompativeis entre A e b.");
  }
  Matriz bmat(A.getRows(), 1);
  for (int i = 0; i < b.getLength(); ++i) bmat.setValue(i, 0, b.getValue(i));
  Matriz xmat = lu(A, bmat);
  Vector x(b.getLength());
  for (int i = 0; i < x.getLength(); ++i) x.setValue(i, xmat.getValue(i, 0));
  return x;
}

Vector solve(const Matriz &A, const Vector &b) {
  int n = b.getLength();
  Matriz b_row(1, n);
  for (int i = 0; i < n; i++)
    b_row.setValue(0, i, b.getValue(i));
  Matriz x_row = lu(A, b_row);
  Vector x(n);
  for (int i = 0; i < n; i++)
    x.setValue(i, x_row.getValue(0, i));
  return x;
}

} // namespace lu
