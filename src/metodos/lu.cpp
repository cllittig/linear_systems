
#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"
#include <stdexcept>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

namespace lu {

Matriz lu(const Matriz &A, const Matriz &b) {
  if (A.getColumns() != A.getRows())
    throw std::runtime_error("matriz precisa ser quadrada");

  int n = A.getRows();
  int bm_rows = b.getRows();
  int bm_cols = b.getColumns();

  if (bm_rows != n && !(bm_rows == 1 && bm_cols == n))
    throw std::runtime_error("Dimensões incompatíveis entre A e b.");

  int rhs = (bm_rows == n) ? bm_cols : 1;
  int cols_aug = n + rhs;

  // Matriz aumentada [A | b]
  Matriz aug(n, cols_aug);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j)
      aug.setValue(i, j, A.getValue(i, j));
    if (bm_rows == n) {
      for (int j = 0; j < bm_cols; ++j)
        aug.setValue(i, n + j, b.getValue(i, j));
    } else {
      aug.setValue(i, n, b.getValue(0, i));
    }
  }

  double global_max = 0.0;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < cols_aug; ++j)
      global_max = std::max(global_max, std::fabs(aug.getValue(i, j)));
  double limiar_singular =
      std::numeric_limits<double>::epsilon() * global_max * n * 10.0;

  // Eliminação gaussiana com pivoteamento parcial
  for (int k = 0; k < n; ++k) {
    int pivo = k;
    double maxv = std::fabs(aug.getValue(k, k));
    for (int i = k + 1; i < n; ++i) {
      double v = std::fabs(aug.getValue(i, k));
      if (v > maxv) { maxv = v; pivo = i; }
    }
    if (maxv < limiar_singular)
      throw std::runtime_error("Matriz singular ou quase singular (criterio relativo).");

    if (pivo != k)
      aug.trocarLinhas(pivo, k);

    double inv_pivo = 1.0 / aug.getValue(k, k);
    for (int i = k + 1; i < n; ++i) {
      double fator = aug.getValue(i, k) * inv_pivo;
      aug.setValue(i, k, 0.0);
      aug.linhaAxpy(i, k, -fator, k + 1);
    }
  }

  // Substituição retroativa
  Matriz X(n, rhs);
  for (int col = 0; col < rhs; ++col) {
    std::vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
      double s = aug.getValue(i, n + col);
      for (int j = i + 1; j < n; ++j) s -= aug.getValue(i, j) * x[j];
      x[i] = s / aug.getValue(i, i);
    }
    for (int i = 0; i < n; ++i) X.setValue(i, col, x[i]);
  }

  if (bm_rows == 1 && bm_cols == n) {
    Matriz out(1, n);
    for (int j = 0; j < n; ++j) out.setValue(0, j, X.getValue(j, 0));
    return out;
  }
  return X;
}

} // namespace lu

namespace lu {
Vector solve(const Matriz &A, const Vector &b) {
  if (A.getRows() != b.getLength())
    throw std::runtime_error("Dimensoes incompativeis entre A e b.");
  Matriz bmat(A.getRows(), 1);
  for (int i = 0; i < b.getLength(); ++i) bmat.setValue(i, 0, b.getValue(i));
  Matriz xmat = lu(A, bmat);
  Vector x(b.getLength());
  for (int i = 0; i < x.getLength(); ++i) x.setValue(i, xmat.getValue(i, 0));
  return x;
}
} // namespace lu
