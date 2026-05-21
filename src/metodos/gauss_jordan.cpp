
#include "algebra_linear/matriz.hpp"
#include <stdexcept>
#include <limits>
#include "algebra_linear/vector.hpp"

namespace gaussjordan {

Matriz gauss_jordan(const Matriz &A, const Matriz &b) {
  if (!A.isSquare()) {
    throw std::runtime_error("A matriz A deve ser quadrada.");
  }

  if (A.getRows() != b.getRows()) {
    throw std::runtime_error("Dimensões incompatíveis entre A e b.");
  }

  int n = A.getRows();
  int m = b.getColumns();

  // [A | b]
  Matriz augmented(n, n + m);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      augmented.setValue(i, j, A.getValue(i, j));
    }
    for (int j = 0; j < m; ++j) {
      augmented.setValue(i, j + n, b.getValue(i, j));
    }
  }

  for (int i = 0; i < n; ++i) {
    double pivot = augmented.getValue(i, i);

    // Relative pivot check: compute row maximum and use machine epsilon scaled threshold
    double row_max = 0.0;
    for (int j = 0; j < n + m; ++j) row_max = std::max(row_max, std::abs(augmented.getValue(i, j)));
    if (row_max == 0.0) {
      throw std::runtime_error("Linha nula detectada — sistema degenerado.");
    }
    double rel_eps = std::numeric_limits<double>::epsilon();
    double pivot_threshold = rel_eps * row_max * 10.0;

    if (std::abs(pivot) < pivot_threshold) {
      bool swapped = false;
      for (int k = i + 1; k < n; ++k) {
        if (std::abs(augmented.getValue(k, i)) > pivot_threshold) {
          augmented.trocarLinhas(i, k);
          swapped = true;
          break;
        }
      }
      if (!swapped) {
        throw std::runtime_error("O sistema não possui solução única (pivô pequeno ou zero).");
      }
      pivot = augmented.getValue(i, i);
    }

    for (int j = 0; j < n + m; ++j) {
      augmented.setValue(i, j, augmented.getValue(i, j) / pivot);
    }

    // Eliminação
    for (int k = 0; k < n; ++k) {
      if (k == i)
        continue;

      double factor = augmented.getValue(k, i);
      for (int j = 0; j < n + m; ++j) {
        augmented.setValue(
            k, j, augmented.getValue(k, j) - factor * augmented.getValue(i, j));
      }
    }
  }

  Matriz x(n, m);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      x.setValue(i, j, augmented.getValue(i, j + n));
    }
  }

  return x;
}

Vector solve(const Matriz &A, const Vector &b) {
  if (A.getRows() != b.getLength()) {
    throw std::runtime_error("Dimensoes incompativeis entre A e b.");
  }
  Matriz bmat(A.getRows(), 1);
  for (int i = 0; i < b.getLength(); ++i) bmat.setValue(i, 0, b.getValue(i));
  Matriz x = gauss_jordan(A, bmat);
  Vector out(b.getLength());
  for (int i = 0; i < out.getLength(); ++i) out.setValue(i, x.getValue(i, 0));
  return out;
}

} // namespace gaussjordan
