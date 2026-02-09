
#include "algebra_linear/matriz.hpp"
#include <stdexcept>

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

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      augmented.setValue(i, j, A.getValue(i, j));
    }
    for (int j = 0; j < m; j++) {
      augmented.setValue(i, j + n, b.getValue(i, j));
    }
  }

  for (int i = 0; i < n; i++) {
    double pivot = augmented.getValue(i, i);

    if (std::abs(pivot) < 1e-9) {
      bool swapped = false;
      for (int k = i + 1; k < n; k++) {
        if (std::abs(augmented.getValue(k, i)) > 1e-9) {
          augmented.trocarLinhas(i, k);
          swapped = true;
          break;
        }
      }
      if (!swapped) {
        throw std::runtime_error("O sistema não possui solução única.");
      }
      pivot = augmented.getValue(i, i);
    }

    for (int j = 0; j < n + m; j++) {
      augmented.setValue(i, j, augmented.getValue(i, j) / pivot);
    }

    // Eliminação
    for (int k = 0; k < n; k++) {
      if (k == i)
        continue;

      double factor = augmented.getValue(k, i);
      for (int j = 0; j < n + m; j++) {
        augmented.setValue(
            k, j, augmented.getValue(k, j) - factor * augmented.getValue(i, j));
      }
    }
  }

  Matriz x(n, m);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      x.setValue(i, j, augmented.getValue(i, j + n));
    }
  }

  return x;
}

}; // namespace gaussjordan
