
#include "algebra_linear/matriz.hpp"
#include <stdexcept>

namespace lu {

// Matriz quadrada!!!
Matriz lu(const Matriz &A, const Matriz &b) {
  // Matriz precisa ser quadrada
  if (A.getColumns() != A.getRows()) {
    throw std::runtime_error("matriz precisa ser quadrada");
  }
  auto n = A.getColumns();
  Matriz lu(n, n);
  auto sum = 0;
  for (auto i = 0; i < n; i++) {
    for (auto j = i; j < n; j++) {
      sum = 0;
      for (auto k = 0; k < i; k++)
        sum += lu.getValue(i, k) * lu.getValue(k, j);
      lu.setValue(i, j, A.getValue(i, j) - sum);
    }
    for (auto j = i + 1; j < n; j++) {
      sum = 0;
      for (auto k = 0; k < i; k++) {
        sum += lu.getValue(j, k) * lu.getValue(k, i);
      }
      lu.setValue(j, i,
                  (1 / lu.getValue(i, i)) * A.getValue(j, i) - sum); // ????
    }
  }
  // lu = L+U-I
  //
  // Ly = b
  Matriz y(1, n);
  for (auto i = 0; i < n; i++) {
    sum = 0;
    for (auto k = 0; k < i; k++)
      sum += lu.getValue(i, k) * y.getValue(0, k);
    y.setValue(0, i, b.getValue(0, i) - sum);
  }
  Matriz x(1, n);
  for (auto i = n - 1; i >= 0; i--) {
    sum = 0;
    for (auto k = i + 1; k < i; k++)
      sum += lu.getValue(i, k) * x.getValue(0, k);
    x.setValue(0, i, (1 / lu.getValue(i, i)) * (y.getValue(0, i) - sum));
  }
  return x;
}
} // namespace lu
