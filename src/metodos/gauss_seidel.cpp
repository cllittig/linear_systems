
#include "algebra_linear/matriz.hpp"
#include <cmath>
#include <stdexcept>

namespace gaussseidel {

Matriz gauss_seidel(const Matriz &A, const Matriz &b) {
  double tol = 1e-6;
  int maxIter = 1000;
  if (!A.isSquare())
    throw std::invalid_argument("A matriz A deve ser quadrada.");

  int n = A.getRows();

  if (b.getRows() != n || b.getColumns() != 1)
    throw std::invalid_argument("O vetor b deve ter dimensão (n x 1).");

  for (int i = 0; i < n; i++) {
    if (std::abs(A.getValue(i, i)) < 1e-12)
      throw std::runtime_error("Pivô nulo na diagonal. Reordene o sistema.");
  }

  Matriz x(n, 1);
  for (int i = 0; i < n; i++)
    x.setValue(i, 0, 0.0);

  for (int iter = 0; iter < maxIter; iter++) {
    double maxErro = 0.0;

    for (int i = 0; i < n; i++) {
      double soma = 0.0;

      for (int j = 0; j < n; j++) {
        if (j != i)
          soma += A.getValue(i, j) * x.getValue(j, 0);
      }

      double xNovo = (b.getValue(i, 0) - soma) / A.getValue(i, i);

      // ???
      double erro = std::abs(xNovo - x.getValue(i, 0));
      if (std::abs(xNovo) > 1e-12)
        erro /= std::abs(xNovo);

      if (erro > maxErro)
        maxErro = erro;

      x.setValue(i, 0, xNovo);
    }

    if (maxErro < tol)
      return x; // convergiu
  }

  throw std::runtime_error("Não convergiu no número máximo de iterações.");
}

} // namespace guassseidel
