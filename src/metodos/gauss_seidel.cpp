
#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"
#include <cmath>
#include <stdexcept>
#include "metodos/solver_utils.hpp"
#include "metodos/gauss_seidel.hpp"

namespace gaussseidel {

Matriz gauss_seidel(const Matriz &A, const Matriz &b, double tol, int maxIter, solver::IterationInfo* info, double omega) {
  if (!A.isSquare())
    throw std::invalid_argument("A matriz A deve ser quadrada.");

  if (!(omega > 0.0 && omega < 2.0))
    throw std::invalid_argument("Parametro omega invalido: deve ser (0,2) para SOR.");

  int n = A.getRows();

  if (b.getRows() != n || b.getColumns() != 1)
    throw std::invalid_argument("O vetor b deve ter dimensão (n x 1).");

  // Checagem de pivô relativa por linha (evita limiar absoluto inseguro)
  for (int i = 0; i < n; ++i) {
    double row_max = 0.0;
    for (int j = 0; j < n; ++j) row_max = std::max(row_max, std::abs(A.getValue(i, j)));
    if (row_max == 0.0) throw std::runtime_error("Linha nula detectada na matriz A; sistema degenerado.");
    if (std::abs(A.getValue(i, i)) < 1e-12 * row_max)
      throw std::runtime_error("Pivô na diagonal muito pequeno (possivel singularidade). Reordene o sistema.");
  }

  Matriz x(n, 1);
  for (int i = 0; i < n; ++i)
    x.setValue(i, 0, 0.0);

  // Inicializa info se fornecido
  if (info) { info->iterations = 0; info->final_residual_norm = 0.0; info->converged = false; }

  // Vetor b como Vector para cálculo de resíduo e norma
  Vector bvec(n);
  for (int i = 0; i < n; ++i) bvec.setValue(i, b.getValue(i, 0));
  double norm_b = bvec.euclidian_length();

  if (norm_b == 0.0) {
    if (info) { info->iterations = 0; info->final_residual_norm = 0.0; info->converged = true; }
    return x;
  }

  for (int iter = 0; iter < maxIter; ++iter) {
    for (int i = 0; i < n; ++i) {
      double soma = 0.0;
      for (int j = 0; j < n; ++j) {
        if (j != i)
          soma += A.getValue(i, j) * x.getValue(j, 0);
      }
      double xNovo = (b.getValue(i, 0) - soma) / A.getValue(i, i);
      // relaxation (SOR)
      xNovo = omega * xNovo + (1.0 - omega) * x.getValue(i, 0);
      x.setValue(i, 0, xNovo);
    }

    // Calcular resíduo relativo ||b - A x|| / ||b||
    Vector xvec(n);
    for (int i = 0; i < n; ++i) xvec.setValue(i, x.getValue(i, 0));
    Vector Ax = multiplicar(A, xvec);
    Vector residual = Ax - bvec;
    double rnorm = residual.euclidian_length();

    if (rnorm <= tol * norm_b) {
      if (info) { info->iterations = iter + 1; info->final_residual_norm = rnorm; info->converged = true; }
      return x; // convergiu
    }
    if (info) info->final_residual_norm = rnorm;
  }

  if (info) { info->iterations = maxIter; info->final_residual_norm = 0.0; info->converged = false; }
  throw std::runtime_error("Não convergiu no número máximo de iterações.");
}

Vector solve(const Matriz &A, const Vector &b, double tol, int maxIter, solver::IterationInfo* info, double omega) {
  if (!A.isSquare())
    throw std::invalid_argument("A matriz A deve ser quadrada.");

  int n = A.getRows();
  if (b.getLength() != n) throw std::invalid_argument("O vetor b deve ter dimensão (n).");

  Matriz bmat(n, 1);
  for (int i = 0; i < n; ++i) bmat.setValue(i, 0, b.getValue(i));
  Matriz xmat = gauss_seidel(A, bmat, tol, maxIter, info, omega);
  Vector x(n);
  for (int i = 0; i < n; ++i) x.setValue(i, xmat.getValue(i, 0));
  return x;
}

} // namespace gaussseidel
