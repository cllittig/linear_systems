#include "metodos/conjugate_gradient.hpp"
#include <cmath>
#include <stdexcept>
#include "metodos/solver_utils.hpp"

namespace conjugate_gradient {

namespace {

void addScaledInPlace(Vector &target, const Vector &direction, double scale) {
  if (target.getLength() != direction.getLength()) {
    throw std::invalid_argument("Vetores com dimensoes incompativeis.");
  }
  for (int i = 0; i < target.getLength(); ++i) {
    target.setValue(i, target.getValue(i) + scale * direction.getValue(i));
  }
}

} // namespace

/**
 * Resolve o sistema linear A*x = b usando o Método do Gradiente Conjugado.
 *
 * @param A             Matriz quadrada e simétrica definida positiva.
 * @param b             Vetor do lado direito.
 * @param tolerance     Critério de parada (norma relativa do resíduo).
 * @param maxIterations Número máximo de iterações permitidas.
 * @param x0            Chute inicial (default: vetor zero).
 * @return              Vetor solução x tal que ||b - A*x|| / ||b|| <=
 * tolerance.
 *
 * @throws std::invalid_argument se os parâmetros forem inválidos.
 * @throws std::runtime_error    se não convergir ou ocorrer falha numérica.
 */
// Helper: retrieve a column vector solution from a Matriz result (if B was empty fallback)
// We'll implement multi-RHS solve and single-RHS by delegation.

// Multi-RHS implementation
Matriz solve(const Matriz &A, const Matriz &B, double tolerance,
             int maxIterations, solver::IterationInfo* info, bool useJacobiPrecond) {
  if (!A.isSquare()) throw std::invalid_argument("A matriz A deve ser quadrada.");
  if (B.getRows() != 0 && A.getRows() != B.getRows())
    throw std::invalid_argument("Dimensoes incompativeis entre A e B.");

  int n = A.getRows();
  int m = (B.getRows() == 0) ? 0 : B.getColumns();
  if (m == 0) {
    // No B provided: this function was called as fallback; return empty
    return Matriz();
  }

  Matriz X(n, m);
  for (int col = 0; col < m; ++col) {
    Vector bcol(n);
    for (int i = 0; i < n; ++i) bcol.setValue(i, B.getValue(i, col));
    // call vector CG with preconditioner flag and optional info ignored per-column
    solver::IterationInfo colInfo;
    Vector x = solve(A, bcol, tolerance, maxIterations, nullptr, &colInfo, useJacobiPrecond);
    for (int i = 0; i < n; ++i) X.setValue(i, col, x.getValue(i));
    if (info) {
      info->iterations += colInfo.iterations;
      info->final_residual_norm = std::max(info->final_residual_norm, colInfo.final_residual_norm);
      info->converged = info->converged && colInfo.converged;
    }
  }
  return X;
}

// New vector solve with diagnostics and optional Jacobi preconditioning
Vector solve(const Matriz &A, const Vector &b, double tolerance,
             int maxIterations, const Vector *x0, solver::IterationInfo* info, bool useJacobiPrecond) {
  if (!A.isSquare()) {
    throw std::invalid_argument("A matriz A deve ser quadrada.");
  }
  if (A.getRows() != b.getLength()) {
    throw std::invalid_argument("Dimensoes incompativeis entre A e b.");
  }
  if (tolerance <= 0.0) throw std::invalid_argument("A tolerancia deve ser positiva.");
  if (maxIterations <= 0) throw std::invalid_argument("maxIterations deve ser positivo.");

  int n = b.getLength();
  Vector x(n);
  if (x0) x = *x0;

  Vector r = b;
  if (x0) {
    Vector Ax0 = multiplicar(A, x);
    addScaledInPlace(r, Ax0, -1.0);
  }

  Vector z(n);
  Vector p(n);
  if (useJacobiPrecond) {
    // Jacobi preconditioner: M = diag(A)
    for (int i = 0; i < n; ++i) {
      double d = A.getValue(i, i);
      if (std::abs(d) < 1e-18) throw std::runtime_error("Jacobi preconditioner: diagonal zero.");
      z.setValue(i, r.getValue(i) / d);
    }
    p = z;
  } else {
    p = r;
  }

  double rzold = (useJacobiPrecond ? r.linear_product(z) : r.linear_product(r));
  double r0norm = std::sqrt(r.linear_product(r));
  if (r0norm <= tolerance) {
    if (info) { info->iterations = 0; info->final_residual_norm = r0norm; info->converged = true; }
    return x;
  }
  double stopThreshold = tolerance * r0norm;

  bool converged = false;
  int iter;
  for (iter = 0; iter < maxIterations; ++iter) {
    Vector Ap = multiplicar(A, p);
    double denom = p.linear_product(Ap);
    if (std::abs(denom) < 1e-18) {
      throw std::runtime_error("Falha numerica no gradiente conjugado: p^T A p ~= 0.");
    }
    double alpha = rzold / denom;
    addScaledInPlace(x, p, alpha);
    addScaledInPlace(r, Ap, -alpha);

    double rnorm = std::sqrt(r.linear_product(r));
    if (rnorm <= stopThreshold) { converged = true; ++iter; break; }

    if (useJacobiPrecond) {
      for (int i = 0; i < n; ++i) z.setValue(i, r.getValue(i) / A.getValue(i, i));
      double rznew = r.linear_product(z);
      double beta = rznew / rzold;
      for (int i = 0; i < n; ++i) p.setValue(i, z.getValue(i) + beta * p.getValue(i));
      rzold = rznew;
    } else {
      double rsnew = r.linear_product(r);
      double beta = rsnew / rzold;
      for (int i = 0; i < n; ++i) p.setValue(i, r.getValue(i) + beta * p.getValue(i));
      rzold = rsnew;
    }
  }

  if (info) {
    info->iterations = iter;
    info->final_residual_norm = std::sqrt(r.linear_product(r));
    info->converged = converged;
  }
  if (!converged) throw std::runtime_error("Gradiente conjugado nao convergiu no numero maximo de iteracoes.");
  return x;
}

} // namespace conjugate_gradient
