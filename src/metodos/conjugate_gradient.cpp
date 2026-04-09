#pragma once
#include "metodos/conjugate_gradient.hpp"
#include <cmath>
#include <stdexcept>

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
 * @return              Vetor solução x tal que ||b - A*x|| / ||b|| <= tolerance.
 *
 * @throws std::invalid_argument se os parâmetros forem inválidos.
 * @throws std::runtime_error    se não convergir ou ocorrer falha numérica.
 */
Vector solve(const Matriz &A, const Vector &b, double tolerance,
             int maxIterations, const Vector *x0 = nullptr) {

  // ------------------------------------------------------------------
  // Validações
  // ------------------------------------------------------------------
  if (!A.isSquare()) {
    throw std::invalid_argument("A matriz A deve ser quadrada.");
  }
  if (!A.isSimetric()) {
    throw std::invalid_argument("A matriz A deve ser simetrica.");
  }
  if (A.getRows() != b.getLength()) {
    throw std::invalid_argument("Dimensoes incompativeis entre A e b.");
  }
  if (tolerance <= 0.0) {
    throw std::invalid_argument("A tolerancia deve ser positiva.");
  }
  if (maxIterations <= 0) {
    throw std::invalid_argument("maxIterations deve ser positivo.");
  }

  const int n = b.getLength();

  // ------------------------------------------------------------------
  // Inicialização: suporta chute inicial x0 arbitrário
  // ------------------------------------------------------------------
  Vector x(n); // x = 0 por padrão
  if (x0 != nullptr) {
    if (x0->getLength() != n) {
      throw std::invalid_argument("x0 tem dimensao incompativel com b.");
    }
    x = *x0;
  }

  // r₀ = b - A·x₀  (se x₀ = 0, r₀ = b)
  Vector r = b;
  if (x0 != nullptr) {
    Vector Ax0 = multiplicar(A, x);
    addScaledInPlace(r, Ax0, -1.0); // r = b - A*x0
  }

  Vector p = r;

  // Norma inicial do resíduo — usada no critério de parada relativo
  double rsold         = r.linear_product(r);
  const double r0norm  = std::sqrt(rsold);

  // Sistema já satisfaz a tolerância no chute inicial
  if (r0norm <= tolerance) {
    return x;
  }

  // Threshold absoluto baseado na norma inicial (critério relativo)
  const double stopThreshold = tolerance * r0norm;

  // ------------------------------------------------------------------
  // Loop principal do Gradiente Conjugado
  // ------------------------------------------------------------------
  for (int iteration = 0; iteration < maxIterations; ++iteration) {

    Vector Ap = multiplicar(A, p);

    double denominator = p.linear_product(Ap); // pᵀ A p

    // Guarda numérica: A deve ser definida positiva → pᵀAp > 0
    if (std::abs(denominator) < 1e-18) {
      throw std::runtime_error(
          "Falha numerica no gradiente conjugado: p^T A p ~= 0. "
          "Verifique se A e definida positiva.");
    }

    double alpha = rsold / denominator;

    addScaledInPlace(x, p,   alpha);  // x = x + α·p
    addScaledInPlace(r, Ap, -alpha);  // r = r - α·A·p

    double rsnew = r.linear_product(r);

    // Critério de parada relativo: ‖r‖ ≤ tol · ‖r₀‖
    if (std::sqrt(rsnew) <= stopThreshold) {
      return x;
    }

    double beta = rsnew / rsold; // β = ‖r_new‖² / ‖r_old‖²

    // p = r + β·p
    for (int i = 0; i < p.getLength(); ++i) {
      p.setValue(i, r.getValue(i) + beta * p.getValue(i));
    }

    rsold = rsnew;
  }

  throw std::runtime_error(
      "Gradiente conjugado nao convergiu no numero maximo de iteracoes.");
}

} // namespace conjugate_gradient
