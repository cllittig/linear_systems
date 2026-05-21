#ifndef CONJUGATE_GRADIENT_HPP
#define CONJUGATE_GRADIENT_HPP

#include <linear_systems/algebra_linear/matriz.hpp>
#include <linear_systems/algebra_linear/vector.hpp>
#include <linear_systems/metodos/solver_utils.hpp>

namespace conjugate_gradient {

Vector solve(const Matriz &A, const Vector &b, double tolerance = 1e-10,
             int maxIterations = 1000, const Vector* x0 = nullptr,
             solver::IterationInfo* info = nullptr, bool useJacobiPrecond = false);

// Solve for multiple RHS (each column of B is a RHS). Returns matrix X with
// solutions in corresponding columns.
Matriz solve(const Matriz &A, const Matriz &B, double tolerance = 1e-10,
             int maxIterations = 1000, solver::IterationInfo* info = nullptr,
             bool useJacobiPrecond = false);

} // namespace conjugate_gradient

#endif // CONJUGATE_GRADIENT_HPP
