#ifndef CONJUGATE_GRADIENT_HPP
#define CONJUGATE_GRADIENT_HPP

#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"

namespace conjugate_gradient {

Vector solve(const Matriz &A, const Vector &b, double tolerance = 1e-10,
             int maxIterations = 1000);

} // namespace conjugate_gradient

#endif // CONJUGATE_GRADIENT_HPP
