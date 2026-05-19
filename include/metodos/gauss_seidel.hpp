#ifndef GAUSSSEIDEL_HPP
#define GAUSSSEIDEL_HPP

#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"

namespace gaussseidel {
    Matriz gauss_seidel(const Matriz& A, const Matriz& b);
    Vector solve(const Matriz& A, const Vector& b);
}

#endif //GAUSSSEIDEL_HPP
