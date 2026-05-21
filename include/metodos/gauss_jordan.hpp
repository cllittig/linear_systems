#ifndef GAUSSJORDAN_HPP
#define GAUSSJORDAN_HPP

#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"


namespace gaussjordan{
    Matriz gauss_jordan(const Matriz &A, const Matriz &b);
    Vector solve(const Matriz &A, const Vector &b);
}

#endif //GAUSSJORDAN_HPP
