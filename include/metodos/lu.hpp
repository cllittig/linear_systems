#ifndef LU_HPP
#define LU_HPP

#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"


namespace lu{
    Matriz lu(const Matriz &A, const Matriz &b);
    Vector solve(const Matriz &A, const Vector &b);
}

#endif //LU_HPP
