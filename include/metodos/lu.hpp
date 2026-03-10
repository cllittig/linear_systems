#ifndef LU_HPP
#define LU_HPP

#include "algebra_linear/matriz.hpp"


namespace lu{
    Matriz lu(const Matriz &A, const Matriz &b);
}

#endif //LU_HPP
