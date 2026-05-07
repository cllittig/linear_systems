#ifndef GAUSSSEIDEL_HPP
#define GAUSSSEIDEL_HPP

#include "algebra_linear/matriz.hpp"


namespace gaussseidel{
    Matriz gauss_seidel(const Matriz& A, const Matriz& b);
}

#endif //GAUSSSEIDEL_HPP
