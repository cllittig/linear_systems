#ifndef GAUSSJORDAN_HPP
#define GAUSSJORDAN_HPP

#include "algebra_linear/matriz.hpp"


namespace gaussjordan{
    Matriz gauss_jordan(const Matriz &A, const Matriz &b);
}

#endif //GAUSSJORDAN_HPP
