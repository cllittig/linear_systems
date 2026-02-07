#ifndef GAUSSJORDAN_HPP
#define GAUSSJORDAN_HPP

#include <stdexcept>
#include "algebra_linear/matriz.hpp"


namespace gaussjordan{
    Matriz gauss_jordan(Matriz data, int current_row =0);
}

#endif //GAUSSJORDAN_HPP
