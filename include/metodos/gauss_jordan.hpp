#ifndef GAUSSJORDAN_HPP
#define GAUSSJORDAN_HPP

#include <vector>
#include <stdexcept>
#include "algebra_linear/matriz.hpp"


namespace gaussjordan{
    std::vector<double> gauss_jordan(Matriz data, int current_row =0);
}

#endif //GAUSSJORDAN_HPP