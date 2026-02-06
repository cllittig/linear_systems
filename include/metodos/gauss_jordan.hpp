#ifndef GAUSSJORDAN_HPP
#define GAUSSJORDAN_HPP

#include <vector>
#include "matriz.hpp"


namespace gaussjordan{
    Matriz gauss_jordan(Matriz data, int current_row =0);
}

#endif //GAUSSJORDAN_HPP
