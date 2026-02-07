
#include <cstdio>
#include "metodos/gauss_jordan.hpp"

int main() {

    Matriz a({3, 3}, { 
        3, 0,  2, 
        2, 0, -2,
        0, 1,  1,
    });

    a.printMatriz();

    auto resultado = gaussjordan::gauss_jordan(a);

    resultado.printMatriz();
    return 0;
}
