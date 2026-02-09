
#include "metodos/gauss_jordan.hpp"
#include <cstdio>

int main() {

    Matriz a({3, 3}, { 
        3, 0,  2, 
        2, 0, -2,
        0, 1,  1,
    });

// 0.2 0.2 0 
// -0.2 0.3 1 
// 0.2 -0.3 0 

    a.print();

    auto resultado = a.inverse();

    printf("\n");

    resultado.print();
    return 0;
}
