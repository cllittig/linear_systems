
#include "metodos/lu.hpp"
#include <cstdio>

int main() {

    Matriz a({2, 2}, { 
        3, 1,
        -6, -4,
    });

    Matriz b({1, 2}, { 5, 11, });


    a.print();

    b.print();

    auto resultado = lu::lu(a, b);

    printf("\n");

    resultado.print();

    return 0;
}
