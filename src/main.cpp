#include "algebra_linear/matriz.hpp"
#include "metodos/gauss_seidel.hpp"

#undef min
#undef max

int main() {
    // 10x - y + 2z = 6
    //  -x + 11y - z = 25
    //  2x -  y + 10z = -11

    // esperado: x=1, y=2, z=-1

    Matriz A(3, 3);
    double dadosA[] = { 10, -1,  2,
                        -1, 11, -1,
                         2, -1, 10 };
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            A.setValue(i, j, dadosA[i * 3 + j]);

    Matriz b(3, 1);
    b.setValue(0, 0,   6);
    b.setValue(1, 0,  25);
    b.setValue(2, 0, -11);

    Matriz x = gaussseidel::gauss_seidel(A, b);
    x.print(); // [1, 2, -1]
    // Precisamos de testes urgentemente!!!!!!
}
