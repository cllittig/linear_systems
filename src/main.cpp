// bibliotecas padrão 
#include <iostream>
#include <vector>
#include <cstdio>
#include <complex> 

// bibliotecas auxiliares externas
#include <FLAME.h>

//classes referentes a algebra linear	
#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"

//metodos numericos
#include "metodos/cholesky.hpp"
#include "metodos/conjugate_gradient.hpp"

// retirando a definição de min e max, evitando conflito com a lib flame
#undef min
#undef max

//função main
int main() {
    try {
        // Entrada do sistema linear Ax = b
        Matriz A(3, 3);
        A.setValue(0, 0, 4);   A.setValue(0, 1, 12);  A.setValue(0, 2, -16);
        A.setValue(1, 0, 12);  A.setValue(1, 1, 37);  A.setValue(1, 2, -43);
        A.setValue(2, 0, -16); A.setValue(2, 1, -43); A.setValue(2, 2, 98);

        Vector b(3);
        b.setValue(0, 0);
        b.setValue(1, 6);
        b.setValue(2, 39);

        const double tolerance = 1e-10;
        const int maxIterations = 100;

        std::cout << "--- Entradas do Gradiente Conjugado ---" << std::endl;
        std::cout << "Matriz A:" << std::endl;
        A.print();
        std::cout << "Vetor b:" << std::endl;
        for (int i = 0; i < b.getLength(); ++i) {
            std::cout << "b[" << i << "] = " << b.getValue(i) << std::endl;
        }
        std::cout << "tolerance = " << tolerance << std::endl;
        std::cout << "maxIterations = " << maxIterations << std::endl;

        std::cout << "\n--- Resolvendo Ax = b com Gradiente Conjugado ---" << std::endl;
        Vector x = conjugate_gradient::solve(A, b, tolerance, maxIterations);

        std::cout << "\n--- Saida (vetor x) ---" << std::endl;
        for (int i = 0; i < x.getLength(); i++) {
            std::cout << "x[" << i << "] = " << x.getValue(i) << std::endl;
        }

        std::cout << "\n--- Verificacao rapida (A*x) ---" << std::endl;
        Vector Ax = multiplicar(A, x);
        for (int i = 0; i < Ax.getLength(); ++i) {
            std::cout << "A*x[" << i << "] = " << Ax.getValue(i)
                      << " | b[" << i << "] = " << b.getValue(i) << std::endl;
        }

    } catch (const std::exception& e) {
        std::cerr << "Erro: " << e.what() << std::endl;
    }

    return 0;
}
