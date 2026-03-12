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

// retirando a definição de min e max, evitando conflito com a lib flame
#undef min
#undef max

//função main
int main() {
	try {
        // 1. Criar uma matriz A (3x3) simétrica e definida positiva
        // Exemplo: 
        // [ 4  12 -16 ]
        // [ 12 37 -43 ]
        // [-16 -43 98 ]
        Matriz A(3, 3);
        A.setValue(0, 0, 4);   A.setValue(0, 1, 12);  A.setValue(0, 2, -16);
        A.setValue(1, 0, 12);  A.setValue(1, 1, 37);  A.setValue(1, 2, -43);
        A.setValue(2, 0, -16); A.setValue(2, 1, -43); A.setValue(2, 2, 98);

        // 2. Criar um vetor b para o sistema Ax = b
        // Se quisermos que a resposta (x) seja [1, 1, 1], o b deve ser:
        Vector b(3);
        b.setValue(0, 0);   // 4(1) + 12(1) - 16(1) = 0
        b.setValue(1, 6);   // 12(1) + 37(1) - 43(1) = 6
        b.setValue(2, 39);  // -16(1) - 43(1) + 98(1) = 39

        std::cout << "--- Matriz A ---" << std::endl;
        A.print(); // Supondo que sua classe matriz tenha o método print

        std::cout << "\n--- Resolvendo Ax = b ---" << std::endl;
        
        // 3. Chamar sua função solve
        Vector x = cholesky::solve(A, b);

        std::cout << "\n--- Resultado do Vetor x ---" << std::endl;
        for (int i = 0; i < x.getLength(); i++) {
            std::cout << "x[" << i << "] = " << x.getValue(i) << std::endl;
        }

    } catch (const std::exception& e) {
        std::cerr << "Erro: " << e.what() << std::endl;
    }

    return 0;
}
