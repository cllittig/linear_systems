#include <cstdio>
#include "../include/metodos/gauss_jordan.hpp"

int main() {

    // Exemplo de uso:// ...existing code...
    Matriz A(3, 4);
    A.setValue(0, 0, 2);  A.setValue(0, 1, 3);  A.setValue(0, 2, -1);  A.setValue(0, 3, -7);
    A.setValue(1, 0, 1);  A.setValue(1, 1, 1);  A.setValue(1, 2, 1);  A.setValue(1, 3, 4);
    A.setValue(2, 0, -1);  A.setValue(2, 1, -2);  A.setValue(2, 2, 3); A.setValue(2, 3, 15);
    // ...existing code...  
    std::vector<double> resultado = gaussjordan::gauss_jordan(A);

    for(double a : resultado)
        printf("resultado: %f", a);
    
    
    return 0;
}