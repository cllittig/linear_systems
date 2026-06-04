#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <vector>
#include <stdexcept>
#include <cstdio>
#include "algebra_linear/matriz.hpp"

class Vector{

private:    
    int length;
    std::vector<double> data;

public:

    //construtor e destrutor
    Vector();
    Vector(int len);
    // Cópia, atribuição, move e destrutor seguem a regra do zero:
    // std::vector gerencia o armazenamento, então os membros especiais
    // gerados pelo compilador são suficientes e corretos.

    //Metodos de acesso
    int getLength() const;
    double getValue(int position) const;

    //metodo de Modificação
    void setValue(int position, double value);

    //operações com vetores
    bool equality(const Vector& other) const;
    Vector operator+(const Vector& other) const;
    Vector operator*(double scale) const;
    Vector operator*(const Vector& other) const;
    Vector operator-(const Vector& other) const;
    Vector axpy(const Vector& y, double alpha) const;
    double euclidian_length() const;
    double linear_product(const Vector& other) const;
};

Vector multiplicar(const Matriz& A, const Vector& v);


#endif //VECTOR_HPP