#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <vector>
#include <stdexcept>
#include <cstdio>

class Vector{

private:    
    int length;
    std::vector<double> data;

public:

    //construtor e destrutor
    Vector();
    Vector(int len);
    Vector(const Vector& other);

    ~Vector();

    //Metodos de acesso
    int getLength() const;
    double getValue(int position) const;

    //metodo de Modificação
    void setValue(int position, double value);

    //operações com vetores
    bool equality(const Vector& other) const;
    Vector assigment(const Vector& other) const;

};


#endif //VECTOR_HPP