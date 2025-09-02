#include "algebra_linear/vetor.hpp"
#include <stdexcept>
#include <vector>

//construtor padrão
Vector::Vector():length(0){}

//construtor parametros
Vector::Vector(int len) : length(len){
    if(len<0)
        throw std::invalid_argument("this length of vector must be non-negative");
    
    data.resize(len);
}

//construtor copia
Vector::Vector(const Vector& other): length(other.getLength()), data(other.data){};

//destrutor
Vector::~Vector(){
    data.clear();
}