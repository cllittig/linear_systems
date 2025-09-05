#include "algebra_linear/vetor.hpp"
#include <stdexcept>
#include <vector>
#include <cmath>

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

//metodos de acesso
int Vector::getLength()const{
    return length;
}

double Vector::getValue(int position)const{
    return data[position];
}

//metodo de modificação
void Vector::setValue(int position, double value){
    data[position] = value;
}

//operações com vetores
bool Vector::equality(const Vector& other)const{

    if(other.getLength() != length){
        return false;
    }

    for(int i =0; i<length; i++){
        if(data[i] != other.getValue(i))
            return false;
    }

    return true;
}

Vector Vector::operator+(const Vector& other) const{
    
    Vector responseVector(length);

    for(int i=0; i<length; i++){
        double sum = other.getValue(i) + getValue(i);
        responseVector.setValue(i, sum);
    }

    return responseVector;

}


Vector Vector::operator*(const int escale) const{
    
    Vector responseVector(length);

    for(int i=0; i<length; i++){
        responseVector.setValue(i, getValue(i)*escale);
    }

    return responseVector;
}


Vector Vector::operator-(const Vector& other) const{
    
    Vector responseVector(length);

    for(int i=0; i<length; i++){
        double sub =getValue(i) - other.getValue(i);
        responseVector.setValue(i, sub);
    }

    return responseVector;

}

Vector Vector::operator*(const Vector& other) const{
    Vector responseVector(length);

    for(int i=0; i<length; i++){
        double mult =getValue(i) * other.getValue(i);
        responseVector.setValue(i, mult);
    }

    return responseVector;
}

Vector Vector::axpy(const Vector& y, double alpha) const{

    return (*this) * alpha + y;

}

double Vector::euclidian_length() const{
    double response =0;

    for(double a: data){
        response += a*a;
    }

    return std::sqrt(response);
}