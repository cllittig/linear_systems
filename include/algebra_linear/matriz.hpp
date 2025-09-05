#ifndef MATRIZ_HPP
#define MATRIZ_HPP

#include <vector>
#include <stdexcept>
#include <cstdio>

class Matriz{

private:

    //atributos privados
    int rows;
    int columns;
    std::vector<std::vector<double>> data;
    bool isSimetricRec(int i, int j) const;

public:
    //construtor e destrutor
    Matriz();
    Matriz(int rows, int columns);
    Matriz(const Matriz& other);
    
    ~Matriz();

    //métodos de acesso
    int getRows() const;
    int getColumns() const;
    double getValue(int row, int column) const;

    //metodo de modificação
    void setValue(int row, int column, double value);

    //operações com matriciais

    Matriz operator+(const Matriz& other) const;
    Matriz operator-(const Matriz& other) const;
    Matriz operator*(double realNumber) const;
    Matriz operator*(const Matriz& other) const;

    //Metodos especificos algebra linear
    bool isSquare() const;
    Matriz transpose() const;
    double determinant() const; 
    Matriz inverse() const;
    Matriz identity() const;  
    bool isSimetric() const;

    //Métodos auxiliares
    void printMatriz() const;
    void trocarLinhas(int linha1, int linha2);
};



#endif //MATRIZ_HPP 