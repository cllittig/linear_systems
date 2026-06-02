#ifndef MATRIZ_HPP
#define MATRIZ_HPP

#include <stdexcept>
#include <vector>

class Matriz{

private:

    //atributos privados
    int rows;
    int columns;
    std::vector<double> data;

public:
    //construtor e destrutor
    Matriz();
    Matriz(int rows, int columns);
    Matriz(std::tuple<int,int>, std::vector<double>);
    // Cópia, atribuição, move e destrutor seguem a regra do zero:
    // std::vector já gerencia os recursos corretamente, então os membros
    // especiais gerados pelo compilador são suficientes (e corretos).

    //métodos de acesso
    int getRows() const;
    int getColumns() const;

    inline double getValue(int row, int column) const {
#ifndef NDEBUG
        if (row < 0 || row >= rows || column < 0 || column >= columns)
            throw std::out_of_range("error in getValue: index out of range");
#endif
        return data[row * columns + column];
    }

    inline void setValue(int row, int column, double value) {
#ifndef NDEBUG
        if (row < 0 || row >= rows || column < 0 || column >= columns)
            throw std::out_of_range("index out of range");
#endif
        data[row * columns + column] = value;
    }

    inline double* rawData() { return data.data(); }
    inline const double* rawData() const { return data.data(); }

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
    void print() const;
    void trocarLinhas(int linha1, int linha2);
    void linhaAxpy(int dest, int src, double alpha, int col_inicio = 0);
    Matriz clonar() const;

    static Matriz identity(int n, int m);

};



#endif //MATRIZ_HPP 
