#include "algebra_linear/matriz.hpp"
#include <stdexcept>
#include <cmath>

//construtor padrão
Matriz::Matriz() : rows(0), columns(0){}

//construtor com parametros
Matriz::Matriz(int rows, int columns) : rows(rows), columns(columns){
    if(rows < 0 || columns < 0){
        throw std::invalid_argument("Number of rows and columns must be non-negative");
    }

    data.resize(rows, std::vector<double>(columns, 0.0));
}

//copia do construtor
Matriz::Matriz(const Matriz& other) : rows(other.getRows()), columns(other.getColumns()), data(other.data){}

//destrutor
Matriz::~Matriz(){
    data.clear();
}

//Metodos de acesso

int Matriz::getRows() const{
    return rows;
}

int Matriz::getColumns() const{
    return columns;
}

double Matriz::getValue(int linha, int coluna) const{
    if(linha < 0 || linha >=rows || coluna <0 || coluna >= columns){
        throw std::out_of_range("error in getValue: index out of range");
    }
    return data[linha][coluna];
}

//Metodos de modificação

void Matriz::setValue(int row, int column, double value){
    if(row < 0 || row >=rows || columns <0 || column >= columns){
        throw std::out_of_range("index out of range");
    }
    data[row][column] = value;
}

//operações matriciais


//soma
Matriz Matriz::operator+(const Matriz& other) const{
    if(rows != other.getRows() || columns != other.getColumns() ){
        throw std::invalid_argument("Number of columns and Rows must be different");
    }

    Matriz responseMatriz(rows, columns);

    for(int i=0; i <rows; i++){
        for(int j =0; j< columns; j++){
            responseMatriz.setValue(i,j,data[i][j] + other.data[i][j]);
        }
    }

    return responseMatriz;
}

//subtração
Matriz Matriz::operator-(const Matriz& other) const{
    if(rows != other.getRows() || columns != other.getColumns() ){
        throw std::invalid_argument("Number of columns and Rows must be different");
    }

    Matriz responseMatriz(rows, columns);

    for(int i=0; i <rows; i++){
        for(int j =0; j< columns; j++){
            responseMatriz.setValue(i,j,data[i][j] - other.data[i][j]);
        }
    }

    return responseMatriz;
}

//multiplicar por real
Matriz Matriz::operator*(double realNumber) const{
    Matriz responseMatriz(rows, columns);

    for(int i=0; i <rows; i++){
        for(int j =0; j< columns; j++){
            responseMatriz.setValue(i,j,data[i][j] * realNumber);
        }
    }

    return responseMatriz;
}

//multiplicação de matriz
Matriz Matriz::operator*(const Matriz& other) const {
    if (columns != other.getRows()) {
        throw std::invalid_argument("Number of columns of the first matrix must be equal to the number of rows of the second matrix");
    }
    Matriz responseMatriz(rows, other.getColumns());

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < other.getColumns(); j++) {
            double tempSum = 0;
            for (int k = 0; k < columns; k++) {
                tempSum += data[i][k] * other.data[k][j];
            }
            responseMatriz.setValue(i, j, tempSum);
        }
    }
    return responseMatriz;
}

//Metodos especificos algebra linear

bool Matriz::isSquare() const{
    if(rows == columns){
        return true;
    }
    return false;
}


//transposta
Matriz Matriz::transpose() const{
    Matriz responseMatriz(columns, rows);

    for(int i=0; i<rows; i++){
        for(int j=0; j<columns; j++){
            responseMatriz.setValue(j,i, data[i][j]);
        }
    }

    return responseMatriz;
}

//simetrica
bool Matriz::isSimetric() const{
    if(rows != columns)
        return false;
    return isSimetricRec(0,0);
}
//Função auxiliar para descobrir simetria na matriz
bool Matriz::isSimetricRec(int i, int j)const{
    if(i>=rows)
        return true;
    if(j>=columns)
        return isSimetricRec(i+1,0);
    if(data[i][j] != data[j][i])
        return false;
    return isSimetricRec(i, j+1);
}


//calculo de determinante  ---- ainda falta terminar
double Matriz::determinant()const{
    if( rows == 1 && columns ==1)
        return data[0][0];
    else if(rows == 2 && columns ==2)
        return (data[0][0]*data[1][1] - data[0][1]*data[1][0]);
    else if(rows ==3 && columns ==3)
        return (
            data[0][0]*data[1][1]*data[2][2]
          + data[0][1]*data[1][2]*data[2][0]
          + data[0][2]*data[1][0]*data[2][1]
          - data[0][2]*data[1][1]*data[2][0]
          - data[0][1]*data[1][0]*data[2][2]
          - data[0][0]*data[1][2]*data[2][1]
        );
    else{
        return 0.0;
    }
}

//Matriz identidade
Matriz Matriz::identity()const{
    if(!isSquare())
        throw std::invalid_argument("Matrix must be square");

    Matriz responseMatriz(rows, columns);

    for(int i=0; i<rows; i++){
        for(int j=0; j<columns; j++){
            if(i == j)
                responseMatriz.setValue(i,j, 1);
            else
                responseMatriz.setValue(i,j, 0);
        }
    }

    return responseMatriz;
}

//matriz inversa -- falta terminar
Matriz Matriz::inverse()const{
    if(!isSquare())
        throw std::invalid_argument("Matrix must be square");
    
    Matriz responseMatriz(rows, columns);



    return responseMatriz;

    
}

// print matriz
void Matriz::printMatriz() const{
    for(int i =0; i< rows; i++){
        for(int j=0; j<columns; j++){
            printf("%f |", data[i][j]);
        }
        printf("\n");
    }
}

void Matriz::trocarLinhas(int linha1, int linha2){
    for(int i =0; i<columns; i++){
        double temp = getValue(linha1, i);
        setValue(linha1, i, getValue(linha2, i));
        setValue(linha2, i, temp);
    }
}
