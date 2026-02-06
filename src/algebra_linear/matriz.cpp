#include "algebra_linear/matriz.hpp"
#include <cblas.h>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <tuple>

// Construtor padrão
Matriz::Matriz() : rows(0), columns(0) {}

Matriz::Matriz(std::tuple<int, int> size, std::vector<double> data)
    : data(data) {
  std::tie(rows, columns) = size;
  if (data.size() % rows != 0) {
    throw "WIP1";
  } else if (data.size() % columns != 0) {
    throw "WIP2";
  }
}

// Construtor com parametros
Matriz::Matriz(int rows, int columns) : rows(rows), columns(columns) {
  if (rows < 0 || columns < 0) {
    throw std::invalid_argument(
        "Number of rows and columns must be non-negative");
  }
  data.resize(rows * columns);
}

// Copia do construtor
Matriz::Matriz(const Matriz &other)
    : rows(other.getRows()), columns(other.getColumns()), data(other.data) {}

// Destrutor
Matriz::~Matriz() { data.clear(); }

// Métodos de acesso

int Matriz::getRows() const { return rows; }

int Matriz::getColumns() const { return columns; }

double Matriz::getValue(int row, int column) const {
  if (row < 0 || row >= rows || column < 0 || column >= columns)
    throw std::out_of_range("error in getValue: index out of range");
  return data[row * columns + column];
}

// Métodos de modificação

void Matriz::setValue(int row, int column, double value) {
  if (row < 0 || row >= rows || columns < 0 || column >= columns) {
    throw std::out_of_range("index out of range");
  }
  data[row * columns + column] = value;
}

// Operações matriciais

// Soma
Matriz Matriz::operator+(const Matriz &other) const {
  if (rows != other.getRows() || columns != other.getColumns()) {
    throw std::invalid_argument("Number of columns and Rows must be different");
  }

  Matriz responseMatriz(rows, columns);

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      responseMatriz.setValue(i, j, getValue(i, j) + other.getValue(i, j));
    }
  }

  return responseMatriz;
}

// Subtração
Matriz Matriz::operator-(const Matriz &other) const {
  if (rows != other.getRows() || columns != other.getColumns()) {
    throw std::invalid_argument("Number of columns and Rows must be different");
  }

  Matriz responseMatriz(rows, columns);

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      responseMatriz.setValue(i, j, getValue(i, j) - other.getValue(i, j));
    }
  }

  return responseMatriz;
}

// Multiplicar por real
Matriz Matriz::operator*(double realNumber) const {
  Matriz responseMatriz(rows, columns);

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      responseMatriz.setValue(i, j, getValue(i, j) * realNumber);
    }
  }

  return responseMatriz;
}

// Multiplicação de matriz
Matriz Matriz::operator*(const Matriz &other) const {
  if (columns != other.getRows()) {
    throw std::invalid_argument(
        "Number of columns of the first matrix must be equal to the number of "
        "rows of the second matrix");
  }
  Matriz responseMatriz(rows, other.getColumns());

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < other.getColumns(); j++) {
      double tempSum = 0;
      for (int k = 0; k < columns; k++) {
        tempSum += getValue(i, j) * other.getValue(k, j);
      }
      responseMatriz.setValue(i, j, tempSum);
    }
  }
  return responseMatriz;
}

// Métodos especificos algebra linear

bool Matriz::isSquare() const {
  if (rows == columns) {
    return true;
  }
  return false;
}

// Transposta
Matriz Matriz::transpose() const {
  Matriz responseMatriz(columns, rows);

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      responseMatriz.setValue(j, i, getValue(i, j));
    }
  }

  return responseMatriz;
}

// Simetrica
bool Matriz::isSimetric() const {
  if (rows != columns)
    return false;
  return isSimetricRec(0, 0);
}
// Função auxiliar para descobrir simetria na matriz
bool Matriz::isSimetricRec(int i, int j) const {
  if (i >= rows)
    return true;
  if (j >= columns)
    return isSimetricRec(i + 1, 0);
  if (getValue(i, j) != getValue(j, i))
    return false;
  return isSimetricRec(i, j + 1);
}

// Cálculo de determinante  ---- ainda falta terminar
double Matriz::determinant() const {
  if (rows == 1 && columns == 1)
    return getValue(0, 0);
  else if (rows == 2 && columns == 2)
    return (getValue(0, 0) * getValue(1, 1) - getValue(0, 1) * getValue(1, 0));
  else if (rows == 3 && columns == 3)
    return (getValue(0, 0) * getValue(1, 1) * getValue(2, 2) +
            getValue(0, 1) * getValue(1, 2) * getValue(2, 0) +
            getValue(0, 2) * getValue(1, 0) * getValue(2, 1) -
            getValue(0, 2) * getValue(1, 1) * getValue(2, 0) -
            getValue(0, 1) * getValue(1, 0) * getValue(2, 2) -
            getValue(0, 0) * getValue(1, 2) * getValue(2, 1));
  else {
    return 0.0;
  }
}

// Matriz identidade
Matriz Matriz::identity() const {
  if (!isSquare())
    throw std::invalid_argument("Matrix must be square");

  Matriz responseMatriz(rows, columns);

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      if (i == j)
        responseMatriz.setValue(i, j, 1);
      else
        responseMatriz.setValue(i, j, 0);
    }
  }

  return responseMatriz;
}

// Matriz inversa -- falta terminar
Matriz Matriz::inverse() const {
  if (!isSquare())
    throw std::invalid_argument("Matrix must be square");

  Matriz responseMatriz(rows, columns);

  return responseMatriz;
}

// Print matriz
void Matriz::printMatriz() const {
  if (rows == 0 && columns == 0) {
    printf("| empty |\n");
    return;
  }
  for (int i = 0; i < rows; i++) {
    printf("| ");
    for (int j = 0; j < columns; j++) {
      printf("%9f ", getValue(i, j));
    }
    printf("|\n");
  }
}

void Matriz::trocarLinhas(int linha1, int linha2) {
  for (int i = 0; i < columns; i++) {
    double temp = getValue(linha1, i);
    setValue(linha1, i, getValue(linha2, i));
    setValue(linha2, i, temp);
  }
}
