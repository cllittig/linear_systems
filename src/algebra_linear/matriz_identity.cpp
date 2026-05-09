// Arquivo separado para evitar conflitos com macros da biblioteca FLAME
#include "algebra_linear/matriz.hpp"

Matriz Matriz::identity(int rows, int cols) {
  Matriz I(rows, cols);
  for (int idx = 0; idx < rows && idx < cols; idx++) {
    I.setValue(idx, idx, 1.0);
  }
  return I;
}
