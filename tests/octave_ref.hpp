#ifndef OCTAVE_REF_HPP
#define OCTAVE_REF_HPP

#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"

/**
 * @file octave_ref.hpp
 * @brief Interface das soluções de referência calculadas via liboctave.
 *
 * Esta interface NÃO inclui nenhum cabeçalho do Octave e expõe apenas os tipos
 * da própria biblioteca (Matriz/Vector). Isso é proposital: os cabeçalhos do
 * Octave declaram, no escopo global, `template<class T> class lu;`, que colide
 * com o `namespace lu` da biblioteca. Mantendo o Octave confinado em
 * octave_ref.cpp, o arquivo de teste pode incluir os métodos da biblioteca
 * (lu::, cholesky::, ...) sem conflito de nomes.
 *
 * Todas as funções resolvem A x = b com a API numérica do Octave
 * (Matrix::solve), que seleciona a fatoração apropriada (LU com pivoteamento
 * parcial para matrizes gerais; Cholesky para SPD) e chama LAPACK diretamente —
 * sem system()/popen()/exec() e sem script Octave.
 */

namespace octave_ref {

Vector solve_lu(const Matriz& A, const Vector& b);
Vector solve_cholesky(const Matriz& A, const Vector& b);
Vector solve_gauss_jordan(const Matriz& A, const Vector& b);
Vector solve_pcg(const Matriz& A, const Vector& b);

} // namespace octave_ref

#endif // OCTAVE_REF_HPP
