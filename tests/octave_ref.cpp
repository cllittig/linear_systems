/**
 * @file octave_ref.cpp
 * @brief Soluções de referência via API C++ do GNU Octave (liboctave).
 *
 * Este arquivo é o ÚNICO lugar do projeto que inclui os cabeçalhos do Octave,
 * isolando o conflito de nomes entre a `template<class T> class lu;` global do
 * Octave e o `namespace lu` da biblioteca (ver octave_ref.hpp).
 *
 * Toda a interação com o Octave é feita no nível da API C++ (Matrix /
 * ColumnVector / Matrix::solve, que internamente chama LAPACK). Não há nenhuma
 * chamada a system(), popen(), exec() nem execução de scripts .m.
 */

#include "octave_ref.hpp"

#include <octave/oct.h>

namespace {

// Matriz (biblioteca) -> Matrix (Octave)
Matrix to_octave(const Matriz& A) {
    Matrix result(A.getRows(), A.getColumns());
    for (int i = 0; i < A.getRows(); ++i)
        for (int j = 0; j < A.getColumns(); ++j)
            result(i, j) = A.getValue(i, j);
    return result;
}

// Vector (biblioteca) -> ColumnVector (Octave)
ColumnVector to_octave(const Vector& b) {
    ColumnVector result(b.getLength());
    for (int i = 0; i < b.getLength(); ++i)
        result(i) = b.getValue(i);
    return result;
}

// ColumnVector (Octave) -> Vector (biblioteca)
Vector from_octave(const ColumnVector& v) {
    Vector result(v.numel());
    for (octave_idx_type i = 0; i < v.numel(); ++i)
        result.setValue(i, v(i));
    return result;
}

// Núcleo comum: resolve A x = b com Matrix::solve (mldivide do Octave).
Vector solve_reference(const Matriz& A, const Vector& b) {
    Matrix oct_A = to_octave(A);
    ColumnVector oct_b = to_octave(b);
    ColumnVector oct_x = oct_A.solve(oct_b);
    return from_octave(oct_x);
}

} // namespace

namespace octave_ref {

// Matrix::solve usa LU com pivoteamento parcial para matrizes gerais.
Vector solve_lu(const Matriz& A, const Vector& b) {
    return solve_reference(A, b);
}

// Para matrizes SPD, Matrix::solve detecta a estrutura e usa Cholesky.
Vector solve_cholesky(const Matriz& A, const Vector& b) {
    return solve_reference(A, b);
}

// O Octave não expõe Gauss-Jordan puro na API C++; o A\b (eliminação direta)
// é a referência usada para comparação.
Vector solve_gauss_jordan(const Matriz& A, const Vector& b) {
    return solve_reference(A, b);
}

// O Octave não expõe Gauss-Seidel na API C++; usamos A\b como referência
// numérica para validar convergência e precisão da implementação iterativa.
Vector solve_gauss_seidel(const Matriz& A, const Vector& b) {
    return solve_reference(A, b);
}

// pcg do Octave não é exposto na API C++ sem o interpretador; para SPD bem
// condicionada, a solução direta serve de referência exata para o CG.
Vector solve_pcg(const Matriz& A, const Vector& b) {
    return solve_reference(A, b);
}

} // namespace octave_ref
