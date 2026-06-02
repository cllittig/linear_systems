#include <gtest/gtest.h>
#include <stdexcept>
#include <tuple>
#include <vector>
#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"

/**
 * @file test_matriz.cpp
 * @brief Testes unitários da classe base Matriz (GoogleTest).
 *
 * Cobre construção, acesso, operadores aritméticos, operações de álgebra
 * linear (transposta, simetria, determinante, inversa, identidade) e os
 * utilitários usados pelos solvers (trocarLinhas, linhaAxpy, clonar,
 * multiplicar matriz-vetor), incluindo os casos de erro.
 *
 * Observação: os bound checks de getValue/setValue só existem quando NDEBUG
 * NÃO está definido (build de debug). Os testes correspondentes são protegidos
 * por #ifndef NDEBUG para não falharem no build release (padrão).
 */

static const double TOL = 1e-12;

namespace {

Matriz make2x2(double a, double b, double c, double d) {
    Matriz m(2, 2);
    m.setValue(0, 0, a); m.setValue(0, 1, b);
    m.setValue(1, 0, c); m.setValue(1, 1, d);
    return m;
}

} // namespace

// --------------------------------------------------------------------------
// Construção e acesso
// --------------------------------------------------------------------------

TEST(MatrizTest, DefaultConstructorIsEmpty) {
    Matriz m;
    EXPECT_EQ(m.getRows(), 0);
    EXPECT_EQ(m.getColumns(), 0);
}

TEST(MatrizTest, DimensionConstructor) {
    Matriz m(3, 4);
    EXPECT_EQ(m.getRows(), 3);
    EXPECT_EQ(m.getColumns(), 4);
    // Inicializada com zeros
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 4; ++j)
            EXPECT_NEAR(m.getValue(i, j), 0.0, TOL);
}

TEST(MatrizTest, TupleDataConstructor) {
    Matriz m(std::make_tuple(2, 2), std::vector<double>{1.0, 2.0, 3.0, 4.0});
    EXPECT_EQ(m.getRows(), 2);
    EXPECT_EQ(m.getColumns(), 2);
    EXPECT_NEAR(m.getValue(0, 0), 1.0, TOL);
    EXPECT_NEAR(m.getValue(1, 1), 4.0, TOL);
}

TEST(MatrizTest, TupleDataConstructorSizeMismatchThrows) {
    EXPECT_THROW(
        Matriz(std::make_tuple(2, 2), std::vector<double>{1.0, 2.0, 3.0}),
        std::invalid_argument);
}

TEST(MatrizTest, TupleDataConstructorNonPositiveDimsThrows) {
    EXPECT_THROW(
        Matriz(std::make_tuple(0, 2), std::vector<double>{}),
        std::invalid_argument);
}

TEST(MatrizTest, GetSetValue) {
    Matriz m(2, 2);
    m.setValue(0, 1, 5.5);
    EXPECT_NEAR(m.getValue(0, 1), 5.5, TOL);
}

#ifndef NDEBUG
TEST(MatrizTest, GetValueOutOfRangeThrows) {
    Matriz m(2, 2);
    EXPECT_THROW(m.getValue(2, 0), std::out_of_range);
    EXPECT_THROW(m.getValue(0, -1), std::out_of_range);
}

TEST(MatrizTest, SetValueOutOfRangeThrows) {
    Matriz m(2, 2);
    EXPECT_THROW(m.setValue(5, 5, 1.0), std::out_of_range);
}
#endif

TEST(MatrizTest, CopyConstructorIsIndependent) {
    Matriz original = make2x2(1, 2, 3, 4);
    Matriz copy = original;            // copy ctor (regra do zero)
    copy.setValue(0, 0, 99.0);
    EXPECT_NEAR(original.getValue(0, 0), 1.0, TOL);  // original intacto
    EXPECT_NEAR(copy.getValue(0, 0), 99.0, TOL);
}

TEST(MatrizTest, CopyAssignmentIsIndependent) {
    Matriz original = make2x2(1, 2, 3, 4);
    Matriz other(2, 2);
    other = original;                  // copy assignment (regra do zero)
    other.setValue(1, 1, -7.0);
    EXPECT_NEAR(original.getValue(1, 1), 4.0, TOL);
    EXPECT_NEAR(other.getValue(1, 1), -7.0, TOL);
}

// --------------------------------------------------------------------------
// Operadores aritméticos
// --------------------------------------------------------------------------

TEST(MatrizTest, Addition) {
    Matriz r = make2x2(1, 2, 3, 4) + make2x2(5, 6, 7, 8);
    EXPECT_NEAR(r.getValue(0, 0), 6.0, TOL);
    EXPECT_NEAR(r.getValue(0, 1), 8.0, TOL);
    EXPECT_NEAR(r.getValue(1, 0), 10.0, TOL);
    EXPECT_NEAR(r.getValue(1, 1), 12.0, TOL);
}

TEST(MatrizTest, AdditionIncompatibleThrows) {
    EXPECT_THROW(make2x2(1, 2, 3, 4) + Matriz(3, 3), std::invalid_argument);
}

TEST(MatrizTest, Subtraction) {
    Matriz r = make2x2(5, 6, 7, 8) - make2x2(1, 2, 3, 4);
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            EXPECT_NEAR(r.getValue(i, j), 4.0, TOL);
}

TEST(MatrizTest, SubtractionIncompatibleThrows) {
    EXPECT_THROW(make2x2(1, 2, 3, 4) - Matriz(3, 3), std::invalid_argument);
}

TEST(MatrizTest, ScalarMultiplication) {
    Matriz r = make2x2(1, 2, 3, 4) * 2.0;
    EXPECT_NEAR(r.getValue(0, 0), 2.0, TOL);
    EXPECT_NEAR(r.getValue(1, 1), 8.0, TOL);
}

TEST(MatrizTest, MatrixMultiplication) {
    // [[1,2],[3,4]] * [[5,6],[7,8]] = [[19,22],[43,50]]
    Matriz r = make2x2(1, 2, 3, 4) * make2x2(5, 6, 7, 8);
    EXPECT_NEAR(r.getValue(0, 0), 19.0, TOL);
    EXPECT_NEAR(r.getValue(0, 1), 22.0, TOL);
    EXPECT_NEAR(r.getValue(1, 0), 43.0, TOL);
    EXPECT_NEAR(r.getValue(1, 1), 50.0, TOL);
}

TEST(MatrizTest, MatrixMultiplicationNonSquare) {
    // (2x3) * (3x2) = (2x2)
    Matriz A(2, 3);
    A.setValue(0, 0, 1); A.setValue(0, 1, 2); A.setValue(0, 2, 3);
    A.setValue(1, 0, 4); A.setValue(1, 1, 5); A.setValue(1, 2, 6);
    Matriz B(3, 2);
    B.setValue(0, 0, 7);  B.setValue(0, 1, 8);
    B.setValue(1, 0, 9);  B.setValue(1, 1, 10);
    B.setValue(2, 0, 11); B.setValue(2, 1, 12);
    Matriz C = A * B;
    EXPECT_EQ(C.getRows(), 2);
    EXPECT_EQ(C.getColumns(), 2);
    EXPECT_NEAR(C.getValue(0, 0), 58.0, TOL);   // 7+18+33
    EXPECT_NEAR(C.getValue(0, 1), 64.0, TOL);   // 8+20+36
    EXPECT_NEAR(C.getValue(1, 0), 139.0, TOL);  // 28+45+66
    EXPECT_NEAR(C.getValue(1, 1), 154.0, TOL);  // 32+50+72
}

TEST(MatrizTest, MatrixMultiplicationIncompatibleThrows) {
    Matriz A(2, 3);
    Matriz B(2, 2);  // colunas de A (3) != linhas de B (2)
    EXPECT_THROW(A * B, std::invalid_argument);
}

// --------------------------------------------------------------------------
// Operações de álgebra linear
// --------------------------------------------------------------------------

TEST(MatrizTest, IsSquare) {
    EXPECT_TRUE(Matriz(3, 3).isSquare());
    EXPECT_FALSE(Matriz(2, 3).isSquare());
}

TEST(MatrizTest, Transpose) {
    Matriz m(2, 3);
    m.setValue(0, 0, 1); m.setValue(0, 1, 2); m.setValue(0, 2, 3);
    m.setValue(1, 0, 4); m.setValue(1, 1, 5); m.setValue(1, 2, 6);
    Matriz t = m.transpose();
    EXPECT_EQ(t.getRows(), 3);
    EXPECT_EQ(t.getColumns(), 2);
    EXPECT_NEAR(t.getValue(0, 0), 1.0, TOL);
    EXPECT_NEAR(t.getValue(2, 1), 6.0, TOL);
    EXPECT_NEAR(t.getValue(1, 0), 2.0, TOL);
}

TEST(MatrizTest, IsSimetric) {
    EXPECT_TRUE(make2x2(1, 2, 2, 1).isSimetric());
    EXPECT_FALSE(make2x2(1, 2, 3, 4).isSimetric());
    EXPECT_FALSE(Matriz(2, 3).isSimetric());  // não quadrada
}

TEST(MatrizTest, Determinant1x1) {
    Matriz m(1, 1);
    m.setValue(0, 0, 7.0);
    EXPECT_NEAR(m.determinant(), 7.0, TOL);
}

TEST(MatrizTest, Determinant2x2) {
    // det([[1,2],[3,4]]) = -2
    EXPECT_NEAR(make2x2(1, 2, 3, 4).determinant(), -2.0, TOL);
}

TEST(MatrizTest, Determinant3x3) {
    // det([[1,2,3],[0,1,4],[5,6,0]]) = 1
    Matriz m(3, 3);
    m.setValue(0, 0, 1); m.setValue(0, 1, 2); m.setValue(0, 2, 3);
    m.setValue(1, 0, 0); m.setValue(1, 1, 1); m.setValue(1, 2, 4);
    m.setValue(2, 0, 5); m.setValue(2, 1, 6); m.setValue(2, 2, 0);
    EXPECT_NEAR(m.determinant(), 1.0, TOL);
}

TEST(MatrizTest, DeterminantNonSquareThrows) {
    EXPECT_THROW(Matriz(2, 3).determinant(), std::invalid_argument);
}

TEST(MatrizTest, DeterminantAboveSupportedSizeThrows) {
    // n > 3 não é suportado: deve lançar, não retornar 0 silenciosamente.
    EXPECT_THROW(Matriz(4, 4).determinant(), std::runtime_error);
}

TEST(MatrizTest, IdentityMember) {
    Matriz I = make2x2(9, 9, 9, 9).identity();
    EXPECT_NEAR(I.getValue(0, 0), 1.0, TOL);
    EXPECT_NEAR(I.getValue(0, 1), 0.0, TOL);
    EXPECT_NEAR(I.getValue(1, 1), 1.0, TOL);
}

TEST(MatrizTest, IdentityMemberNonSquareThrows) {
    EXPECT_THROW(Matriz(2, 3).identity(), std::invalid_argument);
}

TEST(MatrizTest, IdentityStatic) {
    Matriz I = Matriz::identity(3, 3);
    EXPECT_EQ(I.getRows(), 3);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(I.getValue(i, j), (i == j) ? 1.0 : 0.0, TOL);
}

TEST(MatrizTest, InverseTimesOriginalIsIdentity) {
    // A = [[1,2,3],[0,1,4],[5,6,0]] (det = 1, inversível)
    Matriz A(3, 3);
    A.setValue(0, 0, 1); A.setValue(0, 1, 2); A.setValue(0, 2, 3);
    A.setValue(1, 0, 0); A.setValue(1, 1, 1); A.setValue(1, 2, 4);
    A.setValue(2, 0, 5); A.setValue(2, 1, 6); A.setValue(2, 2, 0);
    Matriz inv = A.inverse();
    Matriz prod = A * inv;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(prod.getValue(i, j), (i == j) ? 1.0 : 0.0, 1e-9);
}

TEST(MatrizTest, InverseSingularThrows) {
    // Linhas linearmente dependentes => não inversível
    EXPECT_THROW(make2x2(1, 2, 2, 4).inverse(), std::runtime_error);
}

// --------------------------------------------------------------------------
// Utilitários usados pelos solvers
// --------------------------------------------------------------------------

TEST(MatrizTest, TrocarLinhas) {
    Matriz m(3, 2);
    m.setValue(0, 0, 1); m.setValue(0, 1, 10);
    m.setValue(1, 0, 2); m.setValue(1, 1, 20);
    m.setValue(2, 0, 3); m.setValue(2, 1, 30);
    m.trocarLinhas(0, 2);
    EXPECT_NEAR(m.getValue(0, 0), 3.0, TOL);
    EXPECT_NEAR(m.getValue(0, 1), 30.0, TOL);
    EXPECT_NEAR(m.getValue(2, 0), 1.0, TOL);
    EXPECT_NEAR(m.getValue(2, 1), 10.0, TOL);
}

TEST(MatrizTest, LinhaAxpy) {
    // linha dest += alpha * linha src
    Matriz m(2, 3);
    m.setValue(0, 0, 1); m.setValue(0, 1, 1); m.setValue(0, 2, 1);  // src
    m.setValue(1, 0, 0); m.setValue(1, 1, 0); m.setValue(1, 2, 0);  // dest
    m.linhaAxpy(1, 0, 2.0);  // dest(1) += 2 * src(0)
    EXPECT_NEAR(m.getValue(1, 0), 2.0, TOL);
    EXPECT_NEAR(m.getValue(1, 1), 2.0, TOL);
    EXPECT_NEAR(m.getValue(1, 2), 2.0, TOL);
}

TEST(MatrizTest, LinhaAxpyWithColumnOffset) {
    Matriz m(2, 3);
    m.setValue(0, 0, 5); m.setValue(0, 1, 1); m.setValue(0, 2, 1);  // src
    m.setValue(1, 0, 9); m.setValue(1, 1, 0); m.setValue(1, 2, 0);  // dest
    m.linhaAxpy(1, 0, 1.0, 1);  // só a partir da coluna 1
    EXPECT_NEAR(m.getValue(1, 0), 9.0, TOL);  // coluna 0 inalterada
    EXPECT_NEAR(m.getValue(1, 1), 1.0, TOL);
    EXPECT_NEAR(m.getValue(1, 2), 1.0, TOL);
}

TEST(MatrizTest, Clonar) {
    Matriz original = make2x2(1, 2, 3, 4);
    Matriz clone = original.clonar();
    clone.setValue(0, 0, 42.0);
    EXPECT_NEAR(original.getValue(0, 0), 1.0, TOL);  // cópia profunda
    EXPECT_NEAR(clone.getValue(0, 0), 42.0, TOL);
}

TEST(MatrizTest, MultiplicarMatrizVetor) {
    // [[1,2],[3,4]] * [1,1]^T = [3,7]^T
    Matriz A = make2x2(1, 2, 3, 4);
    Vector v(2);
    v.setValue(0, 1.0);
    v.setValue(1, 1.0);
    Vector r = multiplicar(A, v);
    EXPECT_EQ(r.getLength(), 2);
    EXPECT_NEAR(r.getValue(0), 3.0, TOL);
    EXPECT_NEAR(r.getValue(1), 7.0, TOL);
}

TEST(MatrizTest, MultiplicarDimensionMismatchThrows) {
    Matriz A(2, 3);
    Vector v(2);  // comprimento (2) != colunas de A (3)
    EXPECT_THROW(multiplicar(A, v), std::invalid_argument);
}
