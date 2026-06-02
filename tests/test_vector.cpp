#include <gtest/gtest.h>
#include <stdexcept>
#include "algebra_linear/vector.hpp"

/**
 * @file test_vector.cpp
 * @brief Testes unitários da classe base Vector (GoogleTest).
 *
 * Cobre construção, acesso (com bound checks, que em Vector estão sempre
 * ativos), semântica de cópia, operadores (soma, subtração, produto por
 * escalar, produto componente a componente), produto interno, norma
 * euclidiana, axpy e igualdade, incluindo os casos de erro.
 */

static const double TOL = 1e-12;

namespace {

Vector make3(double a, double b, double c) {
    Vector v(3);
    v.setValue(0, a);
    v.setValue(1, b);
    v.setValue(2, c);
    return v;
}

} // namespace

// --------------------------------------------------------------------------
// Construção e acesso
// --------------------------------------------------------------------------

TEST(VectorTest, DefaultConstructorIsEmpty) {
    Vector v;
    EXPECT_EQ(v.getLength(), 0);
}

TEST(VectorTest, DimensionConstructor) {
    Vector v(5);
    EXPECT_EQ(v.getLength(), 5);
    for (int i = 0; i < 5; ++i)
        EXPECT_NEAR(v.getValue(i), 0.0, TOL);  // inicializado em zero
}

TEST(VectorTest, NegativeLengthThrows) {
    EXPECT_THROW(Vector(-1), std::invalid_argument);
}

TEST(VectorTest, GetSetValue) {
    Vector v(3);
    v.setValue(1, 10.5);
    EXPECT_NEAR(v.getValue(1), 10.5, TOL);
}

TEST(VectorTest, GetValueOutOfRangeThrows) {
    Vector v(3);
    EXPECT_THROW(v.getValue(3), std::out_of_range);
    EXPECT_THROW(v.getValue(-1), std::out_of_range);
}

TEST(VectorTest, SetValueOutOfRangeThrows) {
    Vector v(3);
    EXPECT_THROW(v.setValue(10, 5.0), std::out_of_range);
}

TEST(VectorTest, CopyConstructorIsIndependent) {
    Vector original = make3(1, 2, 3);
    Vector copy = original;        // copy ctor (regra do zero)
    copy.setValue(0, 99.0);
    EXPECT_NEAR(original.getValue(0), 1.0, TOL);
    EXPECT_NEAR(copy.getValue(0), 99.0, TOL);
}

TEST(VectorTest, CopyAssignmentIsIndependent) {
    Vector original = make3(1, 2, 3);
    Vector other(3);
    other = original;              // copy assignment (regra do zero)
    other.setValue(2, -7.0);
    EXPECT_NEAR(original.getValue(2), 3.0, TOL);
    EXPECT_NEAR(other.getValue(2), -7.0, TOL);
}

// --------------------------------------------------------------------------
// Operadores
// --------------------------------------------------------------------------

TEST(VectorTest, Addition) {
    Vector r = make3(1, 2, 3) + make3(4, 5, 6);
    EXPECT_NEAR(r.getValue(0), 5.0, TOL);
    EXPECT_NEAR(r.getValue(1), 7.0, TOL);
    EXPECT_NEAR(r.getValue(2), 9.0, TOL);
}

TEST(VectorTest, AdditionIncompatibleThrows) {
    EXPECT_THROW(make3(1, 2, 3) + Vector(5), std::invalid_argument);
}

TEST(VectorTest, Subtraction) {
    Vector r = make3(4, 5, 6) - make3(1, 2, 3);
    EXPECT_NEAR(r.getValue(0), 3.0, TOL);
    EXPECT_NEAR(r.getValue(1), 3.0, TOL);
    EXPECT_NEAR(r.getValue(2), 3.0, TOL);
}

TEST(VectorTest, ScalarMultiplication) {
    Vector r = make3(1, 2, 3) * 2.0;
    EXPECT_NEAR(r.getValue(0), 2.0, TOL);
    EXPECT_NEAR(r.getValue(1), 4.0, TOL);
    EXPECT_NEAR(r.getValue(2), 6.0, TOL);
}

TEST(VectorTest, ScalarMultiplicationByZero) {
    Vector r = make3(1, 2, 3) * 0.0;
    for (int i = 0; i < 3; ++i)
        EXPECT_NEAR(r.getValue(i), 0.0, TOL);
}

TEST(VectorTest, ComponentWiseProduct) {
    // [1,2,3] .* [4,5,6] = [4,10,18]
    Vector r = make3(1, 2, 3) * make3(4, 5, 6);
    EXPECT_NEAR(r.getValue(0), 4.0, TOL);
    EXPECT_NEAR(r.getValue(1), 10.0, TOL);
    EXPECT_NEAR(r.getValue(2), 18.0, TOL);
}

// --------------------------------------------------------------------------
// Produto interno, norma e axpy
// --------------------------------------------------------------------------

TEST(VectorTest, DotProduct) {
    // 1*4 + 2*5 + 3*6 = 32
    EXPECT_NEAR(make3(1, 2, 3).linear_product(make3(4, 5, 6)), 32.0, TOL);
}

TEST(VectorTest, DotProductWithZeroVector) {
    EXPECT_NEAR(Vector(3).linear_product(make3(1, 2, 3)), 0.0, TOL);
}

TEST(VectorTest, DotProductIncompatibleThrows) {
    EXPECT_THROW(make3(1, 2, 3).linear_product(Vector(2)),
                 std::invalid_argument);
}

TEST(VectorTest, DotProductWithItselfEqualsNormSquared) {
    Vector v = make3(1, 2, 3);
    double norm = v.euclidian_length();
    EXPECT_NEAR(v.linear_product(v), norm * norm, 1e-9);
}

TEST(VectorTest, EuclidianLength) {
    // ||(3,4,0)|| = 5
    EXPECT_NEAR(make3(3, 4, 0).euclidian_length(), 5.0, 1e-9);
}

TEST(VectorTest, EuclidianLengthUnit) {
    EXPECT_NEAR(make3(1, 0, 0).euclidian_length(), 1.0, 1e-9);
}

TEST(VectorTest, Axpy) {
    // axpy(y, alpha) = alpha * this + y : 2*[1,2,3] + [4,5,6] = [6,9,12]
    Vector r = make3(1, 2, 3).axpy(make3(4, 5, 6), 2.0);
    EXPECT_NEAR(r.getValue(0), 6.0, TOL);
    EXPECT_NEAR(r.getValue(1), 9.0, TOL);
    EXPECT_NEAR(r.getValue(2), 12.0, TOL);
}

TEST(VectorTest, AxpyZeroAlpha) {
    // 0 * this + y = y
    Vector r = make3(1, 2, 3).axpy(make3(4, 5, 6), 0.0);
    EXPECT_NEAR(r.getValue(0), 4.0, TOL);
    EXPECT_NEAR(r.getValue(1), 5.0, TOL);
    EXPECT_NEAR(r.getValue(2), 6.0, TOL);
}

// --------------------------------------------------------------------------
// Igualdade e estresse
// --------------------------------------------------------------------------

TEST(VectorTest, Equality) {
    Vector v = make3(1, 2, 3);
    EXPECT_TRUE(v.equality(make3(1, 2, 3)));
    EXPECT_FALSE(v.equality(make3(4, 5, 6)));
    EXPECT_FALSE(v.equality(Vector(2)));  // comprimentos diferentes
}

TEST(VectorTest, LargeVectorOperations) {
    const int n = 1000;
    Vector large(n);
    for (int i = 0; i < n; ++i)
        large.setValue(i, 0.001 * i);
    Vector doubled = large * 2.0;
    for (int i = 0; i < n; ++i)
        EXPECT_NEAR(doubled.getValue(i), 0.002 * i, 1e-10);
}
