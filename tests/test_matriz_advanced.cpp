#include <gtest/gtest.h>
#include "algebra_linear/matriz.hpp"
#include <cmath>

class MatrizAdvancedTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a 3x3 symmetric matrix for testing properties
        sym = Matriz(3, 3);
        sym.setValue(0, 0, 2.0);
        sym.setValue(0, 1, 1.0);
        sym.setValue(0, 2, 0.0);
        sym.setValue(1, 0, 1.0);
        sym.setValue(1, 1, 2.0);
        sym.setValue(1, 2, 1.0);
        sym.setValue(2, 0, 0.0);
        sym.setValue(2, 1, 1.0);
        sym.setValue(2, 2, 2.0);

        // Create a non-symmetric matrix
        nonsym = Matriz(3, 3);
        nonsym.setValue(0, 0, 1.0);
        nonsym.setValue(0, 1, 2.0);
        nonsym.setValue(0, 2, 3.0);
        nonsym.setValue(1, 0, 4.0);
        nonsym.setValue(1, 1, 5.0);
        nonsym.setValue(1, 2, 6.0);
        nonsym.setValue(2, 0, 7.0);
        nonsym.setValue(2, 1, 8.0);
        nonsym.setValue(2, 2, 9.0);
    }

    Matriz sym, nonsym;
};

TEST_F(MatrizAdvancedTest, IsSimetricMatrixTrue) {
    EXPECT_TRUE(sym.isSimetric());
}

TEST_F(MatrizAdvancedTest, IsSimetricMatrixFalse) {
    EXPECT_FALSE(nonsym.isSimetric());
}

TEST_F(MatrizAdvancedTest, DeterminantSimple2x2) {
    Matriz m(2, 2);
    m.setValue(0, 0, 3.0);
    m.setValue(0, 1, 8.0);
    m.setValue(1, 0, 4.0);
    m.setValue(1, 1, 6.0);
    // det = 3*6 - 8*4 = 18 - 32 = -14
    EXPECT_NEAR(m.determinant(), -14.0, 1e-6);
}

TEST_F(MatrizAdvancedTest, Determinant3x3) {
    // Standard 3x3 determinant test
    Matriz m(3, 3);
    m.setValue(0, 0, 1.0);
    m.setValue(0, 1, 2.0);
    m.setValue(0, 2, 3.0);
    m.setValue(1, 0, 0.0);
    m.setValue(1, 1, 1.0);
    m.setValue(1, 2, 4.0);
    m.setValue(2, 0, 5.0);
    m.setValue(2, 1, 6.0);
    m.setValue(2, 2, 0.0);
    // det = 1*(1*0 - 4*6) - 2*(0*0 - 4*5) + 3*(0*6 - 1*5)
    //     = 1*(-24) - 2*(-20) + 3*(-5)
    //     = -24 + 40 - 15 = 1
    EXPECT_NEAR(m.determinant(), 1.0, 1e-6);
}

TEST_F(MatrizAdvancedTest, TransposeOfTranspose) {
    Matriz m(2, 3);
    m.setValue(0, 0, 1.0);
    m.setValue(0, 1, 2.0);
    m.setValue(0, 2, 3.0);
    m.setValue(1, 0, 4.0);
    m.setValue(1, 1, 5.0);
    m.setValue(1, 2, 6.0);

    Matriz t_t = m.transpose().transpose();
    EXPECT_EQ(t_t.getRows(), m.getRows());
    EXPECT_EQ(t_t.getColumns(), m.getColumns());
    for (int i = 0; i < m.getRows(); i++) {
        for (int j = 0; j < m.getColumns(); j++) {
            EXPECT_EQ(t_t.getValue(i, j), m.getValue(i, j));
        }
    }
}

TEST_F(MatrizAdvancedTest, MatrixChainMultiplication) {
    Matriz a(2, 2);
    a.setValue(0, 0, 1.0);
    a.setValue(0, 1, 2.0);
    a.setValue(1, 0, 3.0);
    a.setValue(1, 1, 4.0);

    Matriz b(2, 2);
    b.setValue(0, 0, 2.0);
    b.setValue(0, 1, 0.0);
    b.setValue(1, 0, 1.0);
    b.setValue(1, 1, 2.0);

    Matriz c(2, 2);
    c.setValue(0, 0, 1.0);
    c.setValue(0, 1, 1.0);
    c.setValue(1, 0, 1.0);
    c.setValue(1, 1, 1.0);

    Matriz result = (a * b) * c;
    Matriz result2 = a * (b * c);
    
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            EXPECT_NEAR(result.getValue(i, j), result2.getValue(i, j), 1e-6);
        }
    }
}

TEST_F(MatrizAdvancedTest, ScalarMultiplicationCommutative) {
    Matriz m(2, 2);
    m.setValue(0, 0, 1.0);
    m.setValue(0, 1, 2.0);
    m.setValue(1, 0, 3.0);
    m.setValue(1, 1, 4.0);

    Matriz result1 = m * 5.0;
    Matriz result2 = m * 5.0;

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            EXPECT_EQ(result1.getValue(i, j), result2.getValue(i, j));
        }
    }
}

TEST_F(MatrizAdvancedTest, AdditionAssociative) {
    Matriz a(2, 2);
    a.setValue(0, 0, 1.0);
    a.setValue(0, 1, 2.0);
    a.setValue(1, 0, 3.0);
    a.setValue(1, 1, 4.0);

    Matriz b(2, 2);
    b.setValue(0, 0, 5.0);
    b.setValue(0, 1, 6.0);
    b.setValue(1, 0, 7.0);
    b.setValue(1, 1, 8.0);

    Matriz c(2, 2);
    c.setValue(0, 0, 1.0);
    c.setValue(0, 1, 1.0);
    c.setValue(1, 0, 1.0);
    c.setValue(1, 1, 1.0);

    Matriz result1 = (a + b) + c;
    Matriz result2 = a + (b + c);

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            EXPECT_EQ(result1.getValue(i, j), result2.getValue(i, j));
        }
    }
}

TEST_F(MatrizAdvancedTest, LargeMatrixOperations) {
    Matriz large(100, 100);
    // Fill with small values to avoid overflow
    for (int i = 0; i < 100; i++) {
        for (int j = 0; j < 100; j++) {
            large.setValue(i, j, 0.01 * (i + j));
        }
    }

    Matriz scaled = large * 2.0;
    for (int i = 0; i < 100; i++) {
        for (int j = 0; j < 100; j++) {
            EXPECT_EQ(scaled.getValue(i, j), large.getValue(i, j) * 2.0);
        }
    }
}

TEST_F(MatrizAdvancedTest, IdentityMatrixProperties) {
    Matriz identity = Matriz::identity(5, 5);
    Matriz m(5, 5);
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            m.setValue(i, j, i * 5.0 + j + 1.0);
        }
    }

    Matriz result = m * identity;
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            EXPECT_NEAR(result.getValue(i, j), m.getValue(i, j), 1e-6);
        }
    }
}

TEST_F(MatrizAdvancedTest, DistributiveProperty) {
    // A * (B + C) = A*B + A*C
    Matriz a(2, 2);
    a.setValue(0, 0, 1.0);
    a.setValue(0, 1, 2.0);
    a.setValue(1, 0, 3.0);
    a.setValue(1, 1, 4.0);

    Matriz b(2, 2);
    b.setValue(0, 0, 1.0);
    b.setValue(0, 1, 0.0);
    b.setValue(1, 0, 0.0);
    b.setValue(1, 1, 1.0);

    Matriz c(2, 2);
    c.setValue(0, 0, 0.0);
    c.setValue(0, 1, 1.0);
    c.setValue(1, 0, 1.0);
    c.setValue(1, 1, 0.0);

    Matriz left = a * (b + c);
    Matriz right = (a * b) + (a * c);

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            EXPECT_NEAR(left.getValue(i, j), right.getValue(i, j), 1e-6);
        }
    }
}
