#include <gtest/gtest.h>
#include "algebra_linear/matriz.hpp"
#include <cmath>

class MatrizTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create 2x2 matrix: [[1, 2], [3, 4]]
        m1 = Matriz(2, 2);
        m1.setValue(0, 0, 1.0);
        m1.setValue(0, 1, 2.0);
        m1.setValue(1, 0, 3.0);
        m1.setValue(1, 1, 4.0);

        // Create 2x2 matrix: [[5, 6], [7, 8]]
        m2 = Matriz(2, 2);
        m2.setValue(0, 0, 5.0);
        m2.setValue(0, 1, 6.0);
        m2.setValue(1, 0, 7.0);
        m2.setValue(1, 1, 8.0);
    }

    Matriz m1, m2;
};

TEST_F(MatrizTest, ConstructorDefault) {
    Matriz m;
    EXPECT_EQ(m.getRows(), 0);
    EXPECT_EQ(m.getColumns(), 0);
}

TEST_F(MatrizTest, ConstructorWithDimensions) {
    Matriz m(3, 4);
    EXPECT_EQ(m.getRows(), 3);
    EXPECT_EQ(m.getColumns(), 4);
}

TEST_F(MatrizTest, GetSetValue) {
    Matriz m(2, 2);
    m.setValue(0, 0, 5.5);
    EXPECT_EQ(m.getValue(0, 0), 5.5);
}

TEST_F(MatrizTest, GetValueOutOfRange) {
    Matriz m(2, 2);
    EXPECT_THROW(m.getValue(2, 2), std::out_of_range);
}

TEST_F(MatrizTest, SetValueOutOfRange) {
    Matriz m(2, 2);
    EXPECT_THROW(m.setValue(5, 5, 10.0), std::out_of_range);
}

TEST_F(MatrizTest, CopyConstructor) {
    Matriz copy = m1;
    EXPECT_EQ(copy.getRows(), m1.getRows());
    EXPECT_EQ(copy.getColumns(), m1.getColumns());
    EXPECT_EQ(copy.getValue(0, 0), m1.getValue(0, 0));
}

TEST_F(MatrizTest, MatrizAddition) {
    Matriz result = m1 + m2;
    EXPECT_EQ(result.getValue(0, 0), 6.0);
    EXPECT_EQ(result.getValue(0, 1), 8.0);
    EXPECT_EQ(result.getValue(1, 0), 10.0);
    EXPECT_EQ(result.getValue(1, 1), 12.0);
}

TEST_F(MatrizTest, MatrizSubtraction) {
    Matriz result = m2 - m1;
    EXPECT_EQ(result.getValue(0, 0), 4.0);
    EXPECT_EQ(result.getValue(0, 1), 4.0);
    EXPECT_EQ(result.getValue(1, 0), 4.0);
    EXPECT_EQ(result.getValue(1, 1), 4.0);
}

TEST_F(MatrizTest, MatrizScalarMultiplication) {
    Matriz result = m1 * 2.0;
    EXPECT_EQ(result.getValue(0, 0), 2.0);
    EXPECT_EQ(result.getValue(0, 1), 4.0);
    EXPECT_EQ(result.getValue(1, 0), 6.0);
    EXPECT_EQ(result.getValue(1, 1), 8.0);
}

TEST_F(MatrizTest, MatrizMultiplication) {
    Matriz result = m1 * m2;
    // [[1,2],[3,4]] * [[5,6],[7,8]] = [[19,22],[43,50]]
    EXPECT_EQ(result.getValue(0, 0), 19.0);
    EXPECT_EQ(result.getValue(0, 1), 22.0);
    EXPECT_EQ(result.getValue(1, 0), 43.0);
    EXPECT_EQ(result.getValue(1, 1), 50.0);
}

TEST_F(MatrizTest, IsSquareMatrix) {
    Matriz square(3, 3);
    Matriz notSquare(2, 3);
    EXPECT_TRUE(square.isSquare());
    EXPECT_FALSE(notSquare.isSquare());
}

TEST_F(MatrizTest, Transpose) {
    Matriz m(2, 3);
    m.setValue(0, 0, 1.0);
    m.setValue(0, 1, 2.0);
    m.setValue(0, 2, 3.0);
    m.setValue(1, 0, 4.0);
    m.setValue(1, 1, 5.0);
    m.setValue(1, 2, 6.0);

    Matriz t = m.transpose();
    EXPECT_EQ(t.getRows(), 3);
    EXPECT_EQ(t.getColumns(), 2);
    EXPECT_EQ(t.getValue(0, 0), 1.0);
    EXPECT_EQ(t.getValue(1, 0), 2.0);
    EXPECT_EQ(t.getValue(2, 0), 3.0);
}

TEST_F(MatrizTest, Clone) {
    Matriz clone = m1.clonar();
    EXPECT_EQ(clone.getRows(), m1.getRows());
    EXPECT_EQ(clone.getColumns(), m1.getColumns());
    EXPECT_EQ(clone.getValue(0, 0), m1.getValue(0, 0));
}

TEST_F(MatrizTest, IdentityMatrix) {
    Matriz identity = Matriz::identity(3, 3);
    EXPECT_EQ(identity.getRows(), 3);
    EXPECT_EQ(identity.getColumns(), 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (i == j) {
                EXPECT_EQ(identity.getValue(i, j), 1.0);
            } else {
                EXPECT_EQ(identity.getValue(i, j), 0.0);
            }
        }
    }
}

TEST_F(MatrizTest, TrocaLinhas) {
    Matriz m(3, 2);
    m.setValue(0, 0, 1.0);
    m.setValue(1, 0, 2.0);
    m.setValue(2, 0, 3.0);

    m.trocarLinhas(0, 2);
    EXPECT_EQ(m.getValue(0, 0), 3.0);
    EXPECT_EQ(m.getValue(2, 0), 1.0);
}

TEST_F(MatrizTest, AdditionIncompatibleDimensions) {
    Matriz other(3, 3);
    EXPECT_THROW(m1 + other, std::invalid_argument);
}

TEST_F(MatrizTest, SubtractionIncompatibleDimensions) {
    Matriz other(3, 3);
    EXPECT_THROW(m1 - other, std::invalid_argument);
}

TEST_F(MatrizTest, MultiplicationIncompatibleDimensions) {
    Matriz other(3, 3);
    EXPECT_THROW(m1 * other, std::invalid_argument);
}
