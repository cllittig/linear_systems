#include <gtest/gtest.h>
#include "algebra_linear/matriz.hpp"

class MatrizTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(MatrizTest, ConstructorDefault) {
    Matriz m;
    EXPECT_EQ(m.getRows(), 0);
    EXPECT_EQ(m.getColumns(), 0);
}

TEST_F(MatrizTest, ConstructorWithDimensions) {
    Matriz m(3, 3);
    EXPECT_EQ(m.getRows(), 3);
    EXPECT_EQ(m.getColumns(), 3);
}

TEST_F(MatrizTest, SetAndGetValue) {
    Matriz m(3, 3);
    m.setValue(0, 0, 10.5);
    m.setValue(1, 1, 20.3);
    m.setValue(2, 2, 30.7);
    
    EXPECT_DOUBLE_EQ(m.getValue(0, 0), 10.5);
    EXPECT_DOUBLE_EQ(m.getValue(1, 1), 20.3);
    EXPECT_DOUBLE_EQ(m.getValue(2, 2), 30.7);
}

TEST_F(MatrizTest, CopyConstructor) {
    Matriz m1(2, 2);
    m1.setValue(0, 0, 5.0);
    m1.setValue(0, 1, 6.0);
    
    Matriz m2(m1);
    EXPECT_EQ(m2.getRows(), 2);
    EXPECT_EQ(m2.getColumns(), 2);
    EXPECT_DOUBLE_EQ(m2.getValue(0, 0), 5.0);
    EXPECT_DOUBLE_EQ(m2.getValue(0, 1), 6.0);
}

TEST_F(MatrizTest, Addition) {
    Matriz m1(2, 2);
    m1.setValue(0, 0, 1.0);
    m1.setValue(0, 1, 2.0);
    m1.setValue(1, 0, 3.0);
    m1.setValue(1, 1, 4.0);
    
    Matriz m2(2, 2);
    m2.setValue(0, 0, 1.0);
    m2.setValue(0, 1, 1.0);
    m2.setValue(1, 0, 1.0);
    m2.setValue(1, 1, 1.0);
    
    Matriz m3 = m1 + m2;
    EXPECT_DOUBLE_EQ(m3.getValue(0, 0), 2.0);
    EXPECT_DOUBLE_EQ(m3.getValue(0, 1), 3.0);
    EXPECT_DOUBLE_EQ(m3.getValue(1, 0), 4.0);
    EXPECT_DOUBLE_EQ(m3.getValue(1, 1), 5.0);
}

TEST_F(MatrizTest, Subtraction) {
    Matriz m1(2, 2);
    m1.setValue(0, 0, 5.0);
    m1.setValue(0, 1, 6.0);
    m1.setValue(1, 0, 7.0);
    m1.setValue(1, 1, 8.0);
    
    Matriz m2(2, 2);
    m2.setValue(0, 0, 1.0);
    m2.setValue(0, 1, 1.0);
    m2.setValue(1, 0, 1.0);
    m2.setValue(1, 1, 1.0);
    
    Matriz m3 = m1 - m2;
    EXPECT_DOUBLE_EQ(m3.getValue(0, 0), 4.0);
    EXPECT_DOUBLE_EQ(m3.getValue(0, 1), 5.0);
    EXPECT_DOUBLE_EQ(m3.getValue(1, 0), 6.0);
    EXPECT_DOUBLE_EQ(m3.getValue(1, 1), 7.0);
}

TEST_F(MatrizTest, ScalarMultiplication) {
    Matriz m(2, 2);
    m.setValue(0, 0, 1.0);
    m.setValue(0, 1, 2.0);
    m.setValue(1, 0, 3.0);
    m.setValue(1, 1, 4.0);
    
    Matriz m2 = m * 2.0;
    EXPECT_DOUBLE_EQ(m2.getValue(0, 0), 2.0);
    EXPECT_DOUBLE_EQ(m2.getValue(0, 1), 4.0);
    EXPECT_DOUBLE_EQ(m2.getValue(1, 0), 6.0);
    EXPECT_DOUBLE_EQ(m2.getValue(1, 1), 8.0);
}

TEST_F(MatrizTest, MatrixMultiplication) {
    Matriz A(2, 3);
    A.setValue(0, 0, 1.0);
    A.setValue(0, 1, 2.0);
    A.setValue(0, 2, 3.0);
    A.setValue(1, 0, 4.0);
    A.setValue(1, 1, 5.0);
    A.setValue(1, 2, 6.0);
    
    Matriz B(3, 2);
    B.setValue(0, 0, 7.0);
    B.setValue(0, 1, 8.0);
    B.setValue(1, 0, 9.0);
    B.setValue(1, 1, 10.0);
    B.setValue(2, 0, 11.0);
    B.setValue(2, 1, 12.0);
    
    Matriz C = A * B;
    EXPECT_DOUBLE_EQ(C.getValue(0, 0), 58.0);
    EXPECT_DOUBLE_EQ(C.getValue(0, 1), 64.0);
    EXPECT_DOUBLE_EQ(C.getValue(1, 0), 139.0);
    EXPECT_DOUBLE_EQ(C.getValue(1, 1), 154.0);
}

TEST_F(MatrizTest, Transpose) {
    Matriz m(2, 3);
    m.setValue(0, 0, 1.0);
    m.setValue(0, 1, 2.0);
    m.setValue(0, 2, 3.0);
    m.setValue(1, 0, 4.0);
    m.setValue(1, 1, 5.0);
    m.setValue(1, 2, 6.0);
    
    Matriz mT = m.transpose();
    EXPECT_EQ(mT.getRows(), 3);
    EXPECT_EQ(mT.getColumns(), 2);
    EXPECT_DOUBLE_EQ(mT.getValue(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(mT.getValue(0, 1), 4.0);
    EXPECT_DOUBLE_EQ(mT.getValue(1, 0), 2.0);
    EXPECT_DOUBLE_EQ(mT.getValue(1, 1), 5.0);
}

TEST_F(MatrizTest, IsSquare) {
    Matriz m1(3, 3);
    EXPECT_TRUE(m1.isSquare());
    
    Matriz m2(2, 3);
    EXPECT_FALSE(m2.isSquare());
}

TEST_F(MatrizTest, IsSymmetric) {
    Matriz m1(3, 3);
    m1.setValue(0, 0, 4.0);
    m1.setValue(0, 1, 12.0);
    m1.setValue(0, 2, -16.0);
    m1.setValue(1, 0, 12.0);
    m1.setValue(1, 1, 37.0);
    m1.setValue(1, 2, -43.0);
    m1.setValue(2, 0, -16.0);
    m1.setValue(2, 1, -43.0);
    m1.setValue(2, 2, 98.0);
    
    EXPECT_TRUE(m1.isSimetric());
    
    Matriz m2(2, 2);
    m2.setValue(0, 0, 1.0);
    m2.setValue(0, 1, 2.0);
    m2.setValue(1, 0, 3.0);
    m2.setValue(1, 1, 4.0);
    
    EXPECT_FALSE(m2.isSimetric());
}

TEST_F(MatrizTest, IdentityMatrix) {
    Matriz I = Matriz::identity(3, 3);
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (i == j) {
                EXPECT_DOUBLE_EQ(I.getValue(i, j), 1.0) 
                    << "Identity matrix diagonal at [" << i << "][" << j << "]";
            } else {
                EXPECT_DOUBLE_EQ(I.getValue(i, j), 0.0)
                    << "Identity matrix off-diagonal at [" << i << "][" << j << "]";
            }
        }
    }
}
