#include <gtest/gtest.h>
#include "algebra_linear/vector.hpp"

class VectorTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(VectorTest, ConstructorDefault) {
    Vector v;
    EXPECT_EQ(v.getLength(), 0);
}

TEST_F(VectorTest, ConstructorWithLength) {
    Vector v(5);
    EXPECT_EQ(v.getLength(), 5);
}

TEST_F(VectorTest, SetAndGetValue) {
    Vector v(4);
    v.setValue(0, 1.5);
    v.setValue(1, 2.5);
    v.setValue(2, 3.5);
    v.setValue(3, 4.5);
    
    EXPECT_DOUBLE_EQ(v.getValue(0), 1.5);
    EXPECT_DOUBLE_EQ(v.getValue(1), 2.5);
    EXPECT_DOUBLE_EQ(v.getValue(2), 3.5);
    EXPECT_DOUBLE_EQ(v.getValue(3), 4.5);
}

TEST_F(VectorTest, CopyConstructor) {
    Vector v1(3);
    v1.setValue(0, 10.5);
    v1.setValue(1, 20.0);
    v1.setValue(2, 30.0);
    
    Vector v2(v1);
    EXPECT_EQ(v2.getLength(), 3);
    EXPECT_DOUBLE_EQ(v2.getValue(0), 10.5);
    EXPECT_DOUBLE_EQ(v2.getValue(1), 20.0);
    EXPECT_DOUBLE_EQ(v2.getValue(2), 30.0);
}

TEST_F(VectorTest, Addition) {
    Vector v1(3);
    v1.setValue(0, 1.0);
    v1.setValue(1, 2.0);
    v1.setValue(2, 3.0);
    
    Vector v2(3);
    v2.setValue(0, 4.0);
    v2.setValue(1, 5.0);
    v2.setValue(2, 6.0);
    
    Vector v3 = v1 + v2;
    EXPECT_DOUBLE_EQ(v3.getValue(0), 5.0);
    EXPECT_DOUBLE_EQ(v3.getValue(1), 7.0);
    EXPECT_DOUBLE_EQ(v3.getValue(2), 9.0);
}

TEST_F(VectorTest, Subtraction) {
    Vector v1(3);
    v1.setValue(0, 5.0);
    v1.setValue(1, 7.0);
    v1.setValue(2, 9.0);
    
    Vector v2(3);
    v2.setValue(0, 2.0);
    v2.setValue(1, 3.0);
    v2.setValue(2, 4.0);
    
    Vector v3 = v1 - v2;
    EXPECT_DOUBLE_EQ(v3.getValue(0), 3.0);
    EXPECT_DOUBLE_EQ(v3.getValue(1), 4.0);
    EXPECT_DOUBLE_EQ(v3.getValue(2), 5.0);
}

TEST_F(VectorTest, ScalarMultiplication) {
    Vector v(3);
    v.setValue(0, 2.0);
    v.setValue(1, 3.0);
    v.setValue(2, 4.0);
    
    Vector v2 = v * 3;
    EXPECT_DOUBLE_EQ(v2.getValue(0), 6.0);
    EXPECT_DOUBLE_EQ(v2.getValue(1), 9.0);
    EXPECT_DOUBLE_EQ(v2.getValue(2), 12.0);
}

TEST_F(VectorTest, LinearProduct) {
    Vector v1(3);
    v1.setValue(0, 1.0);
    v1.setValue(1, 2.0);
    v1.setValue(2, 3.0);
    
    Vector v2(3);
    v2.setValue(0, 4.0);
    v2.setValue(1, 5.0);
    v2.setValue(2, 6.0);
    
    double result = v1.linear_product(v2);
    EXPECT_DOUBLE_EQ(result, 32.0);
}

TEST_F(VectorTest, ElementWiseProduct) {
    Vector v1(3);
    v1.setValue(0, 2.0);
    v1.setValue(1, 3.0);
    v1.setValue(2, 4.0);
    
    Vector v2(3);
    v2.setValue(0, 5.0);
    v2.setValue(1, 6.0);
    v2.setValue(2, 7.0);
    
    Vector v3 = v1 * v2;
    EXPECT_DOUBLE_EQ(v3.getValue(0), 10.0);
    EXPECT_DOUBLE_EQ(v3.getValue(1), 18.0);
    EXPECT_DOUBLE_EQ(v3.getValue(2), 28.0);
}

TEST_F(VectorTest, EuclidianLength) {
    Vector v(3);
    v.setValue(0, 3.0);
    v.setValue(1, 4.0);
    v.setValue(2, 0.0);
    
    double length = v.euclidian_length();
    EXPECT_DOUBLE_EQ(length, 5.0);
}

TEST_F(VectorTest, EuclidianLengthSqrt3) {
    Vector v(3);
    v.setValue(0, 1.0);
    v.setValue(1, 1.0);
    v.setValue(2, 1.0);
    
    double length = v.euclidian_length();
    EXPECT_NEAR(length, 1.732050808, 1e-8);
}

TEST_F(VectorTest, AXPYOperation) {
    Vector v1(3);
    v1.setValue(0, 1.0);
    v1.setValue(1, 2.0);
    v1.setValue(2, 3.0);
    
    Vector v2(3);
    v2.setValue(0, 4.0);
    v2.setValue(1, 5.0);
    v2.setValue(2, 6.0);
    
    // AXPY: v1 * alpha + v2 = v1 * 2.0 + v2 = [2,4,6] + [4,5,6] = [6,9,12]
    Vector result = v1.axpy(v2, 2.0);
    EXPECT_DOUBLE_EQ(result.getValue(0), 6.0);
    EXPECT_DOUBLE_EQ(result.getValue(1), 9.0);
    EXPECT_DOUBLE_EQ(result.getValue(2), 12.0);
}

TEST_F(VectorTest, Equality) {
    Vector v1(3);
    v1.setValue(0, 1.0);
    v1.setValue(1, 2.0);
    v1.setValue(2, 3.0);
    
    Vector v2(3);
    v2.setValue(0, 1.0);
    v2.setValue(1, 2.0);
    v2.setValue(2, 3.0);
    
    Vector v3(3);
    v3.setValue(0, 1.0);
    v3.setValue(1, 2.0);
    v3.setValue(2, 4.0);
    
    EXPECT_TRUE(v1.equality(v2));
    EXPECT_FALSE(v1.equality(v3));
}

TEST_F(VectorTest, MatrixVectorMultiplication) {
    Matriz A(2, 3);
    A.setValue(0, 0, 1.0);
    A.setValue(0, 1, 2.0);
    A.setValue(0, 2, 3.0);
    A.setValue(1, 0, 4.0);
    A.setValue(1, 1, 5.0);
    A.setValue(1, 2, 6.0);
    
    Vector v(3);
    v.setValue(0, 1.0);
    v.setValue(1, 2.0);
    v.setValue(2, 3.0);
    
    Vector result = multiplicar(A, v);
    EXPECT_DOUBLE_EQ(result.getValue(0), 14.0);
    EXPECT_DOUBLE_EQ(result.getValue(1), 32.0);
}
