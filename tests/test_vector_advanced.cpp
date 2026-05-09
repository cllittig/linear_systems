#include <gtest/gtest.h>
#include <cmath>
#include "algebra_linear/vector.hpp"

class VectorAdvancedTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(VectorAdvancedTest, VectorZero) {
    Vector v(5);
    
    // All values should be zero by default
    for (int i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ(v.getValue(i), 0.0);
    }
}

TEST_F(VectorAdvancedTest, VectorNegativeValues) {
    Vector v(3);
    v.setValue(0, -1.5);
    v.setValue(1, -2.5);
    v.setValue(2, -3.5);
    
    EXPECT_DOUBLE_EQ(v.getValue(0), -1.5);
    EXPECT_DOUBLE_EQ(v.getValue(1), -2.5);
    EXPECT_DOUBLE_EQ(v.getValue(2), -3.5);
}

TEST_F(VectorAdvancedTest, AdditionWithNegatives) {
    Vector v1(2);
    v1.setValue(0, 5.0);
    v1.setValue(1, -3.0);
    
    Vector v2(2);
    v2.setValue(0, -2.0);
    v2.setValue(1, 3.0);
    
    Vector v3 = v1 + v2;
    EXPECT_DOUBLE_EQ(v3.getValue(0), 3.0);
    EXPECT_DOUBLE_EQ(v3.getValue(1), 0.0);
}

TEST_F(VectorAdvancedTest, SubtractionNegative) {
    Vector v1(3);
    v1.setValue(0, 1.0);
    v1.setValue(1, 2.0);
    v1.setValue(2, 3.0);
    
    Vector v2(3);
    v2.setValue(0, 5.0);
    v2.setValue(1, 6.0);
    v2.setValue(2, 7.0);
    
    Vector v3 = v1 - v2;
    EXPECT_DOUBLE_EQ(v3.getValue(0), -4.0);
    EXPECT_DOUBLE_EQ(v3.getValue(1), -4.0);
    EXPECT_DOUBLE_EQ(v3.getValue(2), -4.0);
}

TEST_F(VectorAdvancedTest, ScalarMultiplicationNegative) {
    Vector v(3);
    v.setValue(0, 1.0);
    v.setValue(1, 2.0);
    v.setValue(2, 3.0);
    
    Vector v2 = v * (-2);
    EXPECT_DOUBLE_EQ(v2.getValue(0), -2.0);
    EXPECT_DOUBLE_EQ(v2.getValue(1), -4.0);
    EXPECT_DOUBLE_EQ(v2.getValue(2), -6.0);
}

TEST_F(VectorAdvancedTest, LinearProductZero) {
    Vector v1(3);
    
    Vector v2(3);
    v2.setValue(0, 5.0);
    v2.setValue(1, 10.0);
    v2.setValue(2, 15.0);
    
    double product = v1.linear_product(v2);
    EXPECT_DOUBLE_EQ(product, 0.0);
}

TEST_F(VectorAdvancedTest, LinearProductNegative) {
    Vector v1(2);
    v1.setValue(0, 1.0);
    v1.setValue(1, 2.0);
    
    Vector v2(2);
    v2.setValue(0, -1.0);
    v2.setValue(1, -2.0);
    
    double product = v1.linear_product(v2);
    EXPECT_DOUBLE_EQ(product, -5.0);
}

TEST_F(VectorAdvancedTest, LinearProductOrthogonal) {
    Vector v1(3);
    v1.setValue(0, 1.0);
    v1.setValue(1, 0.0);
    v1.setValue(2, 0.0);
    
    Vector v2(3);
    v2.setValue(0, 0.0);
    v2.setValue(1, 1.0);
    v2.setValue(2, 0.0);
    
    double product = v1.linear_product(v2);
    EXPECT_DOUBLE_EQ(product, 0.0);
}

TEST_F(VectorAdvancedTest, EuclidianLengthZero) {
    Vector v(4);
    double length = v.euclidian_length();
    EXPECT_DOUBLE_EQ(length, 0.0);
}

TEST_F(VectorAdvancedTest, EuclidianLengthNegativeValues) {
    Vector v(2);
    v.setValue(0, -3.0);
    v.setValue(1, -4.0);
    
    double length = v.euclidian_length();
    EXPECT_DOUBLE_EQ(length, 5.0);
}

TEST_F(VectorAdvancedTest, EuclidianLengthComplexNumbers) {
    Vector v(3);
    v.setValue(0, 2.0);
    v.setValue(1, 2.0);
    v.setValue(2, 1.0);
    
    double length = v.euclidian_length();
    EXPECT_NEAR(length, 3.0, 1e-9);
}

TEST_F(VectorAdvancedTest, AXPYZeroAlpha) {
    Vector v1(3);
    v1.setValue(0, 1.0);
    v1.setValue(1, 2.0);
    v1.setValue(2, 3.0);
    
    Vector v2(3);
    v2.setValue(0, 4.0);
    v2.setValue(1, 5.0);
    v2.setValue(2, 6.0);
    
    // AXPY with alpha = 0: v1 * 0 + v2 = v2
    Vector result = v1.axpy(v2, 0.0);
    
    EXPECT_DOUBLE_EQ(result.getValue(0), 4.0);
    EXPECT_DOUBLE_EQ(result.getValue(1), 5.0);
    EXPECT_DOUBLE_EQ(result.getValue(2), 6.0);
}

TEST_F(VectorAdvancedTest, AXPYNegativeAlpha) {
    Vector v1(2);
    v1.setValue(0, 2.0);
    v1.setValue(1, 3.0);
    
    Vector v2(2);
    v2.setValue(0, 1.0);
    v2.setValue(1, 1.0);
    
    // v1 * (-1) + v2 = [-2, -3] + [1, 1] = [-1, -2]
    Vector result = v1.axpy(v2, -1.0);
    
    EXPECT_DOUBLE_EQ(result.getValue(0), -1.0);
    EXPECT_DOUBLE_EQ(result.getValue(1), -2.0);
}

TEST_F(VectorAdvancedTest, EqualityDifferentOrder) {
    Vector v1(3);
    v1.setValue(0, 1.0);
    v1.setValue(1, 2.0);
    v1.setValue(2, 3.0);
    
    Vector v2(3);
    v2.setValue(0, 1.0);
    v2.setValue(1, 2.0);
    v2.setValue(2, 3.0);
    
    EXPECT_TRUE(v1.equality(v2));
    EXPECT_TRUE(v2.equality(v1));
}

TEST_F(VectorAdvancedTest, EqualitySmallDifference) {
    Vector v1(2);
    v1.setValue(0, 1.0);
    v1.setValue(1, 2.0);
    
    Vector v2(2);
    v2.setValue(0, 1.0000001);
    v2.setValue(1, 2.0);
    
    EXPECT_FALSE(v1.equality(v2));
}

TEST_F(VectorAdvancedTest, ElementWiseProductZero) {
    Vector v1(3);
    
    Vector v2(3);
    v2.setValue(0, 5.0);
    v2.setValue(1, 10.0);
    v2.setValue(2, 15.0);
    
    Vector result = v1 * v2;
    
    EXPECT_DOUBLE_EQ(result.getValue(0), 0.0);
    EXPECT_DOUBLE_EQ(result.getValue(1), 0.0);
    EXPECT_DOUBLE_EQ(result.getValue(2), 0.0);
}

TEST_F(VectorAdvancedTest, ElementWiseProductNegative) {
    Vector v1(2);
    v1.setValue(0, -2.0);
    v1.setValue(1, 3.0);
    
    Vector v2(2);
    v2.setValue(0, 4.0);
    v2.setValue(1, -1.0);
    
    Vector result = v1 * v2;
    
    EXPECT_DOUBLE_EQ(result.getValue(0), -8.0);
    EXPECT_DOUBLE_EQ(result.getValue(1), -3.0);
}

TEST_F(VectorAdvancedTest, LargeVector) {
    Vector v(100);
    for (int i = 0; i < 100; i++) {
        v.setValue(i, static_cast<double>(i));
    }
    
    EXPECT_EQ(v.getLength(), 100);
    EXPECT_DOUBLE_EQ(v.getValue(50), 50.0);
    EXPECT_DOUBLE_EQ(v.getValue(99), 99.0);
}

TEST_F(VectorAdvancedTest, LargeVectorAddition) {
    Vector v1(50);
    Vector v2(50);
    
    for (int i = 0; i < 50; i++) {
        v1.setValue(i, static_cast<double>(i));
        v2.setValue(i, static_cast<double>(i));
    }
    
    Vector v3 = v1 + v2;
    
    for (int i = 0; i < 50; i++) {
        EXPECT_DOUBLE_EQ(v3.getValue(i), static_cast<double>(2 * i));
    }
}

TEST_F(VectorAdvancedTest, VectorCopyModification) {
    Vector v1(3);
    v1.setValue(0, 1.0);
    v1.setValue(1, 2.0);
    v1.setValue(2, 3.0);
    
    Vector v2(v1);
    v2.setValue(0, 999.0);
    
    EXPECT_DOUBLE_EQ(v1.getValue(0), 1.0);
    EXPECT_DOUBLE_EQ(v2.getValue(0), 999.0);
}

TEST_F(VectorAdvancedTest, ArithmeticChain) {
    Vector v1(2);
    v1.setValue(0, 1.0);
    v1.setValue(1, 2.0);
    
    Vector v2(2);
    v2.setValue(0, 3.0);
    v2.setValue(1, 4.0);
    
    Vector v3(2);
    v3.setValue(0, 5.0);
    v3.setValue(1, 6.0);
    
    // (v1 + v2) * 2 - v3
    Vector result = ((v1 + v2) * 2) - v3;
    
    EXPECT_DOUBLE_EQ(result.getValue(0), 3.0);  // (1+3)*2 - 5 = 8 - 5 = 3
    EXPECT_DOUBLE_EQ(result.getValue(1), 6.0);  // (2+4)*2 - 6 = 12 - 6 = 6
}

TEST_F(VectorAdvancedTest, MatrixVectorMultiplicationZero) {
    Matriz A(2, 3);
    
    Vector v(3);
    v.setValue(0, 1.0);
    v.setValue(1, 2.0);
    v.setValue(2, 3.0);
    
    Vector result = multiplicar(A, v);
    
    EXPECT_EQ(result.getLength(), 2);
    EXPECT_DOUBLE_EQ(result.getValue(0), 0.0);
    EXPECT_DOUBLE_EQ(result.getValue(1), 0.0);
}

TEST_F(VectorAdvancedTest, MatrixVectorMultiplicationLarge) {
    Matriz A(3, 4);
    A.setValue(0, 0, 1.0);
    A.setValue(0, 1, 2.0);
    A.setValue(0, 2, 3.0);
    A.setValue(0, 3, 4.0);
    A.setValue(1, 0, 5.0);
    A.setValue(1, 1, 6.0);
    A.setValue(1, 2, 7.0);
    A.setValue(1, 3, 8.0);
    A.setValue(2, 0, 9.0);
    A.setValue(2, 1, 10.0);
    A.setValue(2, 2, 11.0);
    A.setValue(2, 3, 12.0);
    
    Vector v(4);
    v.setValue(0, 1.0);
    v.setValue(1, 1.0);
    v.setValue(2, 1.0);
    v.setValue(3, 1.0);
    
    Vector result = multiplicar(A, v);
    
    EXPECT_DOUBLE_EQ(result.getValue(0), 10.0);  // 1+2+3+4
    EXPECT_DOUBLE_EQ(result.getValue(1), 26.0);  // 5+6+7+8
    EXPECT_DOUBLE_EQ(result.getValue(2), 42.0);  // 9+10+11+12
}
