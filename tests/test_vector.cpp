#include <gtest/gtest.h>
#include "algebra_linear/vector.hpp"
#include <cmath>

class VectorTest : public ::testing::Test {
protected:
    void SetUp() override {
        v1 = Vector(3);
        v1.setValue(0, 1.0);
        v1.setValue(1, 2.0);
        v1.setValue(2, 3.0);

        v2 = Vector(3);
        v2.setValue(0, 4.0);
        v2.setValue(1, 5.0);
        v2.setValue(2, 6.0);
    }

    Vector v1, v2;
};

TEST_F(VectorTest, ConstructorDefault) {
    Vector v;
    EXPECT_EQ(v.getLength(), 0);
}

TEST_F(VectorTest, ConstructorWithDimension) {
    Vector v(5);
    EXPECT_EQ(v.getLength(), 5);
}

TEST_F(VectorTest, GetSetValue) {
    Vector v(3);
    v.setValue(1, 10.5);
    EXPECT_EQ(v.getValue(1), 10.5);
}

TEST_F(VectorTest, GetValueOutOfRange) {
    Vector v(3);
    EXPECT_THROW(v.getValue(5), std::out_of_range);
}

TEST_F(VectorTest, SetValueOutOfRange) {
    Vector v(3);
    EXPECT_THROW(v.setValue(10, 5.0), std::out_of_range);
}

TEST_F(VectorTest, CopyConstructor) {
    Vector copy = v1;
    EXPECT_EQ(copy.getLength(), v1.getLength());
    for (int i = 0; i < v1.getLength(); i++) {
        EXPECT_EQ(copy.getValue(i), v1.getValue(i));
    }
}

TEST_F(VectorTest, VectorAddition) {
    Vector result = v1 + v2;
    EXPECT_EQ(result.getLength(), 3);
    EXPECT_EQ(result.getValue(0), 5.0);
    EXPECT_EQ(result.getValue(1), 7.0);
    EXPECT_EQ(result.getValue(2), 9.0);
}

TEST_F(VectorTest, VectorSubtraction) {
    Vector result = v2 - v1;
    EXPECT_EQ(result.getLength(), 3);
    EXPECT_EQ(result.getValue(0), 3.0);
    EXPECT_EQ(result.getValue(1), 3.0);
    EXPECT_EQ(result.getValue(2), 3.0);
}

TEST_F(VectorTest, ScalarMultiplication) {
    Vector result = v1 * 2;
    EXPECT_EQ(result.getLength(), 3);
    EXPECT_EQ(result.getValue(0), 2.0);
    EXPECT_EQ(result.getValue(1), 4.0);
    EXPECT_EQ(result.getValue(2), 6.0);
}

TEST_F(VectorTest, ScalarMultiplicationByZero) {
    Vector result = v1 * 0;
    EXPECT_EQ(result.getValue(0), 0.0);
    EXPECT_EQ(result.getValue(1), 0.0);
    EXPECT_EQ(result.getValue(2), 0.0);
}

TEST_F(VectorTest, DotProduct) {
    double product = v1.linear_product(v2);
    // 1*4 + 2*5 + 3*6 = 4 + 10 + 18 = 32
    EXPECT_EQ(product, 32.0);
}

TEST_F(VectorTest, DotProductZeroVectors) {
    Vector zero(3);
    zero.setValue(0, 0.0);
    zero.setValue(1, 0.0);
    zero.setValue(2, 0.0);
    double product = zero.linear_product(v1);
    EXPECT_EQ(product, 0.0);
}

TEST_F(VectorTest, EuclidianLength) {
    Vector v(3);
    v.setValue(0, 3.0);
    v.setValue(1, 4.0);
    v.setValue(2, 0.0);
    double norm = v.euclidian_length();
    EXPECT_NEAR(norm, 5.0, 1e-6);
}

TEST_F(VectorTest, EuclidianLengthStandardBasis) {
    Vector e1(3);
    e1.setValue(0, 1.0);
    e1.setValue(1, 0.0);
    e1.setValue(2, 0.0);
    EXPECT_NEAR(e1.euclidian_length(), 1.0, 1e-6);
}

TEST_F(VectorTest, VectorEquality) {
    Vector v3 = v1;
    EXPECT_TRUE(v1.equality(v3));
    EXPECT_FALSE(v1.equality(v2));
}

TEST_F(VectorTest, AdditionIncompatibleDimensions) {
    Vector v3(5);
    EXPECT_THROW(v1 + v3, std::invalid_argument);
}

TEST_F(VectorTest, SubtractionIncompatibleDimensions) {
    Vector v3(5);
    EXPECT_THROW(v1 - v3, std::invalid_argument);
}

TEST_F(VectorTest, Axpy) {
    // result = v1 + alpha * v2
    Vector result = v1.axpy(v2, 2.0);
    // [1,2,3] + 2*[4,5,6] = [1+8, 2+10, 3+12] = [9, 12, 15]
    EXPECT_EQ(result.getValue(0), 9.0);
    EXPECT_EQ(result.getValue(1), 12.0);
    EXPECT_EQ(result.getValue(2), 15.0);
}

TEST_F(VectorTest, AxpyZeroAlpha) {
    Vector result = v1.axpy(v2, 0.0);
    EXPECT_EQ(result.getValue(0), v1.getValue(0));
    EXPECT_EQ(result.getValue(1), v1.getValue(1));
    EXPECT_EQ(result.getValue(2), v1.getValue(2));
}

TEST_F(VectorTest, DotProductWithItself) {
    double product = v1.linear_product(v1);
    double norm = v1.euclidian_length();
    EXPECT_NEAR(product, norm * norm, 1e-6);
}

TEST_F(VectorTest, LargeVectorOperations) {
    Vector large(1000);
    for (int i = 0; i < 1000; i++) {
        large.setValue(i, 0.001 * i);
    }

    Vector doubled = large * 2;
    for (int i = 0; i < 1000; i++) {
        EXPECT_NEAR(doubled.getValue(i), 0.002 * i, 1e-10);
    }
}

TEST_F(VectorTest, ComponentWiseProduct) {
    Vector result = v1 * v2;
    // [1,2,3] * [4,5,6] = [4, 10, 18]
    EXPECT_EQ(result.getValue(0), 4.0);
    EXPECT_EQ(result.getValue(1), 10.0);
    EXPECT_EQ(result.getValue(2), 18.0);
}
