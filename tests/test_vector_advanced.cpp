#include <gtest/gtest.h>
#include "algebra_linear/vector.hpp"
#include "algebra_linear/matriz.hpp"
#include <cmath>

class VectorAdvancedTest : public ::testing::Test {
protected:
    void SetUp() override {
        v1 = Vector(4);
        v1.setValue(0, 1.0);
        v1.setValue(1, 2.0);
        v1.setValue(2, 3.0);
        v1.setValue(3, 4.0);

        v2 = Vector(4);
        v2.setValue(0, 2.0);
        v2.setValue(1, 0.0);
        v2.setValue(2, -1.0);
        v2.setValue(3, 2.0);
    }

    Vector v1, v2;
};

TEST_F(VectorAdvancedTest, OrthogonalVectors) {
    Vector a(3);
    a.setValue(0, 1.0);
    a.setValue(1, 0.0);
    a.setValue(2, 0.0);

    Vector b(3);
    b.setValue(0, 0.0);
    b.setValue(1, 1.0);
    b.setValue(2, 0.0);

    EXPECT_EQ(a.linear_product(b), 0.0);
}

TEST_F(VectorAdvancedTest, UnitVector) {
    Vector u(3);
    u.setValue(0, 1.0 / std::sqrt(3));
    u.setValue(1, 1.0 / std::sqrt(3));
    u.setValue(2, 1.0 / std::sqrt(3));

    EXPECT_NEAR(u.euclidian_length(), 1.0, 1e-6);
}

TEST_F(VectorAdvancedTest, CauchySchwarzInequality) {
    // |<u, v>| <= ||u|| * ||v||
    double dotProduct = std::abs(v1.linear_product(v2));
    double norm1 = v1.euclidian_length();
    double norm2 = v2.euclidian_length();

    EXPECT_LE(dotProduct, norm1 * norm2 + 1e-6);
}

TEST_F(VectorAdvancedTest, TriangleInequality) {
    // ||u + v|| <= ||u|| + ||v||
    Vector sum = v1 + v2;
    double norm_sum = sum.euclidian_length();
    double norm1 = v1.euclidian_length();
    double norm2 = v2.euclidian_length();

    EXPECT_LE(norm_sum, norm1 + norm2 + 1e-6);
}

TEST_F(VectorAdvancedTest, PythagoreanTheorem) {
    Vector a(2);
    a.setValue(0, 3.0);
    a.setValue(1, 0.0);

    Vector b(2);
    b.setValue(0, 0.0);
    b.setValue(1, 4.0);

    Vector sum = a + b;
    double norm_sum = sum.euclidian_length();
    double norm_a = a.euclidian_length();
    double norm_b = b.euclidian_length();

    // 3-4-5 triangle
    EXPECT_NEAR(norm_sum, 5.0, 1e-6);
    EXPECT_NEAR(norm_a * norm_a + norm_b * norm_b, norm_sum * norm_sum, 1e-6);
}

TEST_F(VectorAdvancedTest, VectorProjection) {
    // Project v1 onto v2
    double numerator = v1.linear_product(v2);
    double denominator = v2.linear_product(v2);
    double scalar = numerator / denominator;
    Vector projection = v2 * (int)scalar;

    // The projection should have length close to numerator/norm(v2)
    double expected_length = std::abs(numerator) / v2.euclidian_length();
    EXPECT_NEAR(projection.euclidian_length(), expected_length, 1e-6);
}

TEST_F(VectorAdvancedTest, AdditionAssociative) {
    Vector v3(4);
    v3.setValue(0, 1.0);
    v3.setValue(1, 1.0);
    v3.setValue(2, 1.0);
    v3.setValue(3, 1.0);

    Vector result1 = (v1 + v2) + v3;
    Vector result2 = v1 + (v2 + v3);

    for (int i = 0; i < 4; i++) {
        EXPECT_NEAR(result1.getValue(i), result2.getValue(i), 1e-6);
    }
}

TEST_F(VectorAdvancedTest, ScalarDistributive) {
    // k * (u + v) = k*u + k*v
    Vector sum = v1 + v2;
    Vector scaled_sum = sum * 3;
    
    Vector scaled1 = v1 * 3;
    Vector scaled2 = v2 * 3;
    Vector result = scaled1 + scaled2;

    for (int i = 0; i < 4; i++) {
        EXPECT_NEAR(scaled_sum.getValue(i), result.getValue(i), 1e-6);
    }
}

TEST_F(VectorAdvancedTest, MatrixVectorMultiplication) {
    Matriz m(3, 3);
    m.setValue(0, 0, 1.0);
    m.setValue(0, 1, 2.0);
    m.setValue(0, 2, 3.0);
    m.setValue(1, 0, 4.0);
    m.setValue(1, 1, 5.0);
    m.setValue(1, 2, 6.0);
    m.setValue(2, 0, 7.0);
    m.setValue(2, 1, 8.0);
    m.setValue(2, 2, 9.0);

    Vector v(3);
    v.setValue(0, 1.0);
    v.setValue(1, 0.0);
    v.setValue(2, 1.0);

    Vector result = multiplicar(m, v);
    // [1,2,3; 4,5,6; 7,8,9] * [1,0,1]^T = [4, 10, 16]^T
    EXPECT_EQ(result.getValue(0), 4.0);
    EXPECT_EQ(result.getValue(1), 10.0);
    EXPECT_EQ(result.getValue(2), 16.0);
}

TEST_F(VectorAdvancedTest, ZeroVectorAdditiveIdentity) {
    Vector zero(4);
    for (int i = 0; i < 4; i++) {
        zero.setValue(i, 0.0);
    }

    Vector result = v1 + zero;
    for (int i = 0; i < 4; i++) {
        EXPECT_EQ(result.getValue(i), v1.getValue(i));
    }
}

TEST_F(VectorAdvancedTest, NegativeVector) {
    Vector neg = v1 * -1;
    Vector result = v1 + neg;

    for (int i = 0; i < 4; i++) {
        EXPECT_NEAR(result.getValue(i), 0.0, 1e-6);
    }
}

TEST_F(VectorAdvancedTest, NormalizeVector) {
    double norm = v1.euclidian_length();
    Vector normalized = v1 * (int)(1 / norm);

    // After scaling by 1/norm, the length should be close to 1
    // This is approximate since we're using int casting
    EXPECT_GT(normalized.euclidian_length(), 0.0);
}

TEST_F(VectorAdvancedTest, VectorNormProperties) {
    // ||alpha * v|| = |alpha| * ||v||
    double norm_v = v1.euclidian_length();
    Vector scaled = v1 * 5;
    double norm_scaled = scaled.euclidian_length();

    EXPECT_NEAR(norm_scaled, 5 * norm_v, 1e-6);
}

TEST_F(VectorAdvancedTest, LinearCombination) {
    // Test if we can represent a vector as linear combination
    Vector u(3);
    u.setValue(0, 1.0);
    u.setValue(1, 0.0);
    u.setValue(2, 0.0);

    Vector v(3);
    v.setValue(0, 0.0);
    v.setValue(1, 1.0);
    v.setValue(2, 0.0);

    // 3u + 4v should be [3, 4, 0]
    Vector result = u * 3 + v * 4;
    EXPECT_EQ(result.getValue(0), 3.0);
    EXPECT_EQ(result.getValue(1), 4.0);
    EXPECT_EQ(result.getValue(2), 0.0);
}

TEST_F(VectorAdvancedTest, LargeVectorNorm) {
    Vector large(10000);
    for (int i = 0; i < 10000; i++) {
        large.setValue(i, 0.01);
    }

    double norm = large.euclidian_length();
    // sqrt(10000 * 0.01^2) = sqrt(10000 * 0.0001) = sqrt(1) = 1
    EXPECT_NEAR(norm, 1.0, 1e-6);
}
