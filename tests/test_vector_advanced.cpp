#include "utils/test.hpp"
#include "algebra_linear/vector.hpp"
#include "algebra_linear/matriz.hpp"
#include <cmath>

static Vector make_v1() {
    Vector v(4);
    v.setValue(0, 1.0); v.setValue(1, 2.0); v.setValue(2, 3.0); v.setValue(3, 4.0);
    return v;
}

static Vector make_v2() {
    Vector v(4);
    v.setValue(0, 2.0); v.setValue(1, 0.0); v.setValue(2, -1.0); v.setValue(3, 2.0);
    return v;
}

static void test_orthogonal_vectors() {
    Vector a(3), b(3);
    a.setValue(0, 1.0);
    b.setValue(1, 1.0);
    ASSERT_NEAR(a.linear_product(b), 0.0, 0.0);
}

static void test_unit_vector() {
    double s = 1.0 / std::sqrt(3.0);
    Vector u(3);
    u.setValue(0, s); u.setValue(1, s); u.setValue(2, s);
    ASSERT_NEAR(u.euclidian_length(), 1.0, 1e-6);
}

static void test_cauchy_schwarz() {
    // |<u,v>| <= ||u|| * ||v||
    Vector v1 = make_v1(), v2 = make_v2();
    double dot = std::fabs(v1.linear_product(v2));
    ASSERT_LE(dot, v1.euclidian_length() * v2.euclidian_length() + 1e-6);
}

static void test_triangle_inequality() {
    // ||u + v|| <= ||u|| + ||v||
    Vector v1 = make_v1(), v2 = make_v2();
    ASSERT_LE((v1 + v2).euclidian_length(), v1.euclidian_length() + v2.euclidian_length() + 1e-6);
}

static void test_pythagorean_theorem() {
    // 3-4-5 triangle
    Vector a(2), b(2);
    a.setValue(0, 3.0); b.setValue(1, 4.0);
    ASSERT_NEAR((a + b).euclidian_length(), 5.0, 1e-6);
    ASSERT_NEAR(a.euclidian_length() * a.euclidian_length() +
                b.euclidian_length() * b.euclidian_length(),
                25.0, 1e-6);
}

static void test_vector_projection() {
    // proj_v2(v1) = (v1·v2 / v2·v2) * v2
    Vector v1 = make_v1(), v2 = make_v2();
    double scalar = v1.linear_product(v2) / v2.linear_product(v2);
    Vector proj = v2 * scalar;
    double expected_length = std::fabs(v1.linear_product(v2)) / v2.euclidian_length();
    ASSERT_NEAR(proj.euclidian_length(), expected_length, 1e-6);
}

static void test_addition_associative() {
    Vector v3(4);
    for (int i = 0; i < 4; i++) v3.setValue(i, 1.0);
    Vector v1 = make_v1(), v2 = make_v2();
    Vector r1 = (v1 + v2) + v3;
    Vector r2 = v1 + (v2 + v3);
    ASSERT_VEC_NEAR(r1, r2, 1e-6);
}

static void test_scalar_distributive() {
    // k*(u+v) = k*u + k*v
    Vector v1 = make_v1(), v2 = make_v2();
    ASSERT_VEC_NEAR((v1 + v2) * 3.0, v1 * 3.0 + v2 * 3.0, 1e-6);
}

static void test_matrix_vector_multiplication() {
    // [1,2,3; 4,5,6; 7,8,9] * [1,0,1]^T = [4, 10, 16]
    Matriz m(3, 3);
    m.setValue(0,0,1); m.setValue(0,1,2); m.setValue(0,2,3);
    m.setValue(1,0,4); m.setValue(1,1,5); m.setValue(1,2,6);
    m.setValue(2,0,7); m.setValue(2,1,8); m.setValue(2,2,9);
    Vector v(3);
    v.setValue(0, 1.0); v.setValue(2, 1.0);
    Vector result = multiplicar(m, v);
    ASSERT_NEAR(result.getValue(0),  4.0, 0.0);
    ASSERT_NEAR(result.getValue(1), 10.0, 0.0);
    ASSERT_NEAR(result.getValue(2), 16.0, 0.0);
}

static void test_zero_additive_identity() {
    Vector v1 = make_v1();
    Vector zero(4);
    ASSERT_VEC_NEAR(v1 + zero, v1, 0.0);
}

static void test_negative_vector() {
    Vector v1 = make_v1();
    Vector neg = v1 * -1.0;
    Vector zero(4);
    ASSERT_VEC_NEAR(v1 + neg, zero, 1e-6);
}

static void test_normalize_vector() {
    Vector v1 = make_v1();
    Vector normalized = v1 * (1.0 / v1.euclidian_length());
    ASSERT_NEAR(normalized.euclidian_length(), 1.0, 1e-6);
}

static void test_norm_scaling() {
    // ||alpha*v|| = |alpha| * ||v||
    Vector v1 = make_v1();
    ASSERT_NEAR((v1 * 5.0).euclidian_length(), 5.0 * v1.euclidian_length(), 1e-6);
}

static void test_linear_combination() {
    // 3*e1 + 4*e2 = [3, 4, 0]
    Vector u(3), v(3);
    u.setValue(0, 1.0); v.setValue(1, 1.0);
    Vector result = u * 3.0 + v * 4.0;
    ASSERT_NEAR(result.getValue(0), 3.0, 0.0);
    ASSERT_NEAR(result.getValue(1), 4.0, 0.0);
    ASSERT_NEAR(result.getValue(2), 0.0, 0.0);
}

static void test_large_vector_norm() {
    // sqrt(10000 * 0.01^2) = 1
    Vector large(10000);
    for (int i = 0; i < 10000; i++)
        large.setValue(i, 0.01);
    ASSERT_NEAR(large.euclidian_length(), 1.0, 1e-6);
}

int main() {
    RUN_TEST(test_orthogonal_vectors);
    RUN_TEST(test_unit_vector);
    RUN_TEST(test_cauchy_schwarz);
    RUN_TEST(test_triangle_inequality);
    RUN_TEST(test_pythagorean_theorem);
    RUN_TEST(test_vector_projection);
    RUN_TEST(test_addition_associative);
    RUN_TEST(test_scalar_distributive);
    RUN_TEST(test_matrix_vector_multiplication);
    RUN_TEST(test_zero_additive_identity);
    RUN_TEST(test_negative_vector);
    RUN_TEST(test_normalize_vector);
    RUN_TEST(test_norm_scaling);
    RUN_TEST(test_linear_combination);
    RUN_TEST(test_large_vector_norm);
    TEST_SUMMARY();
}
