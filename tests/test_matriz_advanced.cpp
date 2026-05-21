#include "utils/test.hpp"
#include "algebra_linear/matriz.hpp"

static Matriz make_sym() {
    Matriz m(3, 3);
    m.setValue(0, 0, 2.0); m.setValue(0, 1, 1.0); m.setValue(0, 2, 0.0);
    m.setValue(1, 0, 1.0); m.setValue(1, 1, 2.0); m.setValue(1, 2, 1.0);
    m.setValue(2, 0, 0.0); m.setValue(2, 1, 1.0); m.setValue(2, 2, 2.0);
    return m;
}

static Matriz make_nonsym() {
    Matriz m(3, 3);
    m.setValue(0, 0, 1.0); m.setValue(0, 1, 2.0); m.setValue(0, 2, 3.0);
    m.setValue(1, 0, 4.0); m.setValue(1, 1, 5.0); m.setValue(1, 2, 6.0);
    m.setValue(2, 0, 7.0); m.setValue(2, 1, 8.0); m.setValue(2, 2, 9.0);
    return m;
}

static void test_is_simetric_true() {
    ASSERT_TRUE(make_sym().isSimetric());
}

static void test_is_simetric_false() {
    ASSERT_FALSE(make_nonsym().isSimetric());
}

static void test_determinant_2x2() {
    // det([[3,8],[4,6]]) = 3*6 - 8*4 = -14
    Matriz m(2, 2);
    m.setValue(0, 0, 3.0); m.setValue(0, 1, 8.0);
    m.setValue(1, 0, 4.0); m.setValue(1, 1, 6.0);
    ASSERT_NEAR(m.determinant(), -14.0, 1e-6);
}

static void test_determinant_3x3() {
    // det = -24 + 40 - 15 = 1
    Matriz m(3, 3);
    m.setValue(0, 0, 1.0); m.setValue(0, 1, 2.0); m.setValue(0, 2, 3.0);
    m.setValue(1, 0, 0.0); m.setValue(1, 1, 1.0); m.setValue(1, 2, 4.0);
    m.setValue(2, 0, 5.0); m.setValue(2, 1, 6.0); m.setValue(2, 2, 0.0);
    ASSERT_NEAR(m.determinant(), 1.0, 1e-6);
}

static void test_transpose_of_transpose() {
    Matriz m(2, 3);
    m.setValue(0, 0, 1.0); m.setValue(0, 1, 2.0); m.setValue(0, 2, 3.0);
    m.setValue(1, 0, 4.0); m.setValue(1, 1, 5.0); m.setValue(1, 2, 6.0);
    ASSERT_MAT_NEAR(m.transpose().transpose(), m, 0.0);
}

static void test_chain_multiplication_associative() {
    // (A*B)*C == A*(B*C)
    Matriz a(2, 2), b(2, 2), c(2, 2);
    a.setValue(0,0,1); a.setValue(0,1,2); a.setValue(1,0,3); a.setValue(1,1,4);
    b.setValue(0,0,2); b.setValue(0,1,0); b.setValue(1,0,1); b.setValue(1,1,2);
    c.setValue(0,0,1); c.setValue(0,1,1); c.setValue(1,0,1); c.setValue(1,1,1);
    ASSERT_MAT_NEAR((a * b) * c, a * (b * c), 1e-6);
}

static void test_addition_associative() {
    Matriz a(2,2), b(2,2), c(2,2);
    a.setValue(0,0,1); a.setValue(0,1,2); a.setValue(1,0,3); a.setValue(1,1,4);
    b.setValue(0,0,5); b.setValue(0,1,6); b.setValue(1,0,7); b.setValue(1,1,8);
    c.setValue(0,0,1); c.setValue(0,1,1); c.setValue(1,0,1); c.setValue(1,1,1);
    ASSERT_MAT_NEAR((a + b) + c, a + (b + c), 0.0);
}

static void test_distributive_property() {
    // A*(B+C) == A*B + A*C
    Matriz a(2,2), b(2,2), c(2,2);
    a.setValue(0,0,1); a.setValue(0,1,2); a.setValue(1,0,3); a.setValue(1,1,4);
    b.setValue(0,0,1); b.setValue(0,1,0); b.setValue(1,0,0); b.setValue(1,1,1);
    c.setValue(0,0,0); c.setValue(0,1,1); c.setValue(1,0,1); c.setValue(1,1,0);
    ASSERT_MAT_NEAR(a * (b + c), (a * b) + (a * c), 1e-6);
}

static void test_identity_properties() {
    Matriz I = Matriz::identity(5, 5);
    Matriz m(5, 5);
    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
            m.setValue(i, j, i * 5.0 + j + 1.0);
    ASSERT_MAT_NEAR(m * I, m, 1e-6);
}

static void test_large_matrix_scalar() {
    Matriz large(100, 100);
    for (int i = 0; i < 100; i++)
        for (int j = 0; j < 100; j++)
            large.setValue(i, j, 0.01 * (i + j));
    Matriz scaled = large * 2.0;
    Matriz expected(100, 100);
    for (int i = 0; i < 100; i++)
        for (int j = 0; j < 100; j++)
            expected.setValue(i, j, large.getValue(i, j) * 2.0);
    ASSERT_MAT_NEAR(scaled, expected, 0.0);
}

int main() {
    RUN_TEST(test_is_simetric_true);
    RUN_TEST(test_is_simetric_false);
    RUN_TEST(test_determinant_2x2);
    RUN_TEST(test_determinant_3x3);
    RUN_TEST(test_transpose_of_transpose);
    RUN_TEST(test_chain_multiplication_associative);
    RUN_TEST(test_addition_associative);
    RUN_TEST(test_distributive_property);
    RUN_TEST(test_identity_properties);
    RUN_TEST(test_large_matrix_scalar);
    TEST_SUMMARY();
}
