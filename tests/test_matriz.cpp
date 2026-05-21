#include "utils/test.hpp"
#include "algebra_linear/matriz.hpp"
#include <stdexcept>

static Matriz make_m1() {
    Matriz m(2, 2);
    m.setValue(0, 0, 1.0); m.setValue(0, 1, 2.0);
    m.setValue(1, 0, 3.0); m.setValue(1, 1, 4.0);
    return m;
}

static Matriz make_m2() {
    Matriz m(2, 2);
    m.setValue(0, 0, 5.0); m.setValue(0, 1, 6.0);
    m.setValue(1, 0, 7.0); m.setValue(1, 1, 8.0);
    return m;
}

static void test_constructor_default() {
    Matriz m;
    ASSERT_EQ(m.getRows(), 0);
    ASSERT_EQ(m.getColumns(), 0);
}

static void test_constructor_dimensions() {
    Matriz m(3, 4);
    ASSERT_EQ(m.getRows(), 3);
    ASSERT_EQ(m.getColumns(), 4);
}

static void test_get_set_value() {
    Matriz m(2, 2);
    m.setValue(0, 0, 5.5);
    ASSERT_NEAR(m.getValue(0, 0), 5.5, 0.0);
}

static void test_get_value_out_of_range() {
    Matriz m(2, 2);
    ASSERT_THROW(m.getValue(2, 2), std::out_of_range);
}

static void test_set_value_out_of_range() {
    Matriz m(2, 2);
    ASSERT_THROW(m.setValue(5, 5, 10.0), std::out_of_range);
}

static void test_copy_constructor() {
    Matriz m1 = make_m1();
    Matriz copy = m1;
    ASSERT_EQ(copy.getRows(), m1.getRows());
    ASSERT_EQ(copy.getColumns(), m1.getColumns());
    ASSERT_NEAR(copy.getValue(0, 0), m1.getValue(0, 0), 0.0);
}

static void test_addition() {
    // [[1,2],[3,4]] + [[5,6],[7,8]] = [[6,8],[10,12]]
    Matriz result = make_m1() + make_m2();
    ASSERT_NEAR(result.getValue(0, 0),  6.0, 0.0);
    ASSERT_NEAR(result.getValue(0, 1),  8.0, 0.0);
    ASSERT_NEAR(result.getValue(1, 0), 10.0, 0.0);
    ASSERT_NEAR(result.getValue(1, 1), 12.0, 0.0);
}

static void test_subtraction() {
    // [[5,6],[7,8]] - [[1,2],[3,4]] = [[4,4],[4,4]]
    Matriz result = make_m2() - make_m1();
    ASSERT_NEAR(result.getValue(0, 0), 4.0, 0.0);
    ASSERT_NEAR(result.getValue(0, 1), 4.0, 0.0);
    ASSERT_NEAR(result.getValue(1, 0), 4.0, 0.0);
    ASSERT_NEAR(result.getValue(1, 1), 4.0, 0.0);
}

static void test_scalar_multiplication() {
    Matriz result = make_m1() * 2.0;
    ASSERT_NEAR(result.getValue(0, 0), 2.0, 0.0);
    ASSERT_NEAR(result.getValue(0, 1), 4.0, 0.0);
    ASSERT_NEAR(result.getValue(1, 0), 6.0, 0.0);
    ASSERT_NEAR(result.getValue(1, 1), 8.0, 0.0);
}

static void test_matrix_multiplication() {
    // [[1,2],[3,4]] * [[5,6],[7,8]] = [[19,22],[43,50]]
    Matriz result = make_m1() * make_m2();
    ASSERT_NEAR(result.getValue(0, 0), 19.0, 1e-12);
    ASSERT_NEAR(result.getValue(0, 1), 22.0, 1e-12);
    ASSERT_NEAR(result.getValue(1, 0), 43.0, 1e-12);
    ASSERT_NEAR(result.getValue(1, 1), 50.0, 1e-12);
}

static void test_is_square() {
    ASSERT_TRUE(Matriz(3, 3).isSquare());
    ASSERT_FALSE(Matriz(2, 3).isSquare());
}

static void test_transpose() {
    Matriz m(2, 3);
    m.setValue(0, 0, 1.0); m.setValue(0, 1, 2.0); m.setValue(0, 2, 3.0);
    m.setValue(1, 0, 4.0); m.setValue(1, 1, 5.0); m.setValue(1, 2, 6.0);

    Matriz t = m.transpose();
    ASSERT_EQ(t.getRows(), 3);
    ASSERT_EQ(t.getColumns(), 2);
    ASSERT_NEAR(t.getValue(0, 0), 1.0, 0.0);
    ASSERT_NEAR(t.getValue(1, 0), 2.0, 0.0);
    ASSERT_NEAR(t.getValue(2, 0), 3.0, 0.0);
}

static void test_clone() {
    Matriz m1 = make_m1();
    Matriz clone = m1.clonar();
    ASSERT_EQ(clone.getRows(), m1.getRows());
    ASSERT_EQ(clone.getColumns(), m1.getColumns());
    ASSERT_MAT_NEAR(clone, m1, 0.0);
}

static void test_identity() {
    Matriz I = Matriz::identity(3, 3);
    ASSERT_EQ(I.getRows(), 3);
    ASSERT_EQ(I.getColumns(), 3);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            ASSERT_NEAR(I.getValue(i, j), (i == j) ? 1.0 : 0.0, 0.0);
}

static void test_trocar_linhas() {
    Matriz m(3, 2);
    m.setValue(0, 0, 1.0);
    m.setValue(1, 0, 2.0);
    m.setValue(2, 0, 3.0);
    m.trocarLinhas(0, 2);
    ASSERT_NEAR(m.getValue(0, 0), 3.0, 0.0);
    ASSERT_NEAR(m.getValue(2, 0), 1.0, 0.0);
}

static void test_addition_incompatible() {
    Matriz m1 = make_m1();
    Matriz other(3, 3);
    ASSERT_THROW(m1 + other, std::invalid_argument);
}

static void test_subtraction_incompatible() {
    Matriz m1 = make_m1();
    Matriz other(3, 3);
    ASSERT_THROW(m1 - other, std::invalid_argument);
}

static void test_multiplication_incompatible() {
    Matriz m1 = make_m1();
    Matriz other(3, 3);
    ASSERT_THROW(m1 * other, std::invalid_argument);
}

int main() {
    RUN_TEST(test_constructor_default);
    RUN_TEST(test_constructor_dimensions);
    RUN_TEST(test_get_set_value);
    RUN_TEST(test_get_value_out_of_range);
    RUN_TEST(test_set_value_out_of_range);
    RUN_TEST(test_copy_constructor);
    RUN_TEST(test_addition);
    RUN_TEST(test_subtraction);
    RUN_TEST(test_scalar_multiplication);
    RUN_TEST(test_matrix_multiplication);
    RUN_TEST(test_is_square);
    RUN_TEST(test_transpose);
    RUN_TEST(test_clone);
    RUN_TEST(test_identity);
    RUN_TEST(test_trocar_linhas);
    RUN_TEST(test_addition_incompatible);
    RUN_TEST(test_subtraction_incompatible);
    RUN_TEST(test_multiplication_incompatible);
    TEST_SUMMARY();
}
