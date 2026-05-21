#include "utils/test.hpp"
#include "algebra_linear/vector.hpp"
#include <stdexcept>

static Vector make_v1() {
    Vector v(3);
    v.setValue(0, 1.0); v.setValue(1, 2.0); v.setValue(2, 3.0);
    return v;
}

static Vector make_v2() {
    Vector v(3);
    v.setValue(0, 4.0); v.setValue(1, 5.0); v.setValue(2, 6.0);
    return v;
}

static void test_constructor_default() {
    Vector v;
    ASSERT_EQ(v.getLength(), 0);
}

static void test_constructor_dimension() {
    Vector v(5);
    ASSERT_EQ(v.getLength(), 5);
}

static void test_get_set_value() {
    Vector v(3);
    v.setValue(1, 10.5);
    ASSERT_NEAR(v.getValue(1), 10.5, 0.0);
}

static void test_get_value_out_of_range() {
    Vector v(3);
    ASSERT_THROW(v.getValue(5), std::out_of_range);
}

static void test_set_value_out_of_range() {
    Vector v(3);
    ASSERT_THROW(v.setValue(10, 5.0), std::out_of_range);
}

static void test_copy_constructor() {
    Vector v1 = make_v1();
    Vector copy = v1;
    ASSERT_EQ(copy.getLength(), v1.getLength());
    ASSERT_VEC_NEAR(copy, v1, 0.0);
}

static void test_addition() {
    // [1,2,3] + [4,5,6] = [5,7,9]
    Vector result = make_v1() + make_v2();
    ASSERT_EQ(result.getLength(), 3);
    ASSERT_NEAR(result.getValue(0), 5.0, 0.0);
    ASSERT_NEAR(result.getValue(1), 7.0, 0.0);
    ASSERT_NEAR(result.getValue(2), 9.0, 0.0);
}

static void test_subtraction() {
    // [4,5,6] - [1,2,3] = [3,3,3]
    Vector result = make_v2() - make_v1();
    ASSERT_EQ(result.getLength(), 3);
    ASSERT_NEAR(result.getValue(0), 3.0, 0.0);
    ASSERT_NEAR(result.getValue(1), 3.0, 0.0);
    ASSERT_NEAR(result.getValue(2), 3.0, 0.0);
}

static void test_scalar_multiplication() {
    Vector result = make_v1() * 2.0;
    ASSERT_NEAR(result.getValue(0), 2.0, 0.0);
    ASSERT_NEAR(result.getValue(1), 4.0, 0.0);
    ASSERT_NEAR(result.getValue(2), 6.0, 0.0);
}

static void test_scalar_multiplication_by_zero() {
    Vector result = make_v1() * 0.0;
    ASSERT_NEAR(result.getValue(0), 0.0, 0.0);
    ASSERT_NEAR(result.getValue(1), 0.0, 0.0);
    ASSERT_NEAR(result.getValue(2), 0.0, 0.0);
}

static void test_dot_product() {
    // 1*4 + 2*5 + 3*6 = 32
    ASSERT_NEAR(make_v1().linear_product(make_v2()), 32.0, 0.0);
}

static void test_dot_product_zero_vector() {
    Vector zero(3);
    ASSERT_NEAR(zero.linear_product(make_v1()), 0.0, 0.0);
}

static void test_euclidian_length() {
    // ||(3,4,0)|| = 5
    Vector v(3);
    v.setValue(0, 3.0); v.setValue(1, 4.0); v.setValue(2, 0.0);
    ASSERT_NEAR(v.euclidian_length(), 5.0, 1e-6);
}

static void test_euclidian_length_unit() {
    Vector e1(3);
    e1.setValue(0, 1.0);
    ASSERT_NEAR(e1.euclidian_length(), 1.0, 1e-6);
}

static void test_equality() {
    Vector v1 = make_v1();
    Vector v3 = v1;
    ASSERT_TRUE(v1.equality(v3));
    ASSERT_FALSE(v1.equality(make_v2()));
}

static void test_addition_incompatible() {
    Vector v1 = make_v1();
    Vector v3(5);
    ASSERT_THROW(v1 + v3, std::invalid_argument);
}

static void test_axpy() {
    // axpy(y, alpha) = alpha * this + y: 2*[1,2,3] + [4,5,6] = [6,9,12]
    Vector result = make_v1().axpy(make_v2(), 2.0);
    ASSERT_NEAR(result.getValue(0),  6.0, 0.0);
    ASSERT_NEAR(result.getValue(1),  9.0, 0.0);
    ASSERT_NEAR(result.getValue(2), 12.0, 0.0);
}

static void test_axpy_zero_alpha() {
    // 0*v1 + v2 = v2
    Vector result = make_v1().axpy(make_v2(), 0.0);
    ASSERT_VEC_NEAR(result, make_v2(), 0.0);
}

static void test_dot_product_with_itself() {
    Vector v1 = make_v1();
    double norm = v1.euclidian_length();
    ASSERT_NEAR(v1.linear_product(v1), norm * norm, 1e-6);
}

static void test_large_vector_operations() {
    Vector large(1000);
    for (int i = 0; i < 1000; i++)
        large.setValue(i, 0.001 * i);
    Vector doubled = large * 2.0;
    for (int i = 0; i < 1000; i++)
        ASSERT_NEAR(doubled.getValue(i), 0.002 * i, 1e-10);
}

static void test_component_wise_product() {
    // [1,2,3] * [4,5,6] = [4,10,18]
    Vector result = make_v1() * make_v2();
    ASSERT_NEAR(result.getValue(0),  4.0, 0.0);
    ASSERT_NEAR(result.getValue(1), 10.0, 0.0);
    ASSERT_NEAR(result.getValue(2), 18.0, 0.0);
}

int main() {
    RUN_TEST(test_constructor_default);
    RUN_TEST(test_constructor_dimension);
    RUN_TEST(test_get_set_value);
    RUN_TEST(test_get_value_out_of_range);
    RUN_TEST(test_set_value_out_of_range);
    RUN_TEST(test_copy_constructor);
    RUN_TEST(test_addition);
    RUN_TEST(test_subtraction);
    RUN_TEST(test_scalar_multiplication);
    RUN_TEST(test_scalar_multiplication_by_zero);
    RUN_TEST(test_dot_product);
    RUN_TEST(test_dot_product_zero_vector);
    RUN_TEST(test_euclidian_length);
    RUN_TEST(test_euclidian_length_unit);
    RUN_TEST(test_equality);
    RUN_TEST(test_addition_incompatible);
    RUN_TEST(test_axpy);
    RUN_TEST(test_axpy_zero_alpha);
    RUN_TEST(test_dot_product_with_itself);
    RUN_TEST(test_large_vector_operations);
    RUN_TEST(test_component_wise_product);
    TEST_SUMMARY();
}
