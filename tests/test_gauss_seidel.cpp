#include "utils/test.hpp"
#include "algebra_linear/matriz.hpp"
#include "metodos/gauss_seidel.hpp"

// Sistema 3x3 dominante diagonal: x* = [1, 1, 1]
// A = [4 1 0]   b = [5]
//     [1 3 1]       [5]
//     [0 1 2]       [3]
static Matriz make_A() {
    Matriz A(3, 3);
    A.setValue(0, 0, 4); A.setValue(0, 1, 1); A.setValue(0, 2, 0);
    A.setValue(1, 0, 1); A.setValue(1, 1, 3); A.setValue(1, 2, 1);
    A.setValue(2, 0, 0); A.setValue(2, 1, 1); A.setValue(2, 2, 2);
    return A;
}

// gauss_seidel espera b como Matriz(n, 1)
static Matriz make_b() {
    Matriz b(3, 1);
    b.setValue(0, 0, 5); b.setValue(1, 0, 5); b.setValue(2, 0, 3);
    return b;
}

static void test_solve_solucao_conhecida() {
    Matriz x = gaussseidel::gauss_seidel(make_A(), make_b());

    ASSERT_NEAR(x.getValue(0, 0), 1.0, 1e-5);
    ASSERT_NEAR(x.getValue(1, 0), 1.0, 1e-5);
    ASSERT_NEAR(x.getValue(2, 0), 1.0, 1e-5);
}

static void test_residual() {
    Matriz A = make_A();
    Matriz b = make_b();
    Matriz x = gaussseidel::gauss_seidel(A, b);

    Matriz r = A * x - b;
    double norma_quad = 0.0;
    for (int i = 0; i < r.getRows(); i++)
        norma_quad += r.getValue(i, 0) * r.getValue(i, 0);
    ASSERT_NEAR(norma_quad, 0.0, 1e-10);
}

int main() {
    RUN_TEST(test_solve_solucao_conhecida);
    RUN_TEST(test_residual);
    TEST_SUMMARY();
}
