#include "utils/test.hpp"
#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"
#include "metodos/cholesky.hpp"

// Sistema 3x3 SPD: x* = [1, 1, 1]
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

static Vector make_b() {
    Vector b(3);
    b.setValue(0, 5); b.setValue(1, 5); b.setValue(2, 3);
    return b;
}

static void test_solve_solucao_conhecida() {
    Vector x = cholesky::solve(make_A(), make_b());

    Vector esperado(3);
    esperado.setValue(0, 1.0); esperado.setValue(1, 1.0); esperado.setValue(2, 1.0);

    ASSERT_VEC_NEAR(x, esperado, 1e-10);
}

static void test_residual() {
    Matriz A = make_A();
    Vector b = make_b();
    Vector x = cholesky::solve(A, b);

    Vector r = multiplicar(A, x) - b;
    ASSERT_NEAR(r.euclidian_length(), 0.0, 1e-10);
}

int main() {
    RUN_TEST(test_solve_solucao_conhecida);
    RUN_TEST(test_residual);
    TEST_SUMMARY();
}
