#include "utils/test.hpp"
#include "algebra_linear/matriz.hpp"
#include "metodos/lu.hpp"

// Sistema 3x3: x* = [1, 1, 1]
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

// lu::lu espera b como Matriz(1, n) — vetor linha
static Matriz make_b() {
    Matriz b(1, 3);
    b.setValue(0, 0, 5); b.setValue(0, 1, 5); b.setValue(0, 2, 3);
    return b;
}

static void test_solve_solucao_conhecida() {
    Matriz x = lu::lu(make_A(), make_b());

    ASSERT_NEAR(x.getValue(0, 0), 1.0, 1e-10);
    ASSERT_NEAR(x.getValue(0, 1), 1.0, 1e-10);
    ASSERT_NEAR(x.getValue(0, 2), 1.0, 1e-10);
}

static void test_residual() {
    Matriz A = make_A();
    Matriz b = make_b();
    Matriz x = lu::lu(A, b);

    // r = Ax - b  (aqui x e b são vetores linha Matriz(1,n))
    Matriz r = x * A - b;
    double norma_quad = 0.0;
    for (int j = 0; j < r.getColumns(); j++)
        norma_quad += r.getValue(0, j) * r.getValue(0, j);
    ASSERT_NEAR(norma_quad, 0.0, 1e-10);
}

int main() {
    RUN_TEST(test_solve_solucao_conhecida);
    RUN_TEST(test_residual);
    TEST_SUMMARY();
}
