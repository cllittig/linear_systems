#include <gtest/gtest.h>
#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"
#include "metodos/conjugate_gradient.hpp"

class ConjugateGradientTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(ConjugateGradientTest, SimpleSystem) {
    Matriz A(2, 2);
    A.setValue(0, 0, 2.0);
    A.setValue(0, 1, 1.0);
    A.setValue(1, 0, 1.0);
    A.setValue(1, 1, 2.0);
    
    Vector b(2);
    b.setValue(0, 5.0);
    b.setValue(1, 4.0);
    
    Vector x = conjugate_gradient::solve(A, b, 1e-10, 100);
    
    EXPECT_NEAR(x.getValue(0), 2.0, 1e-6);
    EXPECT_NEAR(x.getValue(1), 1.0, 1e-6);
}

TEST_F(ConjugateGradientTest, ThreeByThreeSystem) {
    Matriz A(3, 3);
    A.setValue(0, 0, 4.0);
    A.setValue(0, 1, 12.0);
    A.setValue(0, 2, -16.0);
    A.setValue(1, 0, 12.0);
    A.setValue(1, 1, 37.0);
    A.setValue(1, 2, -43.0);
    A.setValue(2, 0, -16.0);
    A.setValue(2, 1, -43.0);
    A.setValue(2, 2, 98.0);
    
    Vector b(3);
    b.setValue(0, 0.0);
    b.setValue(1, 6.0);
    b.setValue(2, 39.0);
    
    Vector x = conjugate_gradient::solve(A, b, 1e-10, 100);
    
    Vector Ax = multiplicar(A, x);
    for (int i = 0; i < 3; i++) {
        EXPECT_NEAR(b.getValue(i), Ax.getValue(i), 1e-6)
            << "Verification failed at index " << i;
    }
}

TEST_F(ConjugateGradientTest, IdentityMatrix) {
    Matriz I = Matriz::identity(2, 2);
    
    Vector b(2);
    b.setValue(0, 3.0);
    b.setValue(1, 7.0);
    
    Vector x = conjugate_gradient::solve(I, b, 1e-10, 100);
    
    EXPECT_NEAR(x.getValue(0), 3.0, 1e-9);
    EXPECT_NEAR(x.getValue(1), 7.0, 1e-9);
}

TEST_F(ConjugateGradientTest, ScaledSystem) {
    Matriz A(2, 2);
    A.setValue(0, 0, 3.0);
    A.setValue(0, 1, 1.0);
    A.setValue(1, 0, 1.0);
    A.setValue(1, 1, 3.0);
    
    Vector b(2);
    b.setValue(0, 7.0);
    b.setValue(1, 5.0);
    
    Vector x = conjugate_gradient::solve(A, b, 1e-10, 100);
    
    EXPECT_NEAR(x.getValue(0), 2.0, 1e-6);
    EXPECT_NEAR(x.getValue(1), 1.0, 1e-6);
    
    Vector Ax = multiplicar(A, x);
    EXPECT_NEAR(b.getValue(0), Ax.getValue(0), 1e-6);
    EXPECT_NEAR(b.getValue(1), Ax.getValue(1), 1e-6);
}

TEST_F(ConjugateGradientTest, NumericalStability) {
    Matriz A(4, 4);
    A.setValue(0, 0, 10.0);
    A.setValue(0, 1, 1.0);
    A.setValue(0, 2, 0.0);
    A.setValue(0, 3, 0.0);
    A.setValue(1, 0, 1.0);
    A.setValue(1, 1, 10.0);
    A.setValue(1, 2, 1.0);
    A.setValue(1, 3, 0.0);
    A.setValue(2, 0, 0.0);
    A.setValue(2, 1, 1.0);
    A.setValue(2, 2, 10.0);
    A.setValue(2, 3, 1.0);
    A.setValue(3, 0, 0.0);
    A.setValue(3, 1, 0.0);
    A.setValue(3, 2, 1.0);
    A.setValue(3, 3, 10.0);
    
    Vector b(4);
    b.setValue(0, 11.0);
    b.setValue(1, 12.0);
    b.setValue(2, 12.0);
    b.setValue(3, 11.0);
    
    Vector x = conjugate_gradient::solve(A, b, 1e-10, 100);
    
    Vector Ax = multiplicar(A, x);
    for (int i = 0; i < 4; i++) {
        EXPECT_NEAR(b.getValue(i), Ax.getValue(i), 1e-5)
            << "Numerical stability check failed at index " << i;
    }
}

TEST_F(ConjugateGradientTest, LargerSystem) {
    int n = 5;
    Matriz A(n, n);
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                A.setValue(i, j, 10.0);
            } else if (i == j + 1 || i == j - 1) {
                A.setValue(i, j, 1.0);
            } else {
                A.setValue(i, j, 0.0);
            }
        }
    }
    
    Vector b(n);
    for (int i = 0; i < n; i++) {
        b.setValue(i, 11.0);
    }
    
    Vector x = conjugate_gradient::solve(A, b, 1e-10, 200);
    
    Vector Ax = multiplicar(A, x);
    for (int i = 0; i < n; i++) {
        EXPECT_NEAR(b.getValue(i), Ax.getValue(i), 1e-5)
            << "Larger system check failed at index " << i;
    }
}
