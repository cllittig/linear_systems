#include <gtest/gtest.h>
#include <sstream>
#include "algebra_linear/matriz.hpp"

class MatrizAdvancedTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(MatrizAdvancedTest, Determinant2x2) {
    Matriz m(2, 2);
    m.setValue(0, 0, 4.0);
    m.setValue(0, 1, 7.0);
    m.setValue(1, 0, 2.0);
    m.setValue(1, 1, 6.0);
    
    double det = m.determinant();
    // det = 4*6 - 7*2 = 24 - 14 = 10
    EXPECT_DOUBLE_EQ(det, 10.0);
}

TEST_F(MatrizAdvancedTest, Determinant3x3) {
    Matriz m(3, 3);
    m.setValue(0, 0, 1.0);
    m.setValue(0, 1, 2.0);
    m.setValue(0, 2, 3.0);
    m.setValue(1, 0, 0.0);
    m.setValue(1, 1, 4.0);
    m.setValue(1, 2, 5.0);
    m.setValue(2, 0, 1.0);
    m.setValue(2, 1, 0.0);
    m.setValue(2, 2, 6.0);
    
    double det = m.determinant();
    // Usando Sarrus: 1*4*6 + 2*5*1 + 3*0*0 - 3*4*1 - 1*5*0 - 2*0*6
    // = 24 + 10 + 0 - 12 - 0 - 0 = 22
    EXPECT_DOUBLE_EQ(det, 22.0);
}

TEST_F(MatrizAdvancedTest, Determinant1x1) {
    Matriz m(1, 1);
    m.setValue(0, 0, 5.0);
    
    double det = m.determinant();
    EXPECT_DOUBLE_EQ(det, 5.0);
}

TEST_F(MatrizAdvancedTest, DeterminantZero) {
    Matriz m(2, 2);
    m.setValue(0, 0, 1.0);
    m.setValue(0, 1, 2.0);
    m.setValue(1, 0, 2.0);
    m.setValue(1, 1, 4.0);  // Linhas proporcionais - det = 0
    
    double det = m.determinant();
    EXPECT_NEAR(det, 0.0, 1e-9);
}

TEST_F(MatrizAdvancedTest, Inverse2x2) {
    Matriz m(2, 2);
    m.setValue(0, 0, 4.0);
    m.setValue(0, 1, 7.0);
    m.setValue(1, 0, 2.0);
    m.setValue(1, 1, 6.0);
    
    Matriz inv = m.inverse();
    
    // Verificar: M * M^-1 = I
    Matriz product = m * inv;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            if (i == j) {
                EXPECT_NEAR(product.getValue(i, j), 1.0, 1e-9);
            } else {
                EXPECT_NEAR(product.getValue(i, j), 0.0, 1e-9);
            }
        }
    }
}

TEST_F(MatrizAdvancedTest, InverseIdentity) {
    Matriz I = Matriz::identity(3, 3);
    Matriz invI = I.inverse();
    
    // Inversa da identidade é a própria identidade
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            EXPECT_DOUBLE_EQ(invI.getValue(i, j), I.getValue(i, j));
        }
    }
}

TEST_F(MatrizAdvancedTest, InverseVerification3x3) {
    Matriz m(3, 3);
    m.setValue(0, 0, 1.0);
    m.setValue(0, 1, 2.0);
    m.setValue(0, 2, 0.0);
    m.setValue(1, 0, 0.0);
    m.setValue(1, 1, 1.0);
    m.setValue(1, 2, 3.0);
    m.setValue(2, 0, 4.0);
    m.setValue(2, 1, -1.0);
    m.setValue(2, 2, 1.0);
    
    Matriz inv = m.inverse();
    Matriz identity = m * inv;
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (i == j) {
                EXPECT_NEAR(identity.getValue(i, j), 1.0, 1e-8);
            } else {
                EXPECT_NEAR(identity.getValue(i, j), 0.0, 1e-8);
            }
        }
    }
}

TEST_F(MatrizAdvancedTest, TrocarLinhas) {
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
    
    m.trocarLinhas(0, 2);
    
    EXPECT_DOUBLE_EQ(m.getValue(0, 0), 7.0);
    EXPECT_DOUBLE_EQ(m.getValue(0, 1), 8.0);
    EXPECT_DOUBLE_EQ(m.getValue(0, 2), 9.0);
    EXPECT_DOUBLE_EQ(m.getValue(2, 0), 1.0);
    EXPECT_DOUBLE_EQ(m.getValue(2, 1), 2.0);
    EXPECT_DOUBLE_EQ(m.getValue(2, 2), 3.0);
    EXPECT_DOUBLE_EQ(m.getValue(1, 0), 4.0);
}

TEST_F(MatrizAdvancedTest, TrocarMesmaLinha) {
    Matriz m(2, 2);
    m.setValue(0, 0, 1.0);
    m.setValue(0, 1, 2.0);
    m.setValue(1, 0, 3.0);
    m.setValue(1, 1, 4.0);
    
    m.trocarLinhas(0, 0);
    
    EXPECT_DOUBLE_EQ(m.getValue(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(m.getValue(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(m.getValue(1, 0), 3.0);
    EXPECT_DOUBLE_EQ(m.getValue(1, 1), 4.0);
}

TEST_F(MatrizAdvancedTest, Clonar) {
    Matriz m1(2, 3);
    m1.setValue(0, 0, 1.5);
    m1.setValue(0, 1, 2.5);
    m1.setValue(0, 2, 3.5);
    m1.setValue(1, 0, 4.5);
    m1.setValue(1, 1, 5.5);
    m1.setValue(1, 2, 6.5);
    
    Matriz m2 = m1.clonar();
    
    EXPECT_EQ(m2.getRows(), m1.getRows());
    EXPECT_EQ(m2.getColumns(), m1.getColumns());
    
    for (int i = 0; i < m1.getRows(); i++) {
        for (int j = 0; j < m1.getColumns(); j++) {
            EXPECT_DOUBLE_EQ(m2.getValue(i, j), m1.getValue(i, j));
        }
    }
    
    // Modificar clone não afeta original
    m2.setValue(0, 0, 999.0);
    EXPECT_DOUBLE_EQ(m1.getValue(0, 0), 1.5);
    EXPECT_DOUBLE_EQ(m2.getValue(0, 0), 999.0);
}

TEST_F(MatrizAdvancedTest, IdentityInstanceMethod) {
    Matriz m(3, 3);
    Matriz id = m.identity();
    
    EXPECT_EQ(id.getRows(), 3);
    EXPECT_EQ(id.getColumns(), 3);
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (i == j) {
                EXPECT_DOUBLE_EQ(id.getValue(i, j), 1.0);
            } else {
                EXPECT_DOUBLE_EQ(id.getValue(i, j), 0.0);
            }
        }
    }
}

TEST_F(MatrizAdvancedTest, Print) {
    Matriz m(2, 2);
    m.setValue(0, 0, 1.5);
    m.setValue(0, 1, 2.5);
    m.setValue(1, 0, 3.5);
    m.setValue(1, 1, 4.5);
    
    // Test that print() doesn't throw an exception
    EXPECT_NO_THROW(m.print());
}

TEST_F(MatrizAdvancedTest, MultiplicationByZero) {
    Matriz m(2, 2);
    m.setValue(0, 0, 5.0);
    m.setValue(0, 1, 6.0);
    m.setValue(1, 0, 7.0);
    m.setValue(1, 1, 8.0);
    
    Matriz result = m * 0.0;
    
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            EXPECT_DOUBLE_EQ(result.getValue(i, j), 0.0);
        }
    }
}

TEST_F(MatrizAdvancedTest, DeterminantNotSquare) {
    Matriz m(2, 3);
    m.setValue(0, 0, 1.0);
    m.setValue(0, 1, 2.0);
    m.setValue(0, 2, 3.0);
    m.setValue(1, 0, 4.0);
    m.setValue(1, 1, 5.0);
    m.setValue(1, 2, 6.0);
    
    // Deve retornar 0 ou lançar exceção
    try {
        double det = m.determinant();
        EXPECT_EQ(det, 0.0);
    } catch (const std::exception& e) {
        EXPECT_TRUE(true);
    }
}
