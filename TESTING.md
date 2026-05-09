# 🧪 Guia de Testes Unitários

## Visão Geral

Este projeto implementa uma suite completa de **32 testes unitários** usando **Google Test**, um framework de testes moderno e consolidado para C++.

## Configuração Automática

O CMakeLists.txt faz:
1. Download automático do Google Test v1.12.1
2. Compilação local da biblioteca
3. Linkagem dos testes automaticamente
4. Geração de executável `programa_testes`

Não é necessário instalar Google Test no sistema!

## Executando os Testes

### Executar todos os testes
```bash
cd build
./programa_testes
```

**Saída esperada:**
```
[==========] Running 32 tests from 3 test suites.
[----------] 12 tests from MatrizTest
...
[==========] 32 tests from 3 test suites ran. (X ms total)
[  PASSED  ] 32 tests.
```

### Opções úteis

```bash
# Mostrar tudo
./programa_testes --verbose

# Executor um teste específico
./programa_testes --gtest_filter=MatrizTest.Addition

# Executar todos de uma suite
./programa_testes --gtest_filter=MatrizTest.*

# Listar testes sem executar
./programa_testes --gtest_list_tests

# Repetir testes N vezes
./programa_testes --gtest_repeat=10

# Parar no primeiro erro
./programa_testes --gtest_break_on_failure
```

## Estrutura dos Testes

### 1️⃣ Teste de Matriz (12 testes)

```cpp
TEST_F(MatrizTest, Addition) {
    Matriz m1(2, 2), m2(2, 2);
    // Setup
    m1.setValue(0, 0, 1.0);
    // Teste
    Matriz m3 = m1 + m2;
    EXPECT_DOUBLE_EQ(m3.getValue(0, 0), 2.0);
}
```

**Testes inclusos:**
- `ConstructorDefault` - Verifica construtor padrão
- `ConstructorWithDimensions` - Verifica construtor parametrizado
- `SetAndGetValue` - Operações básicas de acesso
- `CopyConstructor` - Cópia profunda de matrizes
- `Addition` - Adição de matrizes
- `Subtraction` - Subtração de matrizes
- `ScalarMultiplication` - Produto por escalar
- `MatrixMultiplication` - Produto matricial
- `Transpose` - Matriz transposta
- `IsSquare` - Verificação de quadratura
- `IsSymmetric` - Verificação de simetria
- `IdentityMatrix` - Matriz identidade

### 2️⃣ Teste de Vector (14 testes)

```cpp
TEST_F(VectorTest, LinearProduct) {
    Vector v1(3), v2(3);
    v1.setValue(0, 1.0);
    v2.setValue(0, 4.0);
    double result = v1.linear_product(v2);
    EXPECT_DOUBLE_EQ(result, 32.0);
}
```

**Testes inclusos:**
- `ConstructorDefault` e `ConstructorWithLength`
- `SetAndGetValue`
- `CopyConstructor`
- `Addition` e `Subtraction`
- `ScalarMultiplication`
- `LinearProduct` (produto escalar)
- `ElementWiseProduct` (produto elemento a elemento)
- `EuclidianLength` (norma Euclidiana)
- `AXPYOperation` (operação AXPY)
- `Equality` (comparação)
- `MatrixVectorMultiplication`

### 3️⃣ Teste de Gradiente Conjugado (6 testes)

```cpp
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
}
```

**Testes inclusos:**
- `SimpleSystem` - Sistema 2x2 básico
- `ThreeByThreeSystem` - Sistema 3x3 SPD
- `IdentityMatrix` - Sistema com matriz identidade
- `ScaledSystem` - Sistema escalado
- `NumericalStability` - Sistema 4x4 para validar estabilidade
- `LargerSystem` - Sistema 5x5 com banda de dois diagonais

## Casos de Teste por Categoria

### ✅ Operações Básicas
- Set/Get de valores
- Construtores e destrutores
- Cópia profunda

### ✅ Aritmética Linear
- Adição/Subtração
- Multiplicação escalar
- Multiplicação matricial/vetorial

### ✅ Propriedades Matriciais
- Simetria
- Quadratura  
- Identidade
- Norms

### ✅ Métodos Numéricos
- Convergência do Gradiente Conjugado
- Precisão vs tolerância
- Estabilidade numérica
- Sistemas bem-condicionados e mal-condicionados

## Critérios de Sucesso

Cada teste valida:

1. **Correção**: Resultado está correto (EXPECT_EQ, EXPECT_DOUBLE_EQ)
2. **Precisão**: Valores estão dentro da precisão (EXPECT_NEAR)
3. **Limites**: Valores estão dentro de intervalos esperados
4. **Comportamento**: Operações têm comportamento esperado

## Validação da Solução

Os testes de Gradiente Conjugado validam:

```cpp
// Para cada sistema
Vector x = conjugate_gradient::solve(A, b, tolerance, maxIter);
Vector Ax = multiplicar(A, x);

// Verificar: Ax ≈ b
for (int i = 0; i < n; i++) {
    EXPECT_NEAR(b.getValue(i), Ax.getValue(i), precision);
}
```

## Integração com CI/CD

Para adicionar ao CI (GitHub Actions, GitLab CI, etc):

```yaml
test:
  script:
    - cd build
    - cmake ..
    - cmake --build .
    - ./programa_testes --verbose
```

## Adicionando Novos Testes

```cpp
// Adicionar em tests/test_novo.cpp

#include <gtest/gtest.h>
#include "seu/header.hpp"

class MeuTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Setup antes de cada teste
    }
    void TearDown() override {
        // Cleanup após cada teste
    }
};

TEST_F(MeuTest, CasoEspecifico) {
    // Seu teste aqui
    EXPECT_EQ(resultado, esperado);
}
```

Então atualizar CMakeLists.txt para incluir o novo arquivo.

## Troubleshooting

### Teste falha aleatoriamente
- Verifique comparações de ponto flutuante com `EXPECT_NEAR`
- Aumente a tolerância se necessário

### Google Test não compila
```bash
# Limpar cache e reconfigurat
rm -rf build
mkdir build && cd build
cmake -DBUILD_TESTS=ON ..
cmake --build .
```

### Teste muito lento
- Verifique multiplicação matricial para matrizes grandes
- Use `-O3` para otimizações de compilação

## Documentação Google Test

- **Assertions**: EXPECT_EQ, EXPECT_NE, EXPECT_LT, EXPECT_LE, EXPECT_GT, EXPECT_GE
- **Floating Point**: EXPECT_DOUBLE_EQ, EXPECT_NEAR
- **Boolean**: EXPECT_TRUE, EXPECT_FALSE
- **String**: EXPECT_STREQ, EXPECT_STRNE
- **Fixtures**: TEST_F para reutilizar setup/teardown

## Cobertura de Código

Para gerar cobertura de código:

```bash
# Com GCC
cmake -DCMAKE_CXX_FLAGS="--coverage" ..
cmake --build .
./programa_testes
gcov build/CMakeFiles/programa_testes.dir/tests/*.cpp.o
```

## Performance

Tempo de execução típico:
- 12 testes de Matriz: < 1ms
- 14 testes de Vector: < 1ms
- 6 testes de Gradiente Conjugado: < 1ms
- **Total: < 1ms**

---

**Última atualização**: 2 de maio de 2026
