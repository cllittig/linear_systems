# Documentação da biblioteca `linear_systems`

## Visão geral

Esta biblioteca em C++ resolve sistemas lineares da forma

$$
A x = b
$$

Ela reúne estruturas de álgebra linear, métodos numéricos e testes automatizados com validação no GNU Octave.

## O que ela faz

As classes `Matriz` e `Vector` dão suporte às operações básicas do projeto, como criação, acesso a valores, soma, subtração, multiplicação e transposição.

Os métodos implementados são:

- `cholesky`
- `lu`
- `gaussjordan`
- `gaussseidel`
- `conjugate_gradient`

## Como funciona

- **Cholesky**: usado quando `A` é quadrada, simétrica e definida positiva. Faz a decomposição `A = L L^T` e resolve dois sistemas triangulares.
- **LU**: fatora `A` em matrizes triangulares com pivotamento parcial.
- **Gauss-Jordan**: monta a matriz aumentada `[A | b]` e reduz o sistema até a forma final.
- **Gauss-Seidel**: método iterativo com tolerância, número máximo de iterações e fator de relaxação `omega`.
- **Gradiente Conjugado**: método iterativo para matrizes quadradas, simétricas e definidas positivas. Suporta chute inicial, diagnóstico, Jacobi e múltiplos lados direitos.

## Validação

O projeto foi validado com testes unitários em Google Test e com comparação externa em GNU Octave no arquivo [tests/test_compare_with_octave.cpp](tests/test_compare_with_octave.cpp).

## Como usar

### Exemplo com `conjugate_gradient`

```cpp
#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"
#include "metodos/conjugate_gradient.hpp"

int main() {
    Matriz A(3, 3);
    A.setValue(0, 0, 4.0);
    A.setValue(0, 1, 1.0);
    A.setValue(0, 2, 0.0);
    A.setValue(1, 0, 1.0);
    A.setValue(1, 1, 3.0);
    A.setValue(1, 2, 1.0);
    A.setValue(2, 0, 0.0);
    A.setValue(2, 1, 1.0);
    A.setValue(2, 2, 2.0);

    Vector b(3);
    b.setValue(0, 1.0);
    b.setValue(1, 2.0);
    b.setValue(2, 3.0);

    Vector x = conjugate_gradient::solve(A, b);
}
```

### Exemplo com diagnóstico

```cpp
solver::IterationInfo info;
Vector x = conjugate_gradient::solve(A, b, 1e-10, 1000, nullptr, &info, true);
```

## Quando usar cada método

- **Cholesky**: quando `A` for simétrica definida positiva.
- **LU**: para sistemas quadrados em geral.
- **Gauss-Jordan**: para problemas pequenos ou uso didático.
- **Gauss-Seidel**: para sistemas iterativos, se houver convergência.
- **Gradiente Conjugado**: para matrizes simétricas definidas positivas, principalmente em sistemas maiores.

## Como compilar e testar

```bash
cd build
cmake --build .
ctest --output-on-failure -j1
```

## Arquivos principais

- [include/algebra_linear/matriz.hpp](include/algebra_linear/matriz.hpp)
- [include/algebra_linear/vector.hpp](include/algebra_linear/vector.hpp)
- [include/metodos/cholesky.hpp](include/metodos/cholesky.hpp)
- [include/metodos/lu.hpp](include/metodos/lu.hpp)
- [include/metodos/gauss_jordan.hpp](include/metodos/gauss_jordan.hpp)
- [include/metodos/gauss_seidel.hpp](include/metodos/gauss_seidel.hpp)
- [include/metodos/conjugate_gradient.hpp](include/metodos/conjugate_gradient.hpp)
- [tests/test_compare_with_octave.cpp](tests/test_compare_with_octave.cpp)

## Resumo

O projeto resolve sistemas lineares com diferentes métodos numéricos, testes em C++ e validação com Octave.
