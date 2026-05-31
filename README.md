# Linear Systems Solver - Projeto de TCC

Biblioteca em C++ para resolucao de sistemas lineares e operacoes de algebra linear, com implementacao de metodos numericos classicos, testes unitarios proprios e benchmark.

## Visao Geral

- Algebra linear com classes `Matriz` e `Vector`
- Solvers numericos em `src/metodos/`
- Build baseado em `Makefile`
- Suite de testes sem GoogleTest (harness proprio em `include/utils/test.hpp`)
- Integracao opcional com OpenBLAS (ativada por padrao no Makefile)

## Estrutura do Projeto

```text
.
|- Makefile
|- include/
|  |- algebra_linear/
|  |  |- matriz.hpp
|  |  |- sistema_linear.hpp
|  |  `- vector.hpp
|  |- metodos/
|  |  |- cholesky.hpp
|  |  |- conjugate_gradient.hpp
|  |  |- gauss_jordan.hpp
|  |  |- gauss_seidel.hpp
|  |  |- lu.hpp
|  |  `- solver_utils.hpp
|  `- utils/
|     |- read_data.hpp
|     `- test.hpp
|- src/
|  |- algebra_linear/
|  |- metodos/
|  `- utils/
|- tests/
|  |- test_cholesky.cpp
|  |- test_conjugate_gradient.cpp
|  |- test_gauss_jordan.cpp
|  |- test_gauss_seidel.cpp
|  |- test_lu.cpp
|  |- test_matriz.cpp
|  |- test_matriz_advanced.cpp
|  |- test_vector.cpp
|  |- test_vector_advanced.cpp
|  `- benchmark.cpp
`- plot/
   |- plot.py
   `- requirements.txt
```

## Pre-requisitos

- Compilador C++ com suporte a C++17 (`g++` ou `clang++`)
- `make`
- OpenBLAS/LAPACK (quando `BLAS=1`)

### Linux (Ubuntu/Debian)

```bash
sudo apt-get update
sudo apt-get install -y build-essential libopenblas-dev liblapack-dev
```

### macOS (Homebrew)

```bash
brew install gcc openblas lapack
```

## Compilacao

### Build padrao (release)

```bash
make
```

Gera a biblioteca estatica `libls.a` com os objetos compilados em `build/`.

### Build de debug

```bash
make BUILD=debug
```

### Build sem BLAS

```bash
make BLAS=0
```

## Testes

Para compilar e executar todos os testes:

```bash
make test
```

Observacoes:

- Cada arquivo `tests/test_*.cpp` gera um executavel em `build/tests/`.
- Os testes usam macros de assercao definidas em `include/utils/test.hpp`.
- O alvo `test` falha se qualquer binario de teste retornar erro.

## Benchmark

Executa benchmark principal:

```bash
make benchmark
```

Tambem e possivel gerar grafico via alvo de plot (usa `plot/.venv/bin/python`):

```bash
make plot/data.png
```

## Instalacao

Instala biblioteca e headers em `/usr/local`:

```bash
sudo make install
```

Remove a instalacao:

```bash
sudo make uninstall
```

## Metodos Implementados

- `conjugate_gradient` (`src/metodos/conjugate_gradient.cpp`)
- `cholesky` (`src/metodos/cholesky.cpp`)
- `lu` (`src/metodos/lu.cpp`)
- `gaussjordan` (`src/metodos/gauss_jordan.cpp`)
- `gaussseidel` (`src/metodos/gauss_seidel.cpp`)

## Exemplo Rapido

```cpp
#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"
#include "metodos/conjugate_gradient.hpp"

int main() {
    Matriz A(2, 2);
    A.setValue(0, 0, 4.0);
    A.setValue(0, 1, 1.0);
    A.setValue(1, 0, 1.0);
    A.setValue(1, 1, 3.0);

    Vector b(2);
    b.setValue(0, 1.0);
    b.setValue(1, 2.0);

    Vector x = conjugate_gradient::solve(A, b, 1e-10, 1000);
    (void)x;
    return 0;
}
```

## Troubleshooting

### Erro de link com OpenBLAS

- Verifique se `libopenblas-dev` (Linux) ou `openblas` (macOS) esta instalado.
- Teste compilacao sem BLAS: `make BLAS=0`.

### Ambiente sujo de build

```bash
make clean
make test
```

## Licenca

Consulte o arquivo `LICENSE` do repositorio.

## Autor

Projeto desenvolvido como Trabalho de Conclusao de Curso.
