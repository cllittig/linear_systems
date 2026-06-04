# linear_systems

Biblioteca C++ com métodos diretos e iterativos para sistemas lineares densos, desenvolvida como TCC.

**Dependências:** `g++` (C++17), `make`, `libopenblas-dev`

## Comandos

```bash
make              # compila libls.a
make test         # compila e executa todos os testes
make clean        # remove build/ e libls.a
```

```bash
make comparar     # roda benchmark (blas e lapack) e gera plot/comparar.csv
make figuras      # gera as figuras do artigo em artigo2/imagens/
```

## Métodos

| Namespace             | Método                  |
|-----------------------|-------------------------|
| `lu`                  | LU com pivoteamento parcial |
| `cholesky`            | Cholesky ($A = LL^T$)   |
| `gaussjordan`         | Gauss-Jordan            |
| `gaussseidel`         | Gauss-Seidel / SOR      |
| `conjugate_gradient`  | Gradientes Conjugados   |

## Variantes de build

```bash
make BLAS=0              # sem OpenBLAS
make BLAS=1 LAPACK=1     # com LAPACKE (necessário para make comparar)
make BUILD=debug         # sem otimização, com -g
```
