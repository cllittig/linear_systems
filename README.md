# Linear Systems Solver - Projecto de TCC

Sistema em C++ para resolução de sistemas lineares usando métodos numéricos avançados, com suporte a álgebra linear e testes unitários abrangentes.

## 📋 Características

- **Álgebra Linear**: Implementação de classes Matriz e Vector para operações matriciais
- **Métodos Numéricos**: Gradiente Conjugado, Cholesky, LU, Gauss-Jordan
- **Testes Unitários**: Suite completa com Google Test (32 testes)
- **Build Moderno**: CMake com suporte automático para Google Test
- **Bibliotecas Externas**: BLAS/LAPACK para operações otimizadas

## 🛠️ Pré-requisitos

- GCC 11+ ou Clang 13+
- CMake 3.14+
- Make ou Ninja
- OpenBLAS/LAPACK

## 📦 Instalação de Dependências

### Linux (Ubuntu/Debian)
```bash
# Dependências de compilação
sudo apt-get install -y build-essential cmake git

# BLAS/LAPACK
sudo apt-get install -y libopenblas-dev liblapack-dev
```

### macOS
```bash
brew install cmake openblas lapack gcc
```

## 🔨 Compilação

### Estrutura do Projeto

```
.
├── CMakeLists.txt              # Configuração CMake
├── include/                    # Headers
│   ├── algebra_linear/         # Classes Matriz e Vector
│   ├── metodos/                # Métodos numéricos
│   └── utils/                  # Utilitários
├── src/                        # Implementação
│   ├── algebra_linear/
│   ├── metodos/
│   └── main.cpp
└── tests/                      # Testes unitários
    ├── test_matriz.cpp
    ├── test_vector.cpp
    └── test_conjugate_gradient.cpp
```

### Build do Projeto

```bash
# Criar diretório de build
mkdir -p build && cd build

# Configurar com CMake
cmake ..

# Compilar projeto completo
cmake --build .

# Ou usar make diretamente
make -j$(nproc)
```

## 🧪 Executar Testes

```bash
# Depois de compilar, executar testes
./programa_testes

# Com mais verbosidade
./programa_testes --verbose

# Executar teste específico
./programa_testes --gtest_filter=MatrizTest.Addition
```

## 📊 Cobertura de Testes

### Matriz (12 testes)
- ✅ Construtores (default, com dimensões, cópia)
- ✅ Operações (set/get, adição, subtração)
- ✅ Multiplicação (escalar, matricial)
- ✅ Transformações (transposta)
- ✅ Propriedades (quadrada, simétrica)
- ✅ Identidade

### Vector (14 testes)
- ✅ Construtores e operações básicas
- ✅ Aritmética (adição, subtração, multiplicação)
- ✅ Produtos (escalar, elemento a elemento)
- ✅ Normas e operações AXPY
- ✅ Multiplicação matriz-vetor

### Gradiente Conjugado (6 testes)
- ✅ Sistemas 2x2 e 3x3
- ✅ Matriz identidade
- ✅ Sistemas escalados
- ✅ Estabilidade numérica
- ✅ Sistemas maiores (5x5)

**Total: 32 testes, 100% de aprovação**

## 🚀 Executar Programa Principal

```bash
./programa_tcc
```

O programa demonstra o uso do método Gradiente Conjugado para resolver um sistema linear 3x3 simétrico definido positivo.

## 📝 Exemplo de Uso

```cpp
#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"
#include "metodos/conjugate_gradient.hpp"

int main() {
    // Criar sistema Ax = b
    Matriz A(3, 3);
    A.setValue(0, 0, 4.0);   A.setValue(0, 1, 12.0);  A.setValue(0, 2, -16.0);
    A.setValue(1, 0, 12.0);  A.setValue(1, 1, 37.0);  A.setValue(1, 2, -43.0);
    A.setValue(2, 0, -16.0); A.setValue(2, 1, -43.0); A.setValue(2, 2, 98.0);

    Vector b(3);
    b.setValue(0, 0.0);
    b.setValue(1, 6.0);
    b.setValue(2, 39.0);

    // Resolver usando Gradiente Conjugado
    Vector x = conjugate_gradient::solve(A, b, 1e-10, 100);

    // Verificar solução
    Vector Ax = multiplicar(A, x);
    return 0;
}
```

## 🏗️ Arquitetura

### Classes Principais

**Matriz**
```cpp
class Matriz {
    - Dimensões: rows, columns
    - Armazenamento: std::vector<double> data
    - Operações: +, -, *, transpose(), determinant(), inverse()
    - Propriedades: isSquare(), isSimetric()
};
```

**Vector**
```cpp
class Vector {
    - Dimensão: length
    - Armazenamento: std::vector<double> data
    - Operações: +, -, *, (produo escalar, AXPY)
    - Propriedades: euclidian_length(), linear_product()
};
```

**Métodos Numéricos**
```cpp
namespace conjugate_gradient {
    Vector solve(const Matriz &A, const Vector &b, 
                 double tolerance, int maxIterations, const Vector* x0);
}
```

## 🔧 Configuração CMake

O CMakeLists.txt automaticamente:
- Encontra OpenBLAS/LAPACK
- Baixa e compila Google Test (versão 1.12.1)
- Configura flags de otimização (-O3)
- Gera testes executáveis

## 📚 Métodos Numéricos Implementados

| Método | Arquivo | Status |
|--------|---------|--------|
| Gradiente Conjugado | `src/metodos/conjugate_gradient.cpp` | ✅ Testado |
| Decomposição LU | `src/metodos/lu.cpp` | ✅ Implementado |
| Decomposição Cholesky | `src/metodos/cholesky.cpp` | ✅ Implementado |
| Gauss-Jordan | `src/metodos/gauss_jordan.cpp` | ✅ Implementado |

## 🐛 Troubleshooting

### Google Test não encontrado
```bash
# O CMake faz download automático, mas se falhar:
cd build
cmake -DBUILD_TESTS=ON ..
```

### Erros de linkagem OpenBLAS
```bash
# Verificar instalação
pkg-config --cflags --libs openblas
# Se necessário, instalar desenvolvimento
sudo apt-get install libopenblas-dev
```

## 📈 Performance

- Testes executam em < 1ms
- Suporte para múltiplas tolerâncias
- Iterações adaptáveis para Gradiente Conjugado
- Integração com BLAS para operações matriciais otimizadas

## 🎓 Referências

- **Google Test**: https://github.com/google/googletest
- **OpenBLAS**: https://github.com/xianyi/OpenBLAS
- **Métodos Numéricos**: Algorithm Development and Data Transformation (LAPACK/BLAS)

## 📄 Licença

Este projeto segue a licença do arquivo LICENSE no repositório.

## 👨‍💻 Autor

Projeto desenvolvido como Trabalho de Conclusão de Curso.

## 📞 Suporte

Para problemas ou sugestões, adicione issues no repositório.

---

**Última atualização**: 2 de maio de 2026
