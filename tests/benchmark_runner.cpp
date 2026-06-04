// benchmark_runner.cpp
//
// Benchmark automatizado para o TCC. Gera sistemas Ax=b com A simétrica e
// estritamente diagonal dominante (logo, Positiva Definida), resolve o sistema
// pelos 5 métodos numéricos da biblioteca e exporta tempo de execução e norma
// do resíduo para arquivos .csv em data/.
//
// Métodos: LU, Cholesky, Gauss-Seidel, SOR e Gradientes Conjugados.

#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <random>
#include <string>
#include <vector>

#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"
#include "metodos/cholesky.hpp"
#include "metodos/conjugate_gradient.hpp"
#include "metodos/gauss_seidel.hpp"
#include "metodos/lu.hpp"

namespace {

// --- Parâmetros do benchmark ---------------------------------------------
// Granularidade aumentada para observar melhor a curva de complexidade
const std::vector<int> SIZES    = {250, 500, 750, 1000, 1500, 2000, 2500, 3000,3500, 4000,4500, 5000,5500,6000,6500,7000,7500};
const int              NUM_RUNS = 10;    // Quantidade de execuções por tamanho de matriz
const double           TOL      = 1e-8;  // tolerância dos métodos iterativos
const int              MAX_ITER = 10000; // teto de iterações
const double           SOR_OMEGA= 1.5;   // fator de relaxação para o SOR
const unsigned         SEED     = 42;    // semente fixa: execução reprodutível

// Gera A (n x n) simétrica e ESTRITAMENTE diagonal dominante.
Matriz gerarMatrizSPD(int n, std::mt19937& rng) {
    std::uniform_real_distribution<double> off(-1.0, 1.0);
    std::uniform_real_distribution<double> diagBoost(1.0, 2.0);

    Matriz A(n, n);
    std::vector<double> somaAbsLinha(n, 0.0);

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double v = off(rng);
            A.setValue(i, j, v);
            A.setValue(j, i, v); // simetria
            somaAbsLinha[i] += std::abs(v);
            somaAbsLinha[j] += std::abs(v);
        }
    }
    for (int i = 0; i < n; ++i) {
        A.setValue(i, i, somaAbsLinha[i] + diagBoost(rng)); // dominância estrita
    }
    return A;
}

// Vetor exato x preenchido com 1.0.
Vector gerarXExato(int n) {
    Vector x(n);
    for (int i = 0; i < n; ++i) x.setValue(i, 1.0);
    return x;
}

// Norma euclidiana do resíduo ||b - A x||.
double normaResiduo(const Matriz& A, const Vector& b, const Vector& x) {
    Vector r = b - multiplicar(A, x);
    return r.euclidian_length();
}

struct Resultado {
    double tempo_s;
    double residuo;
};

// Cronometra a resolução do sistema
template <typename Solve>
Resultado medir(const Matriz& A, const Vector& b, Solve&& solve) {
    using clock = std::chrono::steady_clock;
    auto t0 = clock::now();
    Vector x = solve(A, b);
    auto t1 = clock::now();
    double tempo_s = std::chrono::duration<double>(t1 - t0).count();
    return {tempo_s, normaResiduo(A, b, x)};
}

} // namespace

int main() {
    std::ofstream tempos("data/benchmark_tempos.csv");
    std::ofstream residuos("data/benchmark_residuos.csv");
    if (!tempos || !residuos) {
        std::fprintf(stderr,
                     "Erro: não foi possível abrir os arquivos em data/.\n"
                     "Execute o binário a partir da raiz do projeto.\n");
        return 1;
    }

    // Adicionada a coluna "Run" para rastrear as múltiplas execuções
    tempos   << "Tamanho_N,Metodo,Run,Tempo_s\n";
    residuos << "Tamanho_N,Metodo,Run,Norma_Residuo\n";

    struct Metodo {
        std::string nome;
        std::function<Vector(const Matriz&, const Vector&)> resolver;
    };
    
    const std::vector<Metodo> metodos = {
        {"LU",       [](const Matriz& A, const Vector& b) { return lu::solve(A, b); }},
        {"Cholesky", [](const Matriz& A, const Vector& b) { return cholesky::solve(A, b); }},
        {"Gauss-Seidel", [](const Matriz& A, const Vector& b) { return gaussseidel::solve(A, b, TOL, MAX_ITER, nullptr, 1.0); }},
        {"SOR", [](const Matriz& A, const Vector& b) { return gaussseidel::solve(A, b, TOL, MAX_ITER, nullptr, SOR_OMEGA); }},
        {"Gradientes Conjugados", [](const Matriz& A, const Vector& b) { return conjugate_gradient::solve(A, b, TOL, MAX_ITER, nullptr, nullptr, false); }},
    };

    std::mt19937 rng(SEED);

    for (int n : SIZES) {
        std::printf("=== Avaliando tamanho N = %d (%d execucoes) ===\n", n, NUM_RUNS);
        std::fflush(stdout);

        for (int run = 1; run <= NUM_RUNS; ++run) {
            // A cada "run", uma nova matriz é gerada para evitar vício nos dados
            Matriz A = gerarMatrizSPD(n, rng);
            Vector xExato = gerarXExato(n);
            Vector b = multiplicar(A, xExato); 

            std::printf("  [Run %02d/%02d] -> ", run, NUM_RUNS);

            for (const Metodo& m : metodos) {
                Resultado r = medir(A, b, m.resolver);

                tempos   << n << ',' << m.nome << ',' << run << ',' << r.tempo_s  << '\n';
                residuos << n << ',' << m.nome << ',' << run << ',' << r.residuo  << '\n';
                
                std::printf("%s ", m.nome.c_str());
            }
            
            tempos.flush();
            residuos.flush();
            std::printf("[OK]\n");
            std::fflush(stdout);
        }
        std::printf("\n");
    }

    std::printf("Concluído. Dados salvos.\n");
    return 0;
}