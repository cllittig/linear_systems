#pragma once

// Macros de teste minimalistas para verificar solvers de sistemas lineares.
//
// Dois tipos de verificação são usados nos testes:
//   - Solução conhecida: compara x retornado pelo solver com x* esperado.
//   - Resíduo: calcula r = Ax - b e verifica que ||r|| < tol.
//     O resíduo é útil porque não depende de conhecer x* de antemão —
//     qualquer solver correto deve produzir ||Ax - b|| ≈ 0.

#include <cmath>
#include <cstdio>

namespace test {
    static int _passed = 0;
    static int _failed = 0;
}

#define RUN_TEST(fn) \
    do { printf("[TESTE] %s\n", #fn); (fn)(); } while (false)

#define ASSERT_TRUE(cond) \
    do { \
        if (cond) { \
            ++test::_passed; \
        } else { \
            printf("  FALHOU  %s:%d  =>  %s\n", __FILE__, __LINE__, #cond); \
            ++test::_failed; \
        } \
    } while (false)

#define ASSERT_NEAR(a, b, tol) \
    do { \
        double _a = (a), _b = (b), _tol = (tol); \
        if (std::fabs(_a - _b) <= _tol) { \
            ++test::_passed; \
        } else { \
            printf("  FALHOU  %s:%d  =>  %.8g != %.8g  (dif=%.2e, tol=%.2e)\n", \
                   __FILE__, __LINE__, _a, _b, std::fabs(_a - _b), _tol); \
            ++test::_failed; \
        } \
    } while (false)

// Compara dois Vector elemento a elemento
#define ASSERT_VEC_NEAR(v, expected, tol) \
    do { \
        bool _ok = true; \
        if ((v).getLength() != (expected).getLength()) { \
            printf("  FALHOU  %s:%d  =>  tamanhos diferentes (%d != %d)\n", \
                   __FILE__, __LINE__, (v).getLength(), (expected).getLength()); \
            ++test::_failed; \
            _ok = false; \
        } else { \
            for (int _i = 0; _i < (v).getLength() && _ok; _i++) { \
                double _vi = (v).getValue(_i), _ei = (expected).getValue(_i); \
                if (std::fabs(_vi - _ei) > (tol)) { \
                    printf("  FALHOU  %s:%d  =>  v[%d]=%.8g != %.8g  (dif=%.2e, tol=%.2e)\n", \
                           __FILE__, __LINE__, _i, _vi, _ei, std::fabs(_vi - _ei), (tol)); \
                    ++test::_failed; \
                    _ok = false; \
                } \
            } \
        } \
        if (_ok) ++test::_passed; \
    } while (false)

#define TEST_SUMMARY() \
    do { \
        printf("\n=== %d passou, %d falhou ===\n", test::_passed, test::_failed); \
        return test::_failed > 0 ? 1 : 0; \
    } while (false)
