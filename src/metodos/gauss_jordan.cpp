
#include "metodos/gauss_jordan.hpp"
#include "algebra_linear/matriz.hpp"

namespace gaussjordan {

void swap(int a, int b) {
}

Matriz gauss_jordan(Matriz matriz, int current_row) {

    // h := 1 /* Initialization of the pivot row */
    // k := 1 /* Initialization of the pivot column */
    //
    // while h ≤ m and k ≤ n
    //     /* Find the k-th pivot: */
    //     i_max := argmax (i = h ... m, abs(A[i, k]))
    //     if A[i_max, k] = 0
    //         /* No pivot in this column, pass to next column */
    //         k := k + 1
    //     else
    //         swap rows(h, i_max)
    //         /* Do for all rows below pivot: */
    //         for i = h + 1 ... m:
    //             f := A[i, k] / A[h, k]
    //             /* Fill with zeros the lower part of pivot column: */
    //             A[i, k] := 0
    //             /* Do for all remaining elements in current row: */
    //             for j = k + 1 ... n:
    //                 A[i, j] := A[i, j] - A[h, j] * f
    //         /* Increase pivot row and column */
    //         h := h + 1
    //         k := k + 1
    // auto h = 1, k = 1;
    // auto m = matriz.getRows(), n = matriz.getColumns();
    //
    // while (h <= m && k <= n) {
    //     auto max = -1;
    //     auto i_max = -1;
    //     for (auto i = 0; i < m; i++) {
    //         if (matriz.getValue(i, k) > max) {
    //             i_max = i;
    //             max = std::abs(matriz.getValue(i, k));
    //        }
    //     }
    //     if (matriz.getValue(i_max, k) == 0) {
    //         k++;
    //     } else {
    //         matriz.data.swap(vector<double> &x);
    //     }
    // }
    
    return Matriz(0, 0);
}

};
