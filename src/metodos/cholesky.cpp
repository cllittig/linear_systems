#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <typeinfo>
#include <limits>

namespace cholesky{

	Matriz cholesky(const Matriz &A) {
    if(!A.isSquare()) throw std::runtime_error("A matriz A deve ser quadrada.");
    if(!A.isSimetric()) throw std::runtime_error("A matriz A deve ser simétrica.");

    const int n = A.getRows();
    Matriz L(n, n); 

	// Compute a relative threshold based on A magnitude to avoid strict absolute checks
	double global_max = 0.0;
	for (int ii = 0; ii < n; ++ii)
	  for (int jj = 0; jj < n; ++jj)
		global_max = std::max(global_max, std::abs(A.getValue(ii, jj)));
	double rel_eps = std::numeric_limits<double>::epsilon();
	double tol = rel_eps * (global_max == 0.0 ? 1.0 : global_max) * 10.0;

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j <= i; ++j) {
			double sum = 0.0;

			if (i == j) {
				for (int k = 0; k < j; ++k) {
					double val = L.getValue(j, k);
					sum += val * val;
				}
				double res = A.getValue(j, j) - sum;
				if (res < -tol) throw std::runtime_error("Matriz não é definida positiva.");
				if (res < tol) res = 0.0; // treat very small negatives/zeros as zero
				double diag = std::sqrt(res);
				if (!(diag > tol)) throw std::runtime_error("Cholesky falhou: diagonal não positiva/ muito pequena.");
				L.setValue(j, j, diag);
			} else {
				for (int k = 0; k < j; ++k) {
					sum += L.getValue(i, k) * L.getValue(j, k);
				}
				double denom = L.getValue(j, j);
				if (std::abs(denom) < tol) throw std::runtime_error("Cholesky falhou: denominador (L[j][j]) muito pequeno.");
				L.setValue(i, j, (A.getValue(i, j) - sum) / denom);
			}
		}
	}
    	return L;
	}

	Vector forwardSubstituition(const Matriz &L, const Vector &b){
		int n = L.getRows();

		Vector y(n);

		for (int i = 0; i < n; ++i) {
			double sum = 0;
			for (int j = 0; j < i; ++j) {
				sum+=L.getValue(i,j) * y.getValue(j);
			}
			y.setValue(i, (b.getValue(i) - sum)/ L.getValue(i,i));
		}
		return y;
	}

	Vector backwardSubstitution(const Matriz &Lt, const Vector &y){
		int n = Lt.getRows();
		Vector x(n);

		for (int i = n - 1; i >= 0; --i) {
			double sum = 0;
			for (int j = i + 1; j < n; ++j) {
				sum += Lt.getValue(i,j) * x.getValue(j);
			}
			x.setValue(i,(y.getValue(i) - sum)/Lt.getValue(i,i));
		}
		return x;
	}

	Vector solve(const Matriz &A, const Vector &b){
		Matriz L = cholesky(A);
		Vector y = forwardSubstituition(L, b);

		Matriz Lt = L.transpose();
		Vector x = backwardSubstitution(Lt, y);

		return x;
	}

	Matriz solve(const Matriz &A, const Matriz &B) {
		// Solve multiple RHS using Cholesky factorization
		if (B.getRows() != A.getRows()) throw std::invalid_argument("Dimensoes incompativeis entre A e B.");
		Matriz L = cholesky(A);
		Matriz Lt = L.transpose();
		int n = A.getRows();
		int m = B.getColumns();
		Matriz X(n, m);
		for (int col = 0; col < m; ++col) {
			Vector bcol(n);
			for (int i = 0; i < n; ++i) bcol.setValue(i, B.getValue(i, col));
			Vector y = forwardSubstituition(L, bcol);
			Vector x = backwardSubstitution(Lt, y);
			for (int i = 0; i < n; ++i) X.setValue(i, col, x.getValue(i));
		}
		return X;
	}
}

