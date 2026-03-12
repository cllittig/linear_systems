#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <typeinfo>

namespace cholesky{

	Matriz cholesky(const Matriz &A) {
    if(!A.isSquare()) throw std::runtime_error("A matriz A deve ser quadrada.");
    if(!A.isSimetric()) throw std::runtime_error("A matriz A deve ser simétrica.");

    const int n = A.getRows();
    Matriz L(n, n); 

    for (int i = 0; i < n; i++) { 
        for (int j = 0; j <= i; j++) {
            double sum = 0;

            if (i == j) {
                for (int k = 0; k < j; k++) {
                    double val = L.getValue(j, k);
                    sum += val * val;
                }
                double res = A.getValue(j, j) - sum;
                if (res < 0) throw std::runtime_error("Matriz não é definida positiva.");
                L.setValue(j, j, std::sqrt(res));
            } else {
                for (int k = 0; k < j; k++) {
                    sum += L.getValue(i, k) * L.getValue(j, k);
                }
                L.setValue(i, j, (A.getValue(i, j) - sum) / L.getValue(j, j));
            }
        }
    }
    	return L;
	}

	Vector forwardSubstituition(const Matriz &L, const Vector &b){
		int n = L.getRows();

		Vector y(n);

		for(int i=0; i <n; i++){
			double sum =0;
			for(int j=0; j<i; j++){
				sum+=L.getValue(i,j) * y.getValue(j);
			}
			y.setValue(i, (b.getValue(i) - sum)/ L.getValue(i,i));
		}
		return y;
	}

	Vector backwardSubstitution(const Matriz &Lt, const Vector &y){
		int n = Lt.getRows();
		Vector x(n);

		for(int i = n-1; i>= 0; i--){
			double sum= 0;
			for(int j = i+1; j< n; j++){
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
}

