#ifndef CHOLESKY_HPP
#define CHOLESKY_HPP

#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"

namespace cholesky{
	Matriz cholesky(const Matriz &A );
	Vector forwardSubstituition(const Matriz &L, const Vector &b);
	Vector backwardSubstitution(const Matriz &Lt, const Vector &y);	
	Vector solve(const Matriz &A, const Vector &b);

}

#endif //CHOLESKY_HPP
