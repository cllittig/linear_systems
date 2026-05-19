#ifndef GAUSSSEIDEL_HPP
#define GAUSSSEIDEL_HPP

#include "algebra_linear/matriz.hpp"
#include "algebra_linear/vector.hpp"
<<<<<<< HEAD
#include "metodos/solver_utils.hpp"


namespace gaussseidel{
    Matriz gauss_seidel(const Matriz& A, const Matriz& b,
                        double tol = 1e-6, int maxIter = 1000,
                        solver::IterationInfo* info = nullptr, double omega = 1.0);
    Vector solve(const Matriz& A, const Vector& b,
                 double tol = 1e-6, int maxIter = 1000,
                 solver::IterationInfo* info = nullptr, double omega = 1.0);
=======

namespace gaussseidel {
    Matriz gauss_seidel(const Matriz& A, const Matriz& b);
    Vector solve(const Matriz& A, const Vector& b);
>>>>>>> a0306ec7ad214946c9df52e0635fa4c496afe0d2
}

#endif //GAUSSSEIDEL_HPP
