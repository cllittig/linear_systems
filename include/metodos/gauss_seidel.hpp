#ifndef GAUSSSEIDEL_HPP
#define GAUSSSEIDEL_HPP

#include <linear_systems/algebra_linear/matriz.hpp>
#include <linear_systems/algebra_linear/vector.hpp>
#include <linear_systems/metodos/solver_utils.hpp>


namespace gaussseidel{
    Matriz gauss_seidel(const Matriz& A, const Matriz& b,
                        double tol = 1e-6, int maxIter = 1000,
                        solver::IterationInfo* info = nullptr, double omega = 1.0);
    Vector solve(const Matriz& A, const Vector& b,
                 double tol = 1e-6, int maxIter = 1000,
                 solver::IterationInfo* info = nullptr, double omega = 1.0);
}

#endif //GAUSSSEIDEL_HPP
