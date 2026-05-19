#ifndef SOLVER_UTILS_HPP
#define SOLVER_UTILS_HPP

#include <cstddef>

namespace solver {
struct IterationInfo {
    int iterations = 0;
    double final_residual_norm = 0.0;
    bool converged = false;
};

} // namespace solver

#endif // SOLVER_UTILS_HPP
