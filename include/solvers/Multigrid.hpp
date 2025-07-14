#ifndef MULTIGRID_HPP
#define MULTIGRID_HPP

#include "utils/SolverLog.hpp"
#include "solvers/config.h"
template<typename Matrix, typename Vector>
struct Multigrid
{
    double      tol       = DEFAULT_TOL;
    int         max_iters = MAX_ITERS;
    std::string name      = "Multigrid";
    SolverLog<Eigen::VectorXd>   log;
    Vector      final_solution;

    Multigrid ()
    {
        log.tolerance      = tol;
        log.max_iterations = max_iters;
    }

    template<typename System>
    void solve(System& system)
    {
        
        for (int k = 0; k < max_iters; k++)
        {
            auto start = std::chrono::high_resolution_clock::now();

            // log.num_of_iterations++;
            // log.res_per_iteration.push_back(r.norm() / b_norm);

            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
        }
        // this->final_solution = u;
    }
};


#endif // MULTIGRID_HPP