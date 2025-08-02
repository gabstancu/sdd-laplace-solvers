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
        log.solver_name    = name;
    }

    template<typename System>
    void solve(System& system)
    {   
        auto& A        = system.A;
        auto& b        = system.b;
        auto& u        = system.u;
        log.system_dim = system.A.rows();
        
        for (int k = 0; k < max_iters; k++)
        {
            auto start = std::chrono::high_resolution_clock::now();

            // log.num_of_iterations++;
            // log.res_per_iteration.push_back(r.norm() / b_norm);

            auto now = std::chrono::high_resolution_clock::now();
            double t = std::chrono::duration<double>(now - start).count();
            if (t > TIMEOUT) {
                log.timed_out = 1;
                break;
            }
        }
        // this->final_solution = u;
    }
};


#endif // MULTIGRID_HPP