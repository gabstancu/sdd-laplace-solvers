#ifndef JACOBI_HPP
#define JACOBI_HPP

#include "utils/SolverLog.hpp"
#include "solvers/config.h"
template<typename Matrix, typename Vector>
struct Jacobi
{
    double      tol             = DEFAULT_TOL;
    int         max_iters       = MAX_ITERS;
    std::string name            = "Jacobi";
    Vector      final_solution;

    SolverLog<Eigen::VectorXd>   log;

    template<typename System>
    Jacobi (System system)
    {
        log.tolerance      = tol;
        log.system_dim     = system.A.rows();
        max_iters          = int(std::min(int(2 * std::pow(log.system_dim, 1)), 50000));
        log.max_iterations = max_iters;
        log.solver_name    = name;
    }

    template<typename System>
    void solve(System& system)
    {   
        auto& A        = system.A;
        auto& b        = system.b;
        auto& u        = system.u;
        Vector u_next(u.size());
        double sum, b_norm, res;

        std::cout << "max_iters: " << max_iters << '\n';

        Vector inv_diag = A.diagonal().cwiseInverse();
        b_norm = b.norm();
        res    = (A * u - b).norm() / b_norm;

        if (res <= tol) 
        {   
            log.converged = 1;
            this->final_solution = u;
            log.final_solution   = this->final_solution;
            return;
        }

        for (int k = 0; k < max_iters; k++)
        {   
            for (int i = 0; i < A.rows(); i++)
            {   
                sum = 0.0;
                for (int j = 0; j < A.cols(); j++)
                {   
                    if (i!=j) 
                        sum += A(i, j) * u[j];
                }
                u_next[i] = inv_diag[i] * (b[i] - sum);
            }

            res = (A * u_next - b).norm() / b_norm;

            log.num_of_iterations++;
            log.res_per_iteration.push_back(res);
            
            if (res <= tol) 
            {   
                log.converged = 1;
                this->final_solution = u_next;
                log.final_solution   = this->final_solution;
                return;
            }

            u.swap(u_next); // now previous
        }
        this->final_solution = u;
        log.final_solution   = this->final_solution;
        return;
    }
};


#endif // JACOBI_HPP