#ifndef SOR_HPP
#define SOR_HPP

#include "utils/SolverLog.hpp"
#include "solvers/config.h"
template<typename Matrix, typename Vector>
struct SOR
{
    double      tol       = DEFAULT_TOL;
    int         max_iters = MAX_ITERS;
    double      omega;
    std::string name      = "SOR";
    SolverLog<Eigen::VectorXd>   log;
    Vector      final_solution;
    
    template<typename System>
    SOR (System system) : omega(system.omega_)
    {   
        log.system_dim     = system.A.rows();
        max_iters          = int(10 * std::sqrt(log.system_dim));
        log.max_iterations = max_iters;
        log.tolerance      = tol;
        log.solver_name    = name;
    }

    template<typename System>
    void solve(System& system)
    {   
        auto& A        = system.A;
        auto& b        = system.b;
        auto& u        = system.u;

        std::cout << "max_iters: " << max_iters << '\n';

        double sum1, sum2;
        Vector inv_diag = A.diagonal().cwiseInverse();
        double b_norm   = b.norm();
        double res      = (A * u - b).norm() / b_norm;

        if (res <= tol) 
        {   
            this->final_solution = u;
            log.final_solution   = this->final_solution;
            log.converged = 1;
            return;
        }

        for (int k = 0; k < max_iters; k++)
        {   
            for (int i = 0; i < A.rows(); i++)
            {
                sum1 = 0.0; sum2 = 0.0;
                for (int j = 0; j < i; j++)
                {
                    sum1 += A(i, j) * u[j];
                }

                for (int j = i+1; j < A.rows(); j++)
                {
                    sum2 += A(i, j) * u[j];
                }
                u[i] = (1 - omega) * u[i] + omega * (inv_diag[i] * (b[i] - sum1 - sum2));
            }
             
            res = (A * u - b).norm() / b.norm();

            log.num_of_iterations++;
            log.res_per_iteration.push_back(res);


            if (res <= tol) 
            {   
                this->final_solution = u;
                log.final_solution   = this->final_solution;
                log.converged = 1;
                return;
            }
        }
        this->final_solution = u;
        log.final_solution   = this->final_solution;
        return;
    }
};


#endif // SOR_HPP