#ifndef GAUSS_SEIDEL_HPP
#define GAUSS_SEIDEL_HPP

#include "utils/SolverLog.hpp"
#include "solvers/config.h"
template<typename Matrix, typename Vector>
struct GaussSeidel
{
    double      tol          = DEFAULT_TOL;
    int         max_iters    = MAX_ITERS;
    std::string name         = "Gauss-Seidel";
    SolverLog<Eigen::VectorXd> log;
    Vector      final_solution;

    template<typename System>
    GaussSeidel (System system)
    {
        log.tolerance      = tol;
        log.system_dim     = system.A.rows();
        max_iters          = int(std::min(int(std::pow(log.system_dim, 1)), 50000));
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
        
        std::cout << "max_iters: " << max_iters << '\n';

        double sum1,   sum2;
        double b_norm, r_norm, res;

        Vector inv_diag = A.diagonal().cwiseInverse();

        b_norm = b.norm();
        res    = (A * u - b).norm() / b_norm;

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
                sum1 = 0; sum2 = 0;
                for (int j = 0; j < i; j++)
                {
                    sum1 += A(i, j) * u[j];
                }

                for (int j = i+1; j < A.rows(); j++)
                {   
                    sum2 += A(i, j) * u[j];
                }
                u[i] = inv_diag[i] * (b[i] - sum1 - sum2);
            }

            res = (A * u - b).norm() / b_norm;
            log.num_of_iterations++;
            log.res_per_iteration.push_back(res);

            if (res <= tol) 
            {   
                log.converged = 1;
                this->final_solution = u;
                log.final_solution   = this->final_solution;
                return;
            }
        }
        this->final_solution = u;
        log.final_solution   = this->final_solution;
        return;
    }
};


#endif // GAUSS_SEIDEL_HPP