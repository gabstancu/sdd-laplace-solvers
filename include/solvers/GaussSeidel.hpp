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

    GaussSeidel ()
    {
        log.tolerance      = tol;
        log.max_iterations = max_iters;
        log.solver_name    = name;
    }

    template<typename System>
    void solve(System& system)
    {   
        auto start = std::chrono::high_resolution_clock::now();
        auto& A        = system.A;
        auto& b        = system.b;
        auto& u        = system.u;
        log.system_dim = system.A.rows();

        double sum1,   sum2;
        double b_norm, r_norm, res;

        b_norm = b.norm();
        res    = (A * u - b).norm() / b_norm;

        if (res < tol) 
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
                    sum1 += A(i, j) * u(j);
                }

                for (int j = i+1; j < A.rows(); j++)
                {   
                    sum2 += A(i, j) * u(j);
                }
                u(i) = (1 / A(i, i)) * (b(i) - sum1 - sum2);
            }

            /* ----------- matrix based ----------- */
            // u = L_inv * (b - U * u);

            res = (A * u - b).norm() / b_norm;
            log.num_of_iterations++;
            log.res_per_iteration.push_back(res);

            if (res < tol) 
            {   
                log.converged = 1;
                this->final_solution = u;
                log.final_solution   = this->final_solution;
                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = end - start;
                log.time_per_iteration.push_back(elapsed);
                return;
            }
            
            auto now = std::chrono::high_resolution_clock::now();
            double elapsed = std::chrono::duration<double>(now - start).count();
            if (elapsed > TIMEOUT) {
                log.timed_out = 1;
                break;
            }

            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            log.time_per_iteration.push_back(elapsed);
        }
        this->final_solution = u;
        log.final_solution   = this->final_solution;
        return;
    }
};


#endif // GAUSS_SEIDEL_HPP