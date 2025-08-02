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
        max_iters          = int(10 * log.system_dim);
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

        double sum1, sum2;
        double b_norm = b.norm();
        double res    = (A * u - b).norm() / b_norm;

        if (res <= tol) 
        {   
            this->final_solution = u;
            log.final_solution   = this->final_solution;
            log.converged = 1;
            return;
        }

        // Matrix L = A.triangularView<Eigen::StrictlyLower>();
        // Matrix U = A.triangularView<Eigen::StrictlyUpper>();
        // Matrix D = A.diagonal().asDiagonal();
        // Matrix K = (D + omega * L).inverse();
        
        // auto start = std::chrono::high_resolution_clock::now();
        for (int k = 0; k < max_iters; k++)
        {   
            Vector u_prev = u;
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
                u(i) = (1 - omega) * u(i) + omega * ((1 / A(i, i)) * (b(i) - sum1 - sum2));
            }

             /* ----------- matrix based ----------- */
             // u = K * (omega * b - (omega * U + (omega - 1) * D) * u);
             
            res = (A * u - b).norm() / b.norm();
            double diff = (u - u_prev).norm() / u_prev.norm();

            log.num_of_iterations++;
            log.res_per_iteration.push_back(res);
            log.diff_per_iteration.push_back(diff);


            if (res <= tol) 
            {   
                this->final_solution = u;
                log.final_solution   = this->final_solution;
                log.converged = 1;
                return;
            }

            // auto now = std::chrono::high_resolution_clock::now();
            // double t = std::chrono::duration<double>(now - start).count();
            // if (t > TIMEOUT) {
            //     log.timed_out = 1;
            //     break;
            // }
        }
        this->final_solution = u;
        log.final_solution   = this->final_solution;
        return;
    }
};


#endif // SOR_HPP