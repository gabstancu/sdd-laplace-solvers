#ifndef JACOBI_HPP
#define JACOBI_HPP

#include "utils/SolverLog.hpp"
#include "solvers/config.h"
template<typename Matrix, typename Vector>
struct Jacobi
{
    double      tol       = DEFAULT_TOL;
    int         max_iters = MAX_ITERS;
    std::string name      = "Jacobi";
    Vector      final_solution;
    SolverLog<Eigen::VectorXd>   log;

    Jacobi ()
    {
        log.tolerance      = tol;
        log.max_iterations = max_iters;
    }

    template<typename System>
    void solve(System& system)
    {   
        auto& A = system.A;
        auto& b = system.b;
        auto& u = system.u;
        Vector u_next = u;
        double sum, b_norm, res;

        b_norm = b.norm();
        res    = (A* u - b).norm() / b_norm;

        if (res < tol) 
        {   
            log.converged = 1;
            this->final_solution = u;
            log.final_solution   = this->final_solution;
            return;
        }

        for (int k = 0; k < max_iters; k++)
        {   
            // std::cout << "--------------------- iter " << k+1 << "---------------------\n";
            for (int i = 0; i < A.rows(); i++)
            {   
                sum = 0;
                for (int j = 0; j < A.cols(); j++)
                {   
                    if (i==j) 
                        continue;
                    sum += A(i, j) * u(j);
                }
                u_next(i) = (1 / A(i, i)) * (b(i) - sum);
            }

            /* ----------- Matrix based ----------- */
            // u = D_inv * (b - (L + U)*u);

            res = (A * u_next - b).norm() / b_norm;
            // std::cout << res << '\n';

            log.num_of_iterations++;
            log.res_per_iteration.push_back(res);
            
            if (res < tol) 
            {   
                log.converged = 1;
                this->final_solution = u;
                log.final_solution   = this->final_solution;
                return;
            }

            u = u_next; // now previous
        }
        this->final_solution = u;
        log.final_solution   = this->final_solution;
        return;
    }
};


#endif // JACOBI_HPP