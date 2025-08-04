#ifndef CONJUGATE_GRADIENT_HPP
#define CONJUGATE_GRADIENT_HPP

#include "utils/SolverLog.hpp"
#include "solvers/config.h"
template<typename Matrix, typename Vector>
struct ConjugateGradient
{
    double      tol       = DEFAULT_TOL;
    int         max_iters = 1e6;
    std::string name      = "CG";
    SolverLog<Eigen::VectorXd>   log;
    Vector      final_solution;

    ConjugateGradient ()
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

        std::cout << "max_iters: " << max_iters << '\n';

        Vector r      = b - A * u; // initial residual
        double b_norm = b.norm();
        double r_norm = r.norm();

        if (r_norm / b_norm <= tol) 
        {   
            this->final_solution = u;
            log.final_solution   = this->final_solution;
            log.converged        = 1;
            return;
        }

        Vector d = r; // initial search direction

        // auto start = std::chrono::high_resolution_clock::now();
        for (int k = 0; k < max_iters; k++)
        {   
            // std::cout << "--------------------- iter. " << k+1 << " ---------------------\n";
            Vector Ad     = A * d;
            Vector u_prev = u;

            double alpha      = ((r.transpose() * r) / (d.transpose() * Ad)).coeff(0); // step size
            double r_prev_dot = (r.transpose() * r).coeff(0); // to calculate beta

            u = u + alpha * d;
            r = r - alpha * Ad;

            r_norm = r.norm();
            
            log.num_of_iterations++;
            log.res_per_iteration.push_back(r_norm / b_norm);

            if (r_norm / b_norm <= tol) 
            {   
                log.converged        = 1;
                this->final_solution = u;
                log.final_solution   = this->final_solution;
                return;
            }

            double beta = r.dot(r) / r_prev_dot;
            d           = r + beta * d; // update direction
        }
        this->final_solution = u;
        log.final_solution   = this->final_solution;
        return;
    }
};


#endif // CONJUGATE_GRADIENT_HPP