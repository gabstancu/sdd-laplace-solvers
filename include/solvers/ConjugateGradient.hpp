#ifndef CONJUGATE_GRADIENT_HPP
#define CONJUGATE_GRADIENT_HPP

#include "utils/SolverLog.hpp"
template<typename Matrix, typename Vector>
struct ConjugateGradient
{
    double tol = DEFAULT_TOL;
    int max_iters = MAX_ITERS;
    std::string name = "CG";
    SolverLog log;

    ConjugateGradient ()
    {
        log.tolerance = tol;
        log.max_iterations = max_iters;
    }

    void solve(LinearSystem<Matrix, Vector>& system)
    {   
        auto start = std::chrono::high_resolution_clock::now();
        auto& A = system.A;
        auto& b = system.b;
        auto& u = system.u;

        Vector r = b - A * u; // initial residual
        double b_norm = b.norm();
        double r_norm = r.norm();

        if (r_norm / b_norm < tol) 
        {   
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            log.time_elapsed = elapsed;
            log.converged = 1;
            return;
        }

        Vector d = r; // initial search direction

        for (int k = 0; k < max_iters; k++)
        {
            // std::cout << "--------------------- iter. " << k+1 << " ---------------------\n";
            Vector Ad = A * d;

            double alpha = ((r.transpose() * r) / (d.transpose() * Ad)).coeff(0); // step size
            double r_prev_dot = (r.transpose() * r).coeff(0); // to calculate beta

            u = u + alpha * d;
            r = r - alpha * Ad;

            r_norm = r.norm();
            
            log.num_of_iterations++;
            log.res_per_iteration.push_back(r.norm() / b_norm);

            if (r_norm / b_norm < tol) 
            {   
                log.converged = 1;
                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = end - start;
                log.time_elapsed = elapsed;
                return;
            }

            double beta = r.dot(r) / r_prev_dot;

            d = r + beta * d; // update direction
        }
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        log.time_elapsed = elapsed;
        return;
    }
};


#endif // CONJUGATE_GRADIENT_HPP