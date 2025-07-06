#ifndef PCG_HPP
#define PCG_HPP

#include "basis/Preconditioners.hpp"
#include "utils/SolverLog.hpp"
template<typename Matrix, typename Vector, typename Precondition>
struct PCG
{
    double tol = DEFAULT_TOL;
    int max_iters = MAX_ITERS;
    std::string name = "PCG";
    Precondition precon;
    SolverLog log;
    

    PCG (const Precondition& p) : precon(p) 
    {
        log.tolerance = tol;
        log.max_iterations = max_iters;
    };

    void solve (LinearSystem<Matrix, Vector>& system)
    {   
        auto start = std::chrono::high_resolution_clock::now();
        auto& A = system.A;
        auto& b = system.b;
        auto& u = system.u;

        Vector r = b - A * u;
        Vector zeta = precon.apply(r);
        Vector d = zeta;

        double b_norm = b.norm();
        if (r.norm() / b_norm < tol) 
        {   
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            log.time_elapsed = elapsed;
            log.converged = 1;
            return;
        }

        for (int k = 0; k < max_iters; k++)
        {   
            // std::cout << "Iteration: " << k+1 << '\n';
            Vector Ad = A * d;
            Vector r_prev = r;
            Vector zeta_prev = zeta;

            double alpha = r.dot(zeta) / d.dot(Ad);
            u = u + alpha * d;
            r = r - alpha * Ad;

            log.num_of_iterations++;
            log.res_per_iteration.push_back(r.norm() / b_norm);

            if (r.norm() / b_norm < tol) 
            {   
                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = end - start;
                log.time_elapsed = elapsed;
                log.converged = 1;
                return;
            }

            zeta = precon.apply(r);
            double beta = r.dot(zeta) / r_prev.dot(zeta_prev);
            d = zeta + beta * d;
        }
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        log.time_elapsed = elapsed;
        return;
    }
};


#endif // PCG_HPP