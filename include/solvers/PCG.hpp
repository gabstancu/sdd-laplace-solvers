#ifndef PCG_HPP
#define PCG_HPP

#include "basis/Preconditioners.hpp"
#include "utils/SolverLog.hpp"
#include "solvers/config.h"
template<typename Matrix, typename Vector, typename Precondition>
struct PCG
{
    double        tol        = DEFAULT_TOL;
    int           max_iters  = 1e6;
    std::string   name       = "PCG";
    Precondition& precon;
    std::string   precon_name;
    SolverLog<Eigen::VectorXd>   log;
    Vector        final_solution;
    
    PCG (Precondition& p, std::string precon_name) : precon(p) 
    {
        log.tolerance      = tol;
        log.max_iterations = max_iters;
        this->precon_name  = precon_name;
        log.solver_name    = name;
        log.precon_name    = precon_name; 
    };

    template<typename System>
    void solve (System& system)
    {   
        auto& A        = system.A;
        auto& b        = system.b;
        auto& u        = system.u;
        log.system_dim = system.A.rows();

        std::cout << "max_iters: " << max_iters << '\n';

        Vector r       = b - A * u;
        Vector zeta    = precon.apply(r);
        Vector d       = zeta;

        double b_norm = b.norm();
        if (r.norm() / b_norm < tol) 
        {   
            log.converged = 1;
            this->final_solution = u;
            log.final_solution   = this->final_solution;
            return;
        }

        // auto start = std::chrono::high_resolution_clock::now();
        for (int k = 0; k < max_iters; k++)
        {   
            // std::cout << "Iteration: " << k+1 << '\n';
            Vector Ad        = A * d;
            Vector u_prev    = u;
            Vector r_prev    = r;
            Vector zeta_prev = zeta;

            double alpha = r.dot(zeta) / d.dot(Ad);
            u = u + alpha * d;
            r = r - alpha * Ad;

            double diff = (u - u_prev).norm() / u.norm();

            log.num_of_iterations++;
            log.res_per_iteration.push_back(r.norm() / b_norm);
            log.diff_per_iteration.push_back(diff);

            if (r.norm() / b_norm <= tol || diff <= tol) 
            {   
                log.converged = 1;
                this->final_solution = u;
                log.final_solution   = this->final_solution;
                return;
            }

            // auto now = std::chrono::high_resolution_clock::now();
            // double t = std::chrono::duration<double>(now - start).count();
            // if (t >= TIMEOUT) {
            //     log.timed_out = 1;
            //     break;
            // }

            zeta        = precon.apply(r);
            double beta = r.dot(zeta) / r_prev.dot(zeta_prev);
            d           = zeta + beta * d;
        }
        this->final_solution = u;
        log.final_solution   = this->final_solution;
        return;
    }
};


#endif // PCG_HPP