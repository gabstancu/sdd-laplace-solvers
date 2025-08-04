#ifndef PCG_HPP
#define PCG_HPP

#include "basis/Preconditioners.hpp"
#include "utils/SolverLog.hpp"
#include "solvers/config.h"
template<typename Vector, typename Precondition>
struct PCG
{
    using Scalar       = typename Vector::Scalar;
    using SparseMatrix = Eigen::SparseMatrix<Scalar, Eigen::RowMajor>;

    double              tol        = DEFAULT_TOL;
    int                 max_iters  = 1e6;
    std::string         name       = "PCG";
    Precondition&       precon;
    std::string         precon_name;
    SolverLog<Vector>   log;
    Vector              final_solution;
    
    PCG (Precondition& p, std::string precon_name) : precon(p) , precon_name(precon_name)
    {
        log.tolerance      = tol;
        log.max_iterations = max_iters;
        log.solver_name    = name;
        log.precon_name    = precon_name; 
    };

    template<typename System>
    void solve (System& system)
    {   
        const auto& A        = system.A;
        const auto& b        = system.b;
              auto& u        = system.u;
        log.system_dim       = system.A.rows();

        std::cout << "max_iters: " << max_iters << '\n';

        Vector r(b.size()), zeta(b.size()), d(b.size());
        Vector Ad(A.rows());

        r.noalias()       = b - A * u;
        zeta.noalias()    = precon.apply(r);
        d.noalias()       = zeta;

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
            Ad.noalias()        = A * d;
            Vector r_prev       = r;
            Vector zeta_prev    = zeta;

            double alpha = r.dot(zeta) / d.dot(Ad);
            u.noalias() += alpha * d;
            r.noalias() -= alpha * Ad;

            log.num_of_iterations++;
            log.res_per_iteration.push_back(r.norm() / b_norm);

            if (r.norm() / b_norm <= tol) 
            {   
                log.converged = 1;
                this->final_solution = u;
                log.final_solution   = this->final_solution;
                return;
            }

            zeta.noalias()        = precon.apply(r);
            double beta           = r.dot(zeta) / r_prev.dot(zeta_prev);
            d.noalias()           = zeta + beta * d;
        }
        this->final_solution = u;
        log.final_solution   = this->final_solution;
        return;
    }
};


#endif // PCG_HPP