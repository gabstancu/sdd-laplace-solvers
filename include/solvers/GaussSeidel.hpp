#ifndef GAUSS_SEIDEL_HPP
#define GAUSS_SEIDEL_HPP

#include "utils/SolverLog.hpp"
#include "solvers/config.h"
template<typename Vector>
struct GaussSeidel
{   
    using Scalar       = typename Vector::Scalar;
    using SparseMatrix = Eigen::SparseMatrix<Scalar, Eigen::RowMajor>;

    double            tol          = DEFAULT_TOL;
    int               max_iters    = MAX_ITERS;
    std::string       name         = "Gauss-Seidel";
    SolverLog<Vector> log;
    Vector            final_solution;

    template<typename System>
    GaussSeidel (System system)
    {
        log.tolerance      = tol;
        log.system_dim     = system.A.rows();
        max_iters          = static_cast<int>(std::min(int(std::pow(log.system_dim, 1)), 50000));
        log.max_iterations = max_iters;
        log.solver_name    = name;
    }

    template<typename System>
    void solve(System& system)
    {   
        const auto& A        = system.A;
        const auto& b        = system.b;
              auto& u        = system.u;

        log.system_dim = system.A.rows();
        
        std::cout << "max_iters: " << max_iters << '\n';

        double sum1,   sum2;
        double b_norm, r_norm, res;

        b_norm = b.norm();
        res    = (A * u - b).norm() / b_norm;

        if (res <= tol) 
        {   
            this->final_solution = u;
            log.final_solution   = this->final_solution;
            log.converged = 1;
            return;
        }

        Vector inv_diag = A.diagonal().cwiseInverse();

        for (int k = 0; k < max_iters; k++)
        {   
            for (int i = 0; i < A.rows(); ++i)
            {
                double sum = 0;
                
                for (typename SparseMatrix::InnerIterator it(A, i); it; ++it)
                {   
                    int j = it.col();
                    if (j != i)
                    {
                        sum += it.value() * u[it.col()];
                    }
                }
                // /* for clarity... */
                // double sum1 = 0, sum2 = 0;
                
                // for (typename SparseMatrix::InnerIterator it(A, i); it; ++it)
                // {   
                //     int j = it.col();
                //     if (j < i)
                //     {
                //         sum1 += it.value() * u[j];
                //     }
                //     else if (j > i)
                //     {
                //         sum2 += it.value() * u[j];
                //     }
                // }

                u[i] = inv_diag[i] * (b[i] - sum);
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