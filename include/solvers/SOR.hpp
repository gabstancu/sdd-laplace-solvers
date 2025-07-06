#ifndef SOR_HPP
#define SOR_HPP

#include "utils/SolverLog.hpp"
template<typename Matrix, typename Vector>
struct SOR
{
    double tol = DEFAULT_TOL;
    int max_iters = MAX_ITERS;
    double omega;
    std::string name = "SOR";
    SolverLog log;
    
    SOR (double omega) : omega(omega)
    {
        log.max_iterations = max_iters;
        log.tolerance = tol;
    }

    void solve(LinearSystem<Matrix, Vector>& system)
    {   
        auto& A = system.A;
        auto& b = system.b;
        auto& u = system.u;

        double sum1, sum2;
        double b_norm = b.norm();
        double res = (A * u - b).norm() / b_norm;

        if (res < tol) 
        {   
            log.converged = 1;
            return;
        }

        // Matrix L = A.triangularView<Eigen::StrictlyLower>();
        // Matrix U = A.triangularView<Eigen::StrictlyUpper>();
        // Matrix D = A.diagonal().asDiagonal();
        // Matrix K = (D + omega * L).inverse();

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
                u(i) = (1 - omega) * u(i) + omega * ((1 / A(i, i)) * (b(i) - sum1 - sum2));
            }

             /* ----------- matrix based ----------- */
             // u = K * (omega * b - (omega * U + (omega - 1) * D) * u);
             
            res = (A * u - b).norm() / b.norm();
            log.num_of_iterations++;
            log.res_per_iteration.push_back(res);
            if (res < tol) 
            {   
                log.converged = 1;
                return;
            }
        }
        return;
    }
};


#endif // SOR_HPP