#ifndef JACOBI_HPP
#define JACOBI_HPP

#include "utils/SolverLog.hpp"
template<typename Matrix, typename Vector>
struct Jacobi
{
    double tol = DEFAULT_TOL;
    int max_iters = MAX_ITERS;

    void solve(LinearSystem<Matrix, Vector>& system) const
    {   
        auto& A = system.A;
        auto& b = system.b;
        auto& u = system.u;
        Vector u_next = u;
        double sum, b_norm, res;

        b_norm = b.norm();
        res = (A* u - b).norm() / b_norm;

        if (res < tol) return;

        for (int k = 0; k < max_iters; k++)
        {
            for (int i = 0; i < A.rows(); i++)
            {   
                sum = 0;
                for (int j = 0; j < A.cols(); j++)
                {   
                    if (i==j) continue;
                    sum += A(i, j) * u(j);
                }
                u_next(i) = (1 / A(i, i)) * (b(i) - sum);
            }

            /* ----------- Matrix based ----------- */
            // u = D_inv * (b - (L + U)*u);

            res = (A * u_next - b).norm() / b_norm;
            if (res < tol) return;

            u = u_next; // now previous
        }
        return;
    }
};


#endif // JACOBI_HPP