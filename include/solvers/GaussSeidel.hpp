#ifndef GAUSS_SEIDEL_HPP
#define GAUSS_SEIDEL_HPP

template<typename Matrix, typename Vector>
struct GaussSeidel
{
    double tol = DEFAULT_TOL;
    int max_iters = MAX_ITERS;

    void solve(LinearSystem<Matrix, Vector>& system) const
    {
        auto& A = system.A;
        auto& b = system.b;
        auto& u = system.u;

        double sum1, sum2;
        double b_norm, r_norm, res;

        b_norm = b.norm();
        res = (A * u - b).norm() / b_norm;

        if (res < tol) return;

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
                u(i) = (1 / A(i, i)) * (b(i) - sum1 - sum2);
            }

            /* ----------- matrix based ----------- */
            // u = L_inv * (b - U * u);

            res = (A * u - b).norm() / b_norm;
            if (res < tol) return;
        }
        return;
    }
};


#endif // GAUSS_SEIDEL_HPP