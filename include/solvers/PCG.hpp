#ifndef PCG_HPP
#define PCG_HPP

#include "basis/Preconditioners.hpp"

template<typename Matrix, typename Vector, typename Precondition>
struct PCG
{
    double tol = DEFAULT_TOL;
    int max_iters = MAX_ITERS;
    Precondition precon;

    PCG (const Precondition& p) : precon(p) {};

    void solve(LinearSystem<Matrix, Vector>& system) const
    {
        auto& A = system.A;
        auto& b = system.b;
        auto& u = system.u;

        Vector r = b - A * u;
        Vector zeta = precon.apply(r);
        Vector d = zeta;

        double b_norm = b.norm();
        if (r.norm() / b_norm < tol) return;

        for (int k = 0; k < max_iters; k++)
        {   
            // std::cout << "Iteration: " << k+1 << '\n';
            Vector Ad = A * d;
            Vector r_prev = r;
            Vector zeta_prev = zeta;

            double alpha = r.dot(zeta) / d.dot(Ad);
            u = u + alpha * d;
            r = r - alpha * Ad;

            if (r.norm() / b_norm < tol) return;

            zeta = precon.apply(r);
            double beta = r.dot(zeta) / r_prev.dot(zeta_prev);
            d = zeta + beta * d;
        }
        return;
    }
};


#endif // PCG_HPP