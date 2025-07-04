#ifndef CONJUGATE_GRADIENT_HPP
#define CONJUGATE_GRADIENT_HPP

template<typename Matrix, typename Vector>
struct ConjugateGradient
{
    double tol = DEFAULT_TOL;
    int max_iters = MAX_ITERS;

    void solve(LinearSystem<Matrix, Vector>& system) const
    {
        auto& A = system.A;
        auto& b = system.b;
        auto& u = system.u;

        Vector r = b - A * u; // initial residual
        double b_norm = b.norm();
        double r_norm = r.norm();

        if (r_norm / b_norm < tol) return;

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
            if (r_norm / b_norm < tol) return;

            double beta = r.dot(r) / r_prev_dot;

            d = r + beta * d; // update direction
        }
        return;
    }
};


#endif // CONJUGATE_GRADIENT_HPP