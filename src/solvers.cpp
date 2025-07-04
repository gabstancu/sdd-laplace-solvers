#include "solvers.hpp"

// TODO: write Gauss-Seidel (omega = 1) and SOR (omega != 1) to be the same 
// TODO: restructure code to templated version

Eigen::MatrixXd Conjugate_Gradient (Eigen::MatrixXd A, Eigen::VectorXd b, Eigen::VectorXd u_0)
{
    Eigen::VectorXd u = u_0;

    // compute initial residual
    Eigen::VectorXd r = b - A * u;

    double r_norm = r.norm();
    if (r_norm / b.norm() < TOL)
    {
        return u;
    }

    // set initial search direction
    Eigen::VectorXd d = r;

    for (int k = 0; k < MAXITERS; k++)
    {
        // std::cout << "--------------------- iter. " << k+1 << " ---------------------\n";
        Eigen::VectorXd Ad = A * d;

        double alpha = ((r.transpose() * r) / (d.transpose() * Ad)).coeff(0); // step size
        double r_prev_dot = (r.transpose() * r).coeff(0); // to calculate beta

        u = u + alpha * d;
        r = r - alpha * Ad;

        r_norm = r.norm();
        
        if (r_norm / b.norm() < TOL)
        {
            // std::cout << u << '\n';
            return u;
        }

        double beta = r.dot(r) / r_prev_dot;

        d = r + beta * d; // update direction
    }

    return u;
}


Eigen::MatrixXd Preconditioned_Conjugate_Gradient (Eigen::MatrixXd A, Eigen::VectorXd b, Eigen::VectorXd u_0, Eigen::MatrixXd precon)
{
    Eigen::VectorXd u = u_0;

    Eigen::VectorXd r = b - A * u;
    Eigen::VectorXd zeta = precon * r;
    Eigen::VectorXd d = zeta;

    double r_norm = r.norm();
    if (r_norm / b.norm() < TOL)
    {
        return u;
    }

    for (int k = 0; k < MAXITERS; k++)
    {   
        // std::cout << "--------------------- iter. " << k+1 << " ---------------------\n";

        Eigen::VectorXd Ad = A * d;
        Eigen::VectorXd r_prev = r;
        Eigen::VectorXd zeta_prev = zeta;

        double alpha = r.dot(zeta) / d.dot(Ad);
        u = u + alpha * d;
        r = r - alpha * Ad;

        r_norm = r.norm();
        if (r_norm / b.norm() < TOL)
        {
            return u;
        }
        
        zeta = precon * r;
        double beta = r.dot(zeta) / r_prev.dot(zeta_prev);
        d = zeta + beta * d;
    }
    
    return u;
}


Eigen::MatrixXd Jacobi (Eigen::MatrixXd A, Eigen::VectorXd b, Eigen::VectorXd u_0)
{
    Eigen::VectorXd u = u_0;
    Eigen::VectorXd u_next = u;
    double sum;

    double res = (A* u - b).norm() / b.norm();
    if (res < TOL)
    {
        return u;
    }

    // Eigen::MatrixXd D = A.diagonal().asDiagonal();
    // // Eigen::MatrixXd D_inv = D.inverse();
    // Eigen::MatrixXd L = A.triangularView<Eigen::StrictlyLower>();
    // Eigen::MatrixXd U = A.triangularView<Eigen::StrictlyUpper>();

    for (int k = 0; k < MAXITERS; k++)
    {
        // std::cout << "--------------------- iter. " << k+1 << " ---------------------\n";
        /* ----------- element based ----------- */
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

        res = (A * u_next - b).norm() / b.norm();
        if (res < TOL)
        {
            return u;
        }
        u = u_next; // now previous
    }

    return u;
}


Eigen::MatrixXd Gauss_Seidel (Eigen::MatrixXd A, Eigen::VectorXd b, Eigen::VectorXd u_0)
{
    Eigen::VectorXd u = u_0;
    double sum1, sum2;
    // Eigen::MatrixXd L = A.triangularView<Eigen::Lower>();
    // Eigen::MatrixXd L_inv = L.inverse();
    // Eigen::MatrixXd U = A.triangularView<Eigen::StrictlyUpper>();

    double res = (A * u - b).norm() / b.norm();
    if (res < TOL)
    {
        return u;
    }

    for (int k = 0; k < MAXITERS; k++)
    {   
        // std::cout << "--------------------- iter. " << k+1 << " ---------------------\n";
        /* ----------- element based ----------- */
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

        res = (A * u - b).norm() / b.norm();
        if (res < TOL)
        {
            return u;
        }
    }

    return u;
}


Eigen::MatrixXd SOR (Eigen::MatrixXd A, Eigen::VectorXd b, Eigen::VectorXd u_0, double omega)
{
    Eigen::VectorXd u = u_0;
    double sum1, sum2;

    // Eigen::MatrixXd L = A.triangularView<Eigen::StrictlyLower>();
    // Eigen::MatrixXd U = A.triangularView<Eigen::StrictlyUpper>();
    // Eigen::MatrixXd D = A.diagonal().asDiagonal();
    // Eigen::MatrixXd K = (D + omega * L).inverse();

    double res = (A * u - b).norm() / b.norm();
    if (res < TOL)
    {
        return u;
    }

    for (int k = 0; k < MAXITERS; k++)
    {
        std::cout << "--------------------- iter. " << k+1 << " ---------------------\n";
        /* ----------- element based ----------- */
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
        if (res < TOL)
        {
            return u;
        }
    }
    return u;
}


Eigen::MatrixXd Multigrid (Eigen::MatrixXd A, Eigen::VectorXd b, Eigen::VectorXd u_0)
{
    Eigen::VectorXd u = u_0;



    return u;
}