#ifndef PRECONDITIONERS_HPP
#define PRECONDITIONERS_HPP

#include <Eigen/IterativeLinearSolvers>

template<typename Matrix, typename Vector>
struct IdentityPreconditioner
{
    Vector apply(Vector r)
    { 
        return r; 
    }
};

template<typename Matrix, typename Vector>
struct DiagonalPreconditioner
{
    Matrix M; // precon. to be applied

    DiagonalPreconditioner(Matrix& A) 
    { 
        M = A.diagonal().asDiagonal().inverse(); 
    }

    Vector apply (Vector r) 
    { 
        return M * r; 
    }
};

template<typename Matrix, typename Vector>
struct SSORPreconditioner 
{
    Matrix M;
    double omega;

    SSORPreconditioner(Matrix& A, double omega_) : omega(omega_) 
    {
        Eigen::MatrixXd D = A.diagonal().asDiagonal();
        Matrix L = -(A.template triangularView<Eigen::StrictlyLower>().toDenseMatrix());
        M = (D - omega * L).inverse() * D * (D - omega * L).transpose().inverse();
    }

    Vector apply(Vector r) 
    {
        return M * r;
    }
};

template<typename Matrix, typename Vector>
struct IncompleteCholeskyPreconditioner 
{
    Eigen::IncompleteCholesky<double> ichol;

    IncompleteCholeskyPreconditioner(Matrix& A) 
    {
        ichol.compute(A);
    }

    // Explicitly delete copy constructor/assignment
    IncompleteCholeskyPreconditioner(const IncompleteCholeskyPreconditioner&) = delete;
    IncompleteCholeskyPreconditioner& operator=(const IncompleteCholeskyPreconditioner&) = delete;

    // Default move constructor/assignment
    IncompleteCholeskyPreconditioner(IncompleteCholeskyPreconditioner&&) = default;
    IncompleteCholeskyPreconditioner& operator=(IncompleteCholeskyPreconditioner&&) = default;

    Vector apply(Vector r) 
    {
        return ichol.solve(r);
    }
};



#endif // PRECONDITIONERS_HPP