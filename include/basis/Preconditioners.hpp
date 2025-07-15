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

    Matrix L, LT;
    Eigen::DiagonalMatrix<typename Vector::Scalar, Eigen::Dynamic> Dinv;
    double omega;

    // SSORPreconditioner(Matrix& A, double omega_) : omega(omega_) 
    // {
    //     Matrix Diag_A  = A.diagonal().asDiagonal();
    //     Matrix Upper_A = A.template triangularView<Eigen::StrictlyUpper>();
    //     Matrix Lower_A = A.template triangularView<Eigen::StrictlyLower>();
    //     this->M        = (omega_ / (2 - omega_)) * ((1/omega_) * Diag_A + Lower_A) * Diag_A.inverse() * ((1/omega_) * Diag_A + Lower_A).transpose();
    // }

    // Vector apply(Vector r) 
    // {
    //     return M.inverse() * r;
    // }

    SSORPreconditioner(Matrix& A, double omega_) : omega(omega_) 
    {
        Eigen::VectorXd Dvec = A.diagonal();
        Dinv = Dvec.cwiseInverse().asDiagonal();

        Matrix D = Dvec.asDiagonal();
        Matrix Lstrict = A.template triangularView<Eigen::StrictlyLower>();
        L = ((1.0 / omega) * D + Lstrict).sparseView();

        LT = L.transpose();
    }

    Vector apply(Vector r) 
    {
        // Solve: L y = r   (forward substitution)
        Vector y = L.template triangularView<Eigen::Lower>().solve(r);

        // Scale: y = D^{-1} y
        y = Dinv * y;

        // Solve: LT z = y  (backward substitution)
        Vector z = LT.template triangularView<Eigen::Upper>().solve(y);

        // Return SSOR preconditioned vector
        return (omega / (2.0 - omega)) * z;
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