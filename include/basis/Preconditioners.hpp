#ifndef PRECONDITIONERS_HPP
#define PRECONDITIONERS_HPP

#include <Eigen/IterativeLinearSolvers>


using SparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;

template<typename Vector>
struct IdentityPreconditioner
{
    Vector apply(Vector r)
    { 
        return r; 
    }
};

template<typename Vector>
struct DiagonalPreconditioner
{
    Vector M; // precon. to be applied

    DiagonalPreconditioner(const SparseMatrix& A) 
    { 
        M = A.diagonal().cwiseInverse(); 
    }

    Vector apply (Vector r) 
    { 
        return M.cwiseProduct(r); 
    }
};

template<typename Vector>
struct SSORPreconditioner 
{
    SparseMatrix L, LT;
    Eigen::DiagonalMatrix<typename Vector::Scalar, Eigen::Dynamic> Dinv;
    double omega;

    SSORPreconditioner(const SparseMatrix& A, double omega_) : omega(omega_) 
    {
        // Get diagonal
        Eigen::VectorXd Dvec = A.diagonal();
        Dinv = Dvec.cwiseInverse().asDiagonal();

        // Create L = (D/ω) + strict lower triangular
        SparseMatrix Ltmp(A.rows(), A.cols());
        Ltmp.reserve(A.nonZeros());
        
        for(int k=0; k<A.outerSize(); ++k) 
        {
            for(SparseMatrix::InnerIterator it(A,k); it; ++it) 
            {
                if(it.row() == it.col()) 
                {
                    Ltmp.insert(it.row(), it.col()) = it.value()/omega;
                }
                else if(it.row() > it.col()) 
                {
                    Ltmp.insert(it.row(), it.col()) = it.value();
                }
            }
        }
        
        L = Ltmp;
        LT = L.transpose();
    }

    Vector apply(Vector r) 
    {
        Vector y = L.triangularView<Eigen::Lower>().solve(r);
        y = Dinv * y;
        Vector z = LT.triangularView<Eigen::Upper>().solve(y);
        return (omega/(2.0-omega)) * z;
    }
};


// template<typename Vector>
// struct IncompleteCholeskyPreconditioner 
// {
//     Eigen::IncompleteCholesky<double> ichol;

//     IncompleteCholeskyPreconditioner(Matrix& A) 
//     {
//         ichol.compute(A);
//     }

//     // Explicitly delete copy constructor/assignment
//     IncompleteCholeskyPreconditioner(const IncompleteCholeskyPreconditioner&) = delete;
//     IncompleteCholeskyPreconditioner& operator=(const IncompleteCholeskyPreconditioner&) = delete;

//     // Default move constructor/assignment
//     IncompleteCholeskyPreconditioner(IncompleteCholeskyPreconditioner&&) = default;
//     IncompleteCholeskyPreconditioner& operator=(IncompleteCholeskyPreconditioner&&) = default;

//     Vector apply(Vector r) 
//     {
//         return ichol.solve(r);
//     }
// };
template<typename Vector>
struct IncompleteCholeskyPreconditioner 
{
    Eigen::IncompleteCholesky<double> ichol;

    IncompleteCholeskyPreconditioner(const SparseMatrix& A) 
    {
        // Convert to column-major if needed
        if(A.IsRowMajor) {
            SparseMatrix colA = A;
            ichol.compute(colA);
        }
        else {
            ichol.compute(A);
        }
    }

    Vector apply(Vector r) 
    {
        return ichol.solve(r);
    }
};


#endif // PRECONDITIONERS_HPP