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
        
        for(int k = 0; k < A.outerSize(); ++k) 
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
        
        L  = Ltmp;
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


template<typename Vector>
struct IncompleteCholeskyPreconditioner 
{   

    using Scalar       = typename Vector::Scalar;
    using Triplet      = Eigen::Triplet<Scalar>;

    SparseMatrix L; // lower triangular factor

    IncompleteCholeskyPreconditioner(const SparseMatrix&  A) 
    {
        int n = A.rows();
        L.resize(n, n);

        std::vector<Triplet> triplets;
        triplets.reserve(A.nonZeros() / 2 + n);

        for (int i = 0; i < A.rows(); ++i)
        {
            for (typename SparseMatrix::InnerIterator it(A, i); it; ++it)
            {
                int j = it.col();

                if (i >= j)
                {
                    triplets.emplace_back(i, j, it.value());
                }
            }
        }
        L.setFromTriplets(triplets.begin(), triplets.end());

        // IC(0) factorisation, no fill-in
        for (int j = 0; j < n; ++j) {
            // Diagonal entry
            Scalar& L_jj = L.coeffRef(j, j);
            L_jj = std::sqrt(L_jj);

            // Update column j
            for (int i = j + 1; i < n; ++i) {
                if (L.coeff(i, j) != 0) 
                {
                    L.coeffRef(i, j) /= L_jj;
                }
            }

            // Update remaining submatrix
            for (int k = j + 1; k < n; ++k) {
                if (L.coeff(k, j) != 0) 
                {
                    const Scalar L_kj = L.coeff(k, j);
                    for (typename SparseMatrix::InnerIterator it(L, k); it && it.col() <= j; ++it) 
                    {
                        const int i = it.col();
                        if (i == j) 
                        {
                            it.valueRef() -= L_kj * L_kj;
                        }
                    }
                }
            }
        }
    }


    Vector apply(const Vector& r) const 
    {
        Vector y = L.triangularView<Eigen::Lower>().solve(r);
        return L.transpose().template triangularView<Eigen::Upper>().solve(y);
    }


    // Explicitly delete copy constructor/assignment
    IncompleteCholeskyPreconditioner(const IncompleteCholeskyPreconditioner&) = delete;
    IncompleteCholeskyPreconditioner& operator=(const IncompleteCholeskyPreconditioner&) = delete;

    // Default move constructor/assignment
    IncompleteCholeskyPreconditioner(IncompleteCholeskyPreconditioner&&) = default;
    IncompleteCholeskyPreconditioner& operator=(IncompleteCholeskyPreconditioner&&) = default;
};


#endif // PRECONDITIONERS_HPP