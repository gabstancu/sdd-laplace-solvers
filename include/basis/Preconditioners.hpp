#ifndef PRECONDITIONERS_HPP
#define PRECONDITIONERS_HPP

template<typename Matrix, typename Vector>
struct IdentityPreconditioner
{
    Vector apply(const Vector& r) const { 
        return r; 
    }
};

template<typename Matrix, typename Vector>
struct DiagonalPreconditioner
{
    Matrix M; // precon. to be applied

    DiagonalPreconditioner(const Matrix& A) { 
        M = A.diagonal().asDiagonal().inverse(); 
    }

    Vector apply (const Vector& r) const { 
        return M * r; 
    }
};



#endif // PRECONDITIONERS_HPP