#ifndef PRECONDITIONERS_HPP
#define PRECONDITIONERS_HPP

#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>

using SparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;

template<typename Vector>
struct IdentityPreconditioner
{
    Vector apply(Vector& r)
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

    Vector apply (Vector& r) 
    { 
        return M.cwiseProduct(r); 
    }
};


template<typename Vector, typename SparseMatrix = Eigen::SparseMatrix<typename Vector::Scalar>>
struct SSORPreconditioner
{
    using Scalar  = typename Vector::Scalar;
    using DiagMat = Eigen::DiagonalMatrix<Scalar, Eigen::Dynamic>;

    // We store (D + ω * strictLower(A)) as Lower triangular,
    // and reuse its transpose as Upper triangular.
    SparseMatrix DL;                 // Lower triangular: diag(D) + ω * A_lower
    DiagMat      D;                  // Diagonal matrix (not inverse)
    double       omega;

    SSORPreconditioner(const SparseMatrix& A, double omega_)
        : omega(omega_)
    {
        if (omega <= 0.0 || omega >= 2.0)
            throw std::invalid_argument("SSOR: omega must be in (0,2)");

        const int n = A.rows();
        if (A.rows() != A.cols())
            throw std::invalid_argument("SSOR: A must be square");

        // Build D (positive for SPD Laplacian)
        Eigen::Matrix<Scalar,Eigen::Dynamic,1> d = A.diagonal();
        if ((d.array() <= Scalar(0)).any())
            throw std::runtime_error("SSOR: nonpositive diagonal in A");
        D = d.asDiagonal();

        // Build DL = D + ω * strictLower(A)  (note: A_lower is negative for Laplacian)
        DL.resize(n,n);
        std::vector<Eigen::Triplet<Scalar>> trips;
        trips.reserve(A.nonZeros()); // safe upper bound

        for (int k = 0; k < A.outerSize(); ++k)
        {
            for (typename SparseMatrix::InnerIterator it(A,k); it; ++it)
            {
                const int i = it.row();
                const int j = it.col();
                if (i == j) {
                    // diagonal = D_ii
                    trips.emplace_back(i, j, d(i));
                } else if (i > j) {
                    // strict lower: add ω * A_ij  (A_ij < 0 for Laplacian)
                    trips.emplace_back(i, j, Scalar(omega) * it.value());
                }
                // do nothing for upper entries; we’ll use DL.transpose() for DU
            }
        }
        DL.setFromTriplets(trips.begin(), trips.end());
        DL.makeCompressed();
    }

    // Apply z = M^{-1} r = ω(2-ω) * (DU)^{-1} * D * (DL)^{-1} * r
    // where DU = DL.transpose() = D + ω * strictUpper(A)
    Vector apply(const Vector& r) const
    {
        // Forward solve (DL) y = r
        Vector y = DL.template triangularView<Eigen::Lower>().solve(r);

        // Multiply by D (NOT D^{-1})
        Vector z_mid = D * y;

        // Backward solve (DU) z = z_mid  with DU = DL^T, use Upper view
        Vector z = DL.transpose()
                       .template triangularView<Eigen::Upper>()
                       .solve(z_mid);

        // Final scaling
        return Scalar(omega * (2.0 - omega)) * z;
    }
};


template<typename Vector, typename Scalar = typename Vector::Scalar>
struct IncompleteCholeskyPreconditioner
{
    using RowMajorSparse = Eigen::SparseMatrix<Scalar, Eigen::RowMajor>;
    using ColMajorSparse = Eigen::SparseMatrix<Scalar>;
    using Triplet        = Eigen::Triplet<Scalar>;

    RowMajorSparse L;                       // lower-triangular factor (row-major)
    std::vector<std::vector<int>> colRows;  // rows i>j that contain (i,j) in pattern
    Scalar damping = Scalar(1e-12);

    // dot( row i, row j ) restricted to columns < j
    static Scalar dot_rows_lt(const RowMajorSparse& M, int i, int j) {
        Scalar s = Scalar(0);
        typename RowMajorSparse::InnerIterator it_i(M, i), it_j(M, j);
        // merge the two sorted column lists while both < j
        while (it_i && it_j && it_i.col() < j && it_j.col() < j) {
            const int ci = it_i.col(), cj = it_j.col();
            if      (ci == cj) { s += it_i.value() * it_j.value(); ++it_i; ++it_j; }
            else if (ci  < cj) { ++it_i; }
            else               { ++it_j; }
        }
        return s;
    }

    IncompleteCholeskyPreconditioner(const ColMajorSparse& A_in,
                                     Scalar damping_eps = Scalar(1e-12))
    : damping(damping_eps)
    {
        const int n = static_cast<int>(A_in.rows());
        L.resize(n, n);

        // copy diagonal + strictly lower from A into row-major L
        std::vector<Triplet> trips;
        trips.reserve(A_in.nonZeros());
        for (int k = 0; k < A_in.outerSize(); ++k) {
            for (typename ColMajorSparse::InnerIterator it(A_in, k); it; ++it) {
                const int i = it.row(), j = it.col();
                if (i >= j) trips.emplace_back(i, j, Scalar(it.value()));
            }
        }
        L.setFromTriplets(trips.begin(), trips.end());
        L.makeCompressed(); // sorted column indices per row

        // precompute rows that appear in each lower column j
        colRows.assign(n, {});
        for (int i = 0; i < n; ++i) {
            for (typename RowMajorSparse::InnerIterator it(L, i); it; ++it) {
                const int j = it.col();
                if (i > j) colRows[j].push_back(i);
            }
        }

        // -------------------- IC(0) factorization --------------------
        for (int j = 0; j < n; ++j) {
            // 1) diagonal update
            const Scalar Ajj0 = L.coeff(j, j);   // original a_jj sitting in L
            Scalar s = Ajj0;
            for (typename RowMajorSparse::InnerIterator it(L, j); it; ++it) {
                const int p = it.col();
                if (p >= j) break;               // only p<j
                const Scalar Ljp = it.value();
                s -= Ljp * Ljp;
            }
            if (s <= Scalar(0)) {                // simple stabilization
                const Scalar base = std::abs(Ajj0);
                s = std::max(base * damping, damping);
            }
            const Scalar Ljj = std::sqrt(s);
            L.coeffRef(j, j) = Ljj;

            // 2) below-diagonal entries in column j that exist in the pattern
            for (int i : colRows[j]) {
                const Scalar Aij0 = L.coeff(i, j);                // original a_ij
                const Scalar dot  = dot_rows_lt(L, i, j);         // ∑_{p<j} L(i,p)L(j,p)
                L.coeffRef(i, j)  = (Aij0 - dot) / Ljj;
            }
        }

        L.makeCompressed();
    }

    // z = M^{-1} r  with M = L L^T
    Vector apply(const Vector& r) const {
        Vector y = L.template triangularView<Eigen::Lower>().solve(r);
        return L.transpose().template triangularView<Eigen::Upper>().solve(y);
    }
};


using SparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using Vector       = Eigen::VectorXd;

struct IC0
{
    Eigen::IncompleteCholesky<double> ic;

    void compute (const SparseMatrix& A)
    {
        ic.setInitialShift(1e-10);
        ic.compute(A);
        if (ic.info() != Eigen::Success) 
        {
            throw std::runtime_error("IC(0) factorisation failed");
        }
    }

    // z = M^{-1} r
    void apply(const Vector& r, Vector& z) const {
        z = ic.solve(r);                  // forward+backward triangular solves
    }
};

#endif // PRECONDITIONERS_HPP