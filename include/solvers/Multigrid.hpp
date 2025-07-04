#ifndef MULTIGRID_HPP
#define MULTIGRID_HPP

template<typename Matrix, typename Vector>
struct Multigrid
{
    double tol = DEFAULT_TOL;
    int max_iters = MAX_ITERS;

    void solve(LinearSystem<Matrix, Vector>& system) const
    {
        
    }
};


#endif // MULTIGRID_HPP