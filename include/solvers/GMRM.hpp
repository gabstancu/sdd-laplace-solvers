#ifndef GMRM_HPP
#define GMRM_HPP

template<typename Matrix, typename Vector>
struct GMRM
{
    double tol = DEFAULT_TOL;
    int max_iters = MAX_ITERS;

    void solve(LinearSystem<Matrix, Vector>& system)
    {
        
    }
};


#endif // GMRM_HPP