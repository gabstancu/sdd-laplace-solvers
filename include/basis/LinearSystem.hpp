#ifndef LINEAR_SYSTEM_HPP
#define LINEAR_SYSTEM_HPP

template<typename Matrix, typename Vector>
struct LinearSystem
{
    Matrix A;
    Vector b;
    Vector u; // solution vector (initial guess)
    int N = A.rows();

    template<typename Solver>
    void solve(Solver& solver)
    {
        solver.solve(*this);
    }

    void print ()
    {
        
    }
};


#endif // LINEAR_SYSTEM_HPP