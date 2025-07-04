#ifndef LINEAR_SYSTEM_HPP
#define LINEAR_SYSTEM_HPP

template<typename Matrix, typename Vector>
struct LinearSystem
{
    Matrix A;
    Vector b;
    Vector u; // solution vector (initial guess)

    template<typename Solver>
    void solve(const Solver& solver)
    {
        solver.solve(*this);
    }
};


#endif // LINEAR_SYSTEM_HPP