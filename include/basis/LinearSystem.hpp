#ifndef LINEAR_SYSTEM_HPP
#define LINEAR_SYSTEM_HPP

#include <chrono>

template<typename Matrix, typename Vector>
struct LinearSystem
{
    Matrix A;
    Vector b;
    Vector u; // solution vector (initial guess)
    int    N = A.rows();

    template<typename Solver>
    void solve(Solver& solver)
    {   
        auto start = std::chrono::high_resolution_clock::now();
        solver.solve(*this);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        solver.log.time_elapsed = elapsed;
    }

    void print ()
    {
        
    }
};


#endif // LINEAR_SYSTEM_HPP