#ifndef LINEAR_SYSTEM_HPP
#define LINEAR_SYSTEM_HPP

#include <chrono>
#include "utils/helper.hpp"
#include <fstream>

template<typename Matrix, typename Vector>
struct LinearSystem
{
    Matrix  A;
    Vector  b;
    Vector  u;  // solution vector (initial guess)
    Vector _u_; // ground truth solution
    int     N ;
    double  omega_;

    

    template<typename Solver>
    void solve(Solver& solver)
    {   
        auto start = std::chrono::high_resolution_clock::now();
        solver.solve(*this);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        solver.log.time_elapsed = elapsed;
    }

    void solve_directly ()
    {
        // this->_u_ = A.llt().solve(this->b);
        this->_u_ = A.lu().solve(this->b);
        
        #ifndef TESTING_MODE
        std::string log_path = get_current_working_directory() + "/logs/direct_solving/";
        std::filesystem::create_directories(log_path);
        std::string filename = log_path + std::to_string(this->N) + ".txt";
        std::ofstream file(filename);
        
        for (double u : this->_u_)
        {
            file << u << " ";
        }
        file << '\n';
        std::cout << "Direct solution saved to " << filename << '\n';
        #endif
    }

    double calc_omega_()
    {
        this->omega_ = 2.0 / (1 + std::sin((M_PI) / (this->N + 1)));
        this->omega_ = std::min(this->omega_, 1.9 + (0.05) / this->N);
        return this->omega_;
    }

    void reset_solution ()
    {
        this->u.setZero();
    }

    void print ()
    {
        
    }
};


#endif // LINEAR_SYSTEM_HPP