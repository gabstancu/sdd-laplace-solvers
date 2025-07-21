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
    Vector  u; 
    Vector _u_; // direct solution
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

    void solve_directly (bool log)
    {
        // this->_u_ = A.llt().solve(this->b);
        this->_u_ = A.lu().solve(this->b);

        if (!log)
            return;

        #ifdef TESTING_MODE
        std::cout << "[TEST] Skipping file logging.\n";
        #endif
        
        #ifndef TESTING_MODE
        std::string log_path = get_current_working_directory() + "/logs/direct/";
        std::filesystem::create_directories(log_path);
        std::string filename = log_path + "inst_" +std::to_string(this->A.rows()) + ".txt";

        std::ofstream file(filename);

        file << "Variables: " << this->A.rows() << '\n';
        file << "omega_: " << this->omega_ << '\n';
        file << "Solution vector: ";
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
        this->omega_ = 2.0 / (1 + std::sin((M_PI) / this->N ));
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