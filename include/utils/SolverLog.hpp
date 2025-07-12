#ifndef SOLVER_LOG_HPP
#define SOLVER_LOG_HPP

#include <string>
#include <chrono>
#include <fstream>
#include <iostream>
template<typename Vector>
struct SolverLog
{   
    int    num_of_iterations = 0;
    int    max_iterations    = 0;
    double tolerance;
    int    converged         = 0;

    std::vector<double>           res_per_iteration;
    std::chrono::duration<double> time_elapsed{0};
    Vector                        final_solution;

    void print(std::string solver_name = "")
    {
        if (!solver_name.empty())
            std::cout << "Solver: " << solver_name << '\n';
        std::cout << "Iterations performed: " << num_of_iterations << '\n';
        std::cout << "Converged: " <<  converged << '\n';
        std::cout << "Time elapsed: " << time_elapsed.count() << " seconds\n";
    }
    
    void log_to_file (std::string filename, std::string solver_name = "")
    {
        std::ofstream file(filename);

        if (!solver_name.empty())
            file << "C Solver: " << solver_name << '\n';
        file << "C Iterations performed: " << num_of_iterations << '\n';
        file << "C Time elapsed: " << time_elapsed.count() << '\n';
        file << "C Converged: " << converged << '\n';
        file << "C Residual per iteration: ";
        for (double res : res_per_iteration)
        {
            file << res << " ";
        }
        file << '\n';
        file << "C Final solution: ";
        for (double u_ : final_solution)
        {
            file << u_ << " ";
        }
        file << '\n';
        std::cout << "Logs saved to: " << filename << '\n';
    }
};

#endif // SOLVER_LOG_HPP