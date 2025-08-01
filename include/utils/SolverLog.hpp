#ifndef SOLVER_LOG_HPP
#define SOLVER_LOG_HPP

#include <string>
#include <chrono>
#include <fstream>
#include <iostream>
#include "helper.hpp"
template<typename Vector>
struct SolverLog
{   
    int    num_of_iterations = 0;
    int    max_iterations    = 0;
    double tolerance;
    int    converged         = 0;   
    int    timed_out         = 0; 
    int    system_dim;     

    std::vector<double>           res_per_iteration  = {};
    std::vector<double>           time_per_iteration = {}; 
    std::chrono::duration<double> time_elapsed{0};
    std::chrono::duration<double> precon_init_time{0};
    Vector                        final_solution;
    Vector                        direct_solution;

    std::string solver_name;
    std::string precon_name;


    void print ()
    {
        if (!solver_name.empty())
            std::cout << "Solver: " << solver_name << '\n';
        if (!precon_name.empty())
            std::cout << "Precon.: " << precon_name << '\n';
        std::cout << "Dim.: " << system_dim << '\n';
        std::cout << "Iterations performed: " << num_of_iterations << '\n';
        std::cout << "Converged: " <<  converged << '\n';
        std::cout << "Timed out: " <<  timed_out << '\n';
        std::cout << "Time elapsed: " << time_elapsed.count() << " seconds\n";
    }

    void log_to_file ()
    {   
        std::string log_path = get_current_working_directory() + "/logs/";

        if (std::filesystem::create_directories(log_path))
        {
            std::cout << "Directory created: " << log_path << '\n';
        }
        else
        {
            std::cout << "Directory already exists or failed to create.\n";
        }

        std::string filename = log_path + solver_name + "_";
        if (!precon_name.empty())
            filename += precon_name + "_";
        filename += std::to_string(system_dim) + ".txt";

        std::ofstream file(filename);

        if (!solver_name.empty())
            file << "C Solver: " << solver_name << '\n';
        if (!precon_name.empty())
            file << "C Precon.: " << precon_name << '\n';
        file << "C Variables: " << system_dim << '\n';
        file << "C Iterations performed: " << this->num_of_iterations << '\n';
        file << "C Time elapsed: " << this->time_elapsed.count() << '\n';
        file << "C Converged: " << this->converged << '\n';
        file << "Timed out: " <<  timed_out << '\n';
        file << "C Residual per iteration: ";
        for (double res : this->res_per_iteration)
        {
            file << res << " ";
        }
        file << '\n';
        file << "C Final solution: ";
        for (double u_ : this->final_solution)
        {
            file << u_ << " ";
        }
        file << '\n';
        std::cout << "Logs saved to: " << filename << '\n';

    }
    
    void log_to_file (std::string filename)
    {   
        std::ofstream file(filename);

        if (!solver_name.empty())
            file << "C Solver: " << solver_name << '\n';
        if (!precon_name.empty())
        {
            file << "C Precon.: " << precon_name << '\n';
            if (precon_init_time.count())
                file << "C Precon. init. time: " << precon_init_time.count() << '\n';
        }
        file << "C Variables: " << system_dim << '\n';
        file << "C Iterations performed: " << this->num_of_iterations << '\n';
        file << "C Time elapsed: " << this->time_elapsed.count() << '\n';
        file << "C Converged: " << this->converged << '\n';
        file << "Timed out: " <<  timed_out << '\n';
        file << "C Residual per iteration: ";
        for (double res : this->res_per_iteration)
        {
            file << res << " ";
        }
        file << '\n';
        file << "C Final solution: ";
        for (double u_ : this->final_solution)
        {
            file << u_ << " ";
        }
        file << '\n';
        std::cout << "Logs saved to: " << filename << '\n';
    }
};

#endif // SOLVER_LOG_HPP