#ifndef EVALUATE_HPP
#define EVALUATE_HPP

#include <cstdlib> // for std::getenv
#include <filesystem>
#include <chrono>
#include "utils/helper.hpp"

#define START_GRID_DIMENSION 12
#define MAX_GRID_DIMENSION   100
#define STEP_SIZE            8


template<typename System, typename Solver>
void evaluate (System& system, Solver& solver, int grid_dim)
{
    std::string log_path = get_current_working_directory() + "/logs/" + solver.name +"/";

    if (std::filesystem::create_directories(log_path)) 
    {
        std::cout << "Directory created: " << log_path << '\n';
    } 
    else 
    {
        std::cout << "Directory already exists or failed to create.\n";
    }
    std::string filename = solver.name + "_" + std::to_string(grid_dim) + ".txt";

    system.solve(solver);
    solver.log.print();
    solver.log.log_to_file(log_path+filename, solver.name);
}

#endif // EVALUATE_HPP