#ifndef EVALUATE_HPP
#define EVALUATE_HPP

#include <cstdlib> // for std::getenv
#include <filesystem>
#include <chrono>
#include <chrono>

#define START_GRID_DIMENSION 12
#define MAX_GRID_DIMENSION 100
#define STEP_SIZE 8

inline std::string get_home_directory () 
{
    const char* home = std::getenv("HOME");
    return home ? std::string(home) : "";
}


inline std::string get_current_working_directory() {
    return std::filesystem::current_path().string();
}


template<typename System, typename Solver>
void evaluate (System& system, Solver& solver, int grid_dim)
{
    std::string log_path = get_current_working_directory() + "/results/" + solver.name +"/";
    if (std::filesystem::create_directories(log_path)) {
        std::cout << "Directory created: " << log_path << '\n';
    } else {
        std::cout << "Directory already exists or failed to create.\n";
    }
    std::string filename = solver.name + "_" + std::to_string(grid_dim) + ".txt";

    std::cout << log_path+filename << '\n';
    system.solve(solver);
    solver.log.log_to_file(log_path+filename, solver.name);
    solver.log.print();
}

#endif // EVALUATE_HPP