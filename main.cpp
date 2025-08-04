#include "Eigen/Eigen"
#include <string>
#include <ginac/ginac.h>
#include "basis/basis.hpp"
#include "solvers/solvers.hpp"
#include "Laplace2D.hpp"
#include "utils/helper.hpp"
#include <iostream>

#define START_GRID_DIMENSION 12
#define MAX_GRID_DIMENSION   100
#define STEP_SIZE            8
 

/* capstone project part */
void evaluate_loop ()
{
    std::vector<GiNaC::symbol> variables;
    GiNaC::symbol              x_("x"), y_("y"); 

    variables.push_back(x_); 
    variables.push_back(y_);
    std::pair<std::pair<double, double>, std::pair<double, double>> domain;
    domain = std::make_pair(std::make_pair(0.0, 1.0), std::make_pair(0.0, 1.0));

    GiNaC::ex top    = 4 * GiNaC::sin(3 * GiNaC::Pi * x_);
    GiNaC::ex bottom = GiNaC::sin(GiNaC::Pi * x_);
    double left      = 0;
    double right     = 0;

    GiNaC::ex analytical_solution = (1 / GiNaC::sinh(GiNaC::Pi)) * GiNaC::sinh(GiNaC::Pi * (1 - y_)) * GiNaC::sin(GiNaC::Pi * x_)
                                  + (4 / GiNaC::sinh(3 * GiNaC::Pi)) * GiNaC::sinh(3 * GiNaC::Pi * y_) * GiNaC::sin(3 * GiNaC::Pi * x_);

    for (int dim = START_GRID_DIMENSION; dim <= MAX_GRID_DIMENSION; dim+=STEP_SIZE)
    {   

        std::cout << "===================================== GRID DIMENSION " << dim << " =====================================\n";
        Laplace2D<Eigen::MatrixXd, Eigen::VectorXd> laplace(dim, variables, domain);
        laplace.bc.top.expr           = top;
        laplace.bc.bottom.expr        = bottom;
        laplace.bc.left.expr          = left;
        laplace.bc.right.expr         = right;
        laplace.analytical_expression = analytical_solution;
        laplace.initialise_grid();
        laplace.evaluate_analytical_solution();

        std::string gridfile = "analytical_"+std::to_string(dim)+".dat";
        laplace.save_grid(gridfile, true, "capstone/grids/");


        LinearSystem<Eigen::MatrixXd, Eigen::VectorXd> system = laplace.construct_system();
        double omega_ = system.calc_omega_();
        std::cout << "omega_: " << omega_ << '\n';

        ConjugateGradient<Eigen::MatrixXd, Eigen::VectorXd> CG;
        std::cout << "----------------------- " << CG.name << " -----------------------\n";

        std::string log_path = get_current_working_directory() + "/capstone/eval_results/" + CG.name +"/";
        std::filesystem::create_directories(log_path);
        std::string filename = CG.name + "_approx_" + std::to_string(dim) + ".txt";

        system.solve(CG);
        CG.log.print();
        CG.log.log_to_file(log_path+filename);

        laplace.fill_grid(system.u);
        system.reset_solution();


        SSORPreconditioner<Eigen::MatrixXd, Eigen::VectorXd>    precon(system.A, omega_);
        PCG<Eigen::MatrixXd, Eigen::VectorXd, decltype(precon)> PCG(precon, "SSOR");
        std::cout << "----------------------- " << PCG.name << " -----------------------\n";

        log_path = get_current_working_directory() + "/capstone/eval_results/" + PCG.name +"/";
        std::filesystem::create_directories(log_path);
        filename = PCG.name + "_approx_" + std::to_string(dim) + ".txt";

        system.solve(PCG);
        PCG.log.print();
        PCG.log.log_to_file(log_path+filename);

        laplace.fill_grid(system.u);
        system.reset_solution();

        GaussSeidel<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::VectorXd> GS(system);
        std::cout << "----------------------- " << GS.name << " -----------------------\n";

        log_path = get_current_working_directory() + "/capstone/eval_results/" + GS.name +"/";
        std::filesystem::create_directories(log_path);
        filename = GS.name + "_approx_" + std::to_string(dim) + ".txt";

        system.solve(GS);
        GS.log.print();
        GS.log.log_to_file(log_path+filename);

        laplace.fill_grid(system.u);
        system.reset_solution();

        Jacobi<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::VectorXd> Jacobi(system);
        std::cout << "----------------------- " << Jacobi.name << " -----------------------\n";

        log_path = get_current_working_directory() + "/capstone/eval_results/" + Jacobi.name +"/";
        std::filesystem::create_directories(log_path);
        filename = Jacobi.name + "_approx_" + std::to_string(dim) + ".txt";

        system.solve(Jacobi);
        Jacobi.log.print();
        Jacobi.log.log_to_file(log_path+filename);

        laplace.fill_grid(system.u);
        system.reset_solution();

        SOR<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>, Eigen::VectorXd> SOR(system);
        std::cout << "----------------------- " << SOR.name << " -----------------------\n";

        log_path = get_current_working_directory() + "/capstone/eval_results/" + SOR.name +"/";
        std::filesystem::create_directories(log_path);
        filename = SOR.name + "_approx_" + std::to_string(dim) + ".txt";

        system.solve(SOR);
        SOR.log.print();
        SOR.log.log_to_file(log_path+filename);

        laplace.fill_grid(system.u);
        system.reset_solution();
        // break;
    }
}


void evaluate_preconditioners ()
{
    std::vector<GiNaC::symbol> variables;
    GiNaC::symbol              x_("x"), y_("y"); 

    variables.push_back(x_); 
    variables.push_back(y_);
    std::pair<std::pair<double, double>, std::pair<double, double>> domain;
    domain = std::make_pair(std::make_pair(0.0, 1.0), std::make_pair(0.0, 1.0));

    GiNaC::ex top    = 4 * GiNaC::sin(3 * GiNaC::Pi * x_);
    GiNaC::ex bottom = GiNaC::sin(GiNaC::Pi * x_);
    double left      = 0;
    double right     = 0;

    GiNaC::ex analytical_solution = (1 / GiNaC::sinh(GiNaC::Pi)) * GiNaC::sinh(GiNaC::Pi * (1 - y_)) * GiNaC::sin(GiNaC::Pi * x_)
                                  + (4 / GiNaC::sinh(3 * GiNaC::Pi)) * GiNaC::sinh(3 * GiNaC::Pi * y_) * GiNaC::sin(3 * GiNaC::Pi * x_);


    for (int dim = START_GRID_DIMENSION; dim <= MAX_GRID_DIMENSION; dim += STEP_SIZE)
    {
        std::cout << "===================================== GRID DIMENSION " << dim << " =====================================\n";
        Laplace2D<Eigen::MatrixXd, Eigen::VectorXd> laplace(dim, variables, domain);
        laplace.bc.top.expr           = top;
        laplace.bc.bottom.expr        = bottom;
        laplace.bc.left.expr          = left;
        laplace.bc.right.expr         = right;
        laplace.analytical_expression = analytical_solution;
        laplace.initialise_grid();
        laplace.evaluate_analytical_solution();

        std::string gridfile = "analytical_"+std::to_string(dim)+".dat";
        laplace.save_grid(gridfile, true, "capstone/grids/");

        LinearSystem<Eigen::MatrixXd, Eigen::VectorXd> system = laplace.construct_system();
        double omega_ = system.calc_omega_();
        std::cout << "omega_: " << omega_ << '\n';

        std::cout << "------------- Identity -------------\n";
        auto start = std::chrono::high_resolution_clock::now();
        IdentityPreconditioner<Eigen::MatrixXd, Eigen::VectorXd>    precon_I;
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        PCG<Eigen::MatrixXd, Eigen::VectorXd, decltype(precon_I)> PCG_I(precon_I, "Identity");
        PCG_I.log.precon_init_time = elapsed;

        std::string log_path = get_current_working_directory() + "/capstone/eval_precon/" + PCG_I.precon_name +"/";
        std::filesystem::create_directories(log_path);
        std::string filename = PCG_I.name + "_approx_" + std::to_string(dim) + ".txt";

        system.solve(PCG_I);
        PCG_I.log.print();
        PCG_I.log.log_to_file(log_path+filename);

        laplace.fill_grid(system.u);
        system.reset_solution();


        std::cout << "------------- Diagonal -------------\n";
        start = std::chrono::high_resolution_clock::now();
        DiagonalPreconditioner<Eigen::MatrixXd, Eigen::VectorXd>    precon_D(system.A);
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        PCG<Eigen::MatrixXd, Eigen::VectorXd, decltype(precon_D)> PCG_D(precon_D, "Diagonal");
        PCG_D.log.precon_init_time = elapsed;

        log_path = get_current_working_directory() + "/capstone/eval_precon/" + PCG_D.precon_name +"/";
        std::filesystem::create_directories(log_path);
        filename = PCG_D.name + "_approx_" + std::to_string(dim) + ".txt";

        system.solve(PCG_D);
        PCG_D.log.print();
        PCG_D.log.log_to_file(log_path+filename);

        laplace.fill_grid(system.u);
        system.reset_solution();


        std::cout << "------------- SSOR -------------\n";
        start = std::chrono::high_resolution_clock::now();
        SSORPreconditioner<Eigen::MatrixXd, Eigen::VectorXd>    precon_SSOR(system.A, omega_);
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        PCG<Eigen::MatrixXd, Eigen::VectorXd, decltype(precon_SSOR)> PCG_SSOR(precon_SSOR, "SSOR");
        PCG_SSOR.log.precon_init_time = elapsed;

        log_path = get_current_working_directory() + "/capstone/eval_precon/" + PCG_SSOR.precon_name +"/";
        std::filesystem::create_directories(log_path);
        filename = PCG_SSOR.name + "_approx_" + std::to_string(dim) + ".txt";

        system.solve(PCG_SSOR);
        PCG_SSOR.log.print();
        PCG_SSOR.log.log_to_file(log_path+filename);

        laplace.fill_grid(system.u);
        system.reset_solution();

        std::cout << "------------- Incomplete Cholesky -------------\n";
        Eigen::SparseMatrix<double> A_sparse = system.A.sparseView();
        start = std::chrono::high_resolution_clock::now();
        auto precon_IC = IncompleteCholeskyPreconditioner<Eigen::SparseMatrix<double>, Eigen::VectorXd>(A_sparse);
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        PCG<Eigen::MatrixXd, Eigen::VectorXd, decltype(precon_IC)> PCG_IC(precon_IC, "IncompleteCholesky");
        PCG_IC.log.precon_init_time = elapsed;

        log_path = get_current_working_directory() + "/capstone/eval_precon/" + PCG_IC.precon_name +"/";
        std::filesystem::create_directories(log_path);
        filename = PCG_IC.name + "_approx_" + std::to_string(dim) + ".txt";

        system.solve(PCG_IC);
        PCG_IC.log.print();
        PCG_IC.log.log_to_file(log_path+filename);

        laplace.fill_grid(system.u);
        system.reset_solution();
        // break;
    }
}


int main ()
{   
    evaluate_loop();
    evaluate_preconditioners();
    // std::vector<GiNaC::symbol> variables;
    // GiNaC::symbol              x_("x"), y_("y"); 
    // int GRID_SIZE = 10;

    // variables.push_back(x_); 
    // variables.push_back(y_);
    // std::pair<std::pair<double, double>, std::pair<double, double>> domain;
    // domain = std::make_pair(std::make_pair(0.0, 1.0), std::make_pair(0.0, 1.0));

    // GiNaC::ex top    = GiNaC::sin(GiNaC::Pi * x_) * GiNaC::sinh(GiNaC::Pi);
    // GiNaC::ex bottom = 0;
    // double left      = 0;
    // double right     = 0;
    // // double left      = 4;
    // // double right     = 5;

    // GiNaC::ex analytical_solution = GiNaC::sin(GiNaC::Pi * x_) * GiNaC::sinh(GiNaC::Pi * y_);

    // Laplace2D<Eigen::MatrixXd, Eigen::VectorXd> laplace(GRID_SIZE, variables, domain);
    // laplace.bc.top.expr    = top;
    // laplace.bc.bottom.expr = bottom;
    // laplace.bc.left.expr   = left;
    // laplace.bc.right.expr  = right;
    // laplace.initialise_grid();
    // laplace.analytical_expression = analytical_solution;

    // // std::cout << laplace.grid << "\n\n";

    // laplace.evaluate_analytical_solution();
    // // std::cout << laplace.analytical_solution << "\n\n";

    // LinearSystem<Eigen::MatrixXd, Eigen::VectorXd> system = laplace.construct_system();
    // double omega_ = system.calc_omega_();

    // // std::cout << laplace.grid << "\n";
    // // std::cout << "Coefficient matrix:\n" << system.A << "\n";

    // SSORPreconditioner<Eigen::MatrixXd, Eigen::VectorXd> precon(system.A, omega_);
    // PCG<Eigen::MatrixXd, Eigen::VectorXd, decltype(precon)> solver(precon, "SSOR");
    // system.solve(solver);
    // laplace.fill_grid(system.u);
    // std::cout << "Approximation:\n" << laplace.grid                << "\n\n";
    // std::cout << "Analytical:\n"    << laplace.analytical_solution << "\n\n";

    // std::cout << "Analytical / approximation diff.:\n" << laplace.analytical_solution - laplace.grid << "\n";

    // // std::string filename = solver.name + "_approx_" + std::to_string(GRID_SIZE) + ".txt";
    // // laplace.save_grid(filename, false);
    // // filename = solver.name + "_analytical_" + std::to_string(GRID_SIZE) + ".txt";
    // // laplace.save_grid(filename, true);

    // /* ------------------------- Usage example -------------------------*/
    // // int N = 5;

    // // Eigen::MatrixXd A(N, N);
    // // A <<  -2,  1,  0,  0,  0,
    // //        1, -2,  1,  0,  0, 
    // //        0,  1, -2,  1,  0, 
    // //        0,  0,  1, -2,  1, 
    // //        0,  0,  0,  1, -2;

    // // double omega_ = 1.60;

    // // Eigen::VectorXd x_exact(N);
    // // x_exact << 1, 2, 3, 4, 5;

    // // Eigen::VectorXd  b  = A * x_exact;
    // // Eigen::VectorXd u_0 = Eigen::VectorXd::Zero(N); // initial guess
    // // LinearSystem<Eigen::MatrixXd, Eigen::VectorXd> system{A, b, u_0};

    // // SSORPreconditioner<Eigen::MatrixXd, Eigen::VectorXd> precon(A, omega_);
    // // PCG<Eigen::MatrixXd, Eigen::VectorXd, decltype(precon)> solver(precon, "SSOR");
    // // system.solve(solver);
    // // std::cout << system.u << "\n\n";
    // // solver.log.log_to_file();

    // // system.reset_solution();

    // ConjugateGradient<Eigen::MatrixXd, Eigen::VectorXd> CG;
    // system.solve(CG);
    // laplace.fill_grid(system.u);
    // std::cout << "Appoximation:\n" << laplace.grid << "\n\n";
    // std::cout << "Analytical:\n" << laplace.analytical_solution << "\n\n";

    // std::cout << "Analytical / approximation diff.:\n" << laplace.analytical_solution - laplace.grid << "\n";

    // std::string filename = CG.name + "_approx_" + std::to_string(GRID_SIZE) + ".txt";
    // laplace.save_grid(filename, false);
    // filename = CG.name + "_analytical_" + std::to_string(GRID_SIZE) + ".txt";
    // laplace.save_grid(filename, true);

    // solver.log.log_to_file();

    // system.reset_solution();

    // system.solve_directly(false);

    return 0;
}