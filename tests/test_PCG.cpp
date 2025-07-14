#include <iostream>
#include <Eigen/Dense>

#include "solvers/PCG.hpp"
#include "basis/LinearSystem.hpp"
#include "Laplace2D.hpp"
#include "basis/Preconditioners.hpp"

using SparseMatrix = Eigen::SparseMatrix<double>;
using Matrix       = Eigen::MatrixXd;
using Vector       = Eigen::VectorXd;

void test_PCG_direct (const std::string name)
{
    int N = 5;

    std::cout << "Testing PCG + " << name << " on a 5x5 system...\n";

    Eigen::MatrixXd A(N, N);
    A <<  2, -1,  0,  0,  0,
         -1,  2, -1,  0,  0, 
          0, -1,  2, -1,  0, 
          0,  0, -1,  2, -1, 
          0,  0,  0, -1,  2;

    Eigen::VectorXd x_exact(N);
    x_exact << 1, 2, 3, 4, 5;

    Eigen::VectorXd  b  = A * x_exact;
    Eigen::VectorXd  u = Eigen::VectorXd::Zero(N); // initial guess

    // initialising the linear system
    LinearSystem<Eigen::MatrixXd, Eigen::VectorXd> system{A, b, u};

    if (name == "Identity")
    {
        IdentityPreconditioner<Eigen::MatrixXd, Eigen::VectorXd> precon;
        PCG<Eigen::MatrixXd, Eigen::VectorXd, decltype(precon)> solver(precon, "Identity");
        system.solve(solver);

        system.solve_directly();

        solver.log.print();
        std::cout << "Approximation accuracy: " << (system.u - system._u_).norm() << '\n';

        system.reset_solution();
    }
    else if (name == "Diagonal")
    {
        DiagonalPreconditioner<Eigen::MatrixXd, Eigen::VectorXd> precon(system.A);
        PCG<Eigen::MatrixXd, Eigen::VectorXd, decltype(precon)> solver(precon, "Diagonal");
        system.solve(solver);

        system.solve_directly();

        solver.log.print();
        std::cout << "Approximation accuracy: " << (system.u - system._u_).norm() << '\n';

        system.reset_solution();
    }
    else if (name == "SSOR")
    {
        double omega_ = system.calc_omega_();
        SSORPreconditioner<Eigen::MatrixXd, Eigen::VectorXd> precon(system.A, omega_);
        PCG<Eigen::MatrixXd, Eigen::VectorXd, decltype(precon)> solver(precon, "SSOR");
        system.solve(solver);

        system.solve_directly();

        solver.log.print();
        std::cout << "Approximation accuracy: " << (system.u - system._u_).norm() << '\n';

        system.reset_solution();
    }
    else if (name == "IncompleteCholesky")
    {
        Eigen::SparseMatrix<double> A_sparse = system.A.sparseView();
        auto precon = IncompleteCholeskyPreconditioner<SparseMatrix, Vector>(A_sparse);
        PCG<Eigen::MatrixXd, Eigen::VectorXd, decltype(precon)> solver(precon, "IncompleteCholesky");
        system.solve(solver);

        system.solve_directly();

        solver.log.print();
        std::cout << "Approximation accuracy: " << (system.u - system._u_).norm() << '\n';

        system.reset_solution();
    }
}


void test_PCG_Laplace2D (const std::string name)
{
    int GRID_SIZE = 10;

    std::cout << "Testing PCG + " << name << " on a 10x10 grid...\n";

    /* symbolic setup */
    std::vector<GiNaC::symbol> variables;
    GiNaC::symbol              x_("x"), y_("y"); 
    variables.push_back(x_); 
    variables.push_back(y_);

    /* domain (unit rectangle by default) */
    std::pair<std::pair<double, double>, std::pair<double, double>> domain;
    domain = std::make_pair(std::make_pair(0.0, 1.0), std::make_pair(0.0, 1.0));

    /* init. boundary conditions */
    GiNaC::ex top    = GiNaC::sin(GiNaC::Pi * x_);
    GiNaC::ex bottom = GiNaC::pow(x_, 2) + GiNaC::pow(y_, 2);
    double left      = 4;
    double right     = 5;

    Laplace2D<Matrix, Vector> laplace2D(GRID_SIZE, variables, domain);
    laplace2D.bc.top.expr    = top;
    laplace2D.bc.bottom.expr = bottom;
    laplace2D.bc.left.expr   = left;
    laplace2D.bc.right.expr  = right;
    laplace2D.initialise_grid();

    LinearSystem<Eigen::MatrixXd, Eigen::VectorXd> system = laplace2D.construct_system();

    
    if (name == "Identity")
    {
        IdentityPreconditioner<Eigen::MatrixXd, Eigen::VectorXd> precon;
        PCG<Eigen::MatrixXd, Eigen::VectorXd, decltype(precon)> solver(precon, "Identity");
        system.solve(solver);

        system.solve_directly();

        solver.log.print();
        std::cout << "Approximation accuracy: " << (system.u - system._u_).norm() << '\n';

        system.reset_solution();
    }
    else if (name == "Diagonal")
    {
        DiagonalPreconditioner<Eigen::MatrixXd, Eigen::VectorXd> precon(system.A);
        PCG<Eigen::MatrixXd, Eigen::VectorXd, decltype(precon)> solver(precon, "Diagonal");
        system.solve(solver);

        system.solve_directly();

        solver.log.print();
        std::cout << "Approximation accuracy: " << (system.u - system._u_).norm() << '\n';

        system.reset_solution();
    }
    else if (name == "SSOR")
    {
        double omega_ = system.calc_omega_();
        SSORPreconditioner<Eigen::MatrixXd, Eigen::VectorXd> precon(system.A, omega_);
        PCG<Eigen::MatrixXd, Eigen::VectorXd, decltype(precon)> solver(precon, "SSOR");
        system.solve(solver);

        system.solve_directly();

        solver.log.print();
        std::cout << "Approximation accuracy: " << (system.u - system._u_).norm() << '\n';

        system.reset_solution();
    }
    else if (name == "IncompleteCholesky")
    {   
        Eigen::SparseMatrix<double> A_sparse = system.A.sparseView();
        auto precon = IncompleteCholeskyPreconditioner<SparseMatrix, Vector>(A_sparse);
        PCG<Eigen::MatrixXd, Eigen::VectorXd, decltype(precon)> solver(precon, "IncompleteCholesky");
        system.solve(solver);

        system.solve_directly();

        solver.log.print();
        std::cout << "Approximation accuracy: " << (system.u - system._u_).norm() << '\n';

        system.reset_solution();
    }
}


int main ()
{   
    std::cout << "======================= PCD with Identity preconditioner =======================\n";
    test_PCG_direct    ("Identity");
    test_PCG_Laplace2D ("Identity");


    std::cout << "======================= PCD with Diagonal preconditioner =======================\n";
    test_PCG_direct    ("Diagonal");
    test_PCG_Laplace2D ("Diagonal");


    std::cout << "======================= PCD with SSOR preconditioner =======================\n";
    test_PCG_direct    ("SSOR");
    test_PCG_Laplace2D ("SSOR");


    std::cout << "======================= PCD with Incomplete Cholesky preconditioner =======================\n";
    test_PCG_direct    ("IncompleteCholesky");
    test_PCG_Laplace2D ("IncompleteCholesky");

    return 0;
}