// #include "solvers.hpp"
// #include "PDE.hpp"
#include "Eigen/Eigen"
#include <string>
#include <ginac/ginac.h>
#include <gnuplot-iostream.h>
#include <cstdlib> // for std::getenv
#include "basis/basis.hpp"
#include "solvers/solvers.hpp"


std::string get_home_directory () 
{
    const char* home = std::getenv("HOME");
    return home ? std::string(home) : "";
}


int main ()
{   
    // int GRID_SIZE = 48;

    // std::vector<GiNaC::symbol> variables;
    // GiNaC::symbol x_("x"), y_("y"); 
    // variables.push_back(x_); variables.push_back(y_);

    // GiNaC::ex  top   = GiNaC::sin(GiNaC::Pi * x_);
    // GiNaC::ex bottom = GiNaC::pow(x_, 2) + GiNaC::pow(y_, 2);
    // double left  = 4;
    // double right = 2;
    
    // BoundaryConditions bc;
    // bc.top = top; 
    // bc.bottom = bottom;
    // bc.left = left;
    // bc.right = right;

    // Eigen::MatrixXd grid = initialise_grid(GRID_SIZE, bc, variables);
    // std::cout << grid << "\n\n";

    // Eigen::MatrixXd A; Eigen::VectorXd b; Eigen::VectorXd u;

    // construct_system(grid, GRID_SIZE, A, b, u);

    // Eigen::VectorXd u_0 = u; // initial guess
    // u = Conjugate_Gradient(A, b, u_0);
    // std::cout << "A:\n" << A << "\n\n";
    // std::cout << "b:\n" << b << "\n\n";
    // std::cout << "u:\n" << u << "\n\n";

    // fill_grid(grid, u);
    // std::cout << "grid:\n" << grid << "\n\n";

    // for (int dim = 8; dim <= 96; dim+=8)
    // {
    //     std::cout << "dim: " << dim << '\n';
    // }

    

    /* ---------------------- test solvers ---------------------- */
    int N = 5;

    Eigen::MatrixXd A(N, N);
    A <<  2, -1,  0,  0,  0,
         -1,  2, -1,  0,  0, 
          0, -1,  2, -1,  0, 
          0,  0, -1,  2, -1, 
          0,  0,  0, -1,  2;

    Eigen::VectorXd x_exact(N);
    x_exact << 1, 2, 3, 4, 5;

    Eigen::VectorXd b = A * x_exact;
    Eigen::VectorXd u_0 = Eigen::VectorXd::Zero(N); // initial guess
    LinearSystem<Eigen::MatrixXd, Eigen::VectorXd> system{A, b, u_0};

    DiagonalPreconditioner<Eigen::MatrixXd, Eigen::VectorXd> precon(A);
    PCG<Eigen::MatrixXd, Eigen::VectorXd, decltype(precon)> solver(precon);
    system.solve(solver);
    std::cout << system.u << '\n';

    // Eigen::VectorXd u = Conjugate_Gradient(A, b, u_0);
    // std::cout << u << '\n';

    /* test with other preconditioners */
    // Eigen::MatrixXd precon = A.diagonal().asDiagonal().inverse();
    // std::cout << precon << '\n';
    // Eigen::VectorXd u = Preconditioned_Conjugate_Gradient(A, b, u_0, precon); 
    // std::cout << u << '\n';

    // Eigen::VectorXd u = Jacobi(A, b, u_0); 
    // std::cout << "u\n" << u << "\n\n";

    // Eigen::VectorXd u = Gauss_Seidel(A, b, u_0); 
    // std::cout << "u\n" << u << "\n\n";

    // Eigen::VectorXd u = SOR(A, b, u_0, 1.05); 
    // std::cout << "u\n" << u << "\n\n";

    return 0;
}