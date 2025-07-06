#include "Eigen/Eigen"
#include <string>
#include <ginac/ginac.h>
#include <gnuplot-iostream.h>
#include <cstdlib> // for std::getenv
#include "basis/basis.hpp"
#include "solvers/solvers.hpp"
#include "Laplace2D.hpp"


std::string get_home_directory () 
{
    const char* home = std::getenv("HOME");
    return home ? std::string(home) : "";
}


int main ()
{   
    int GRID_SIZE = 8;

    std::vector<GiNaC::symbol> variables;
    GiNaC::symbol x_("x"), y_("y"); 
    variables.push_back(x_); variables.push_back(y_);

    GiNaC::ex top    = GiNaC::sin(GiNaC::Pi * x_);
    GiNaC::ex bottom = GiNaC::pow(x_, 2) + GiNaC::pow(y_, 2);
    double left  = 4;
    double right = 2;

    Laplace2D<Eigen::MatrixXd, Eigen::VectorXd> laplace(GRID_SIZE, variables);
    laplace.bc.top.expr = top;
    laplace.bc.bottom.expr = bottom;
    laplace.bc.left.expr = left;
    laplace.bc.right.expr = right;
    laplace.initialise_grid();
    std::cout << laplace.grid << '\n';

    LinearSystem<Eigen::MatrixXd, Eigen::VectorXd> system = laplace.construct_system();
    DiagonalPreconditioner<Eigen::MatrixXd, Eigen::VectorXd> precon(system.A);
    PCG<Eigen::MatrixXd, Eigen::VectorXd, decltype(precon)> solver(precon);
    system.solve(solver);
    laplace.fill_grid(system.u);

    std::cout << system.u << '\n';
    std::cout << laplace.grid << '\n';

    // for (int dim = 8; dim <= 96; dim+=8)
    // {
    //     std::cout << "dim: " << dim << '\n';
    // }


    /* ---------------------- test solvers ---------------------- */
    // int N = 5;

    // Eigen::MatrixXd A(N, N);
    // A <<  2, -1,  0,  0,  0,
    //      -1,  2, -1,  0,  0, 
    //       0, -1,  2, -1,  0, 
    //       0,  0, -1,  2, -1, 
    //       0,  0,  0, -1,  2;

    // Eigen::VectorXd x_exact(N);
    // x_exact << 1, 2, 3, 4, 5;

    // Eigen::VectorXd  b  = A * x_exact;
    // Eigen::VectorXd u_0 = Eigen::VectorXd::Zero(N); // initial guess
    // LinearSystem<Eigen::MatrixXd, Eigen::VectorXd> system{A, b, u_0};

    // // DiagonalPreconditioner<Eigen::MatrixXd, Eigen::VectorXd> precon(A);
    // // PCG<Eigen::MatrixXd, Eigen::VectorXd, decltype(precon)> solver(precon);
    // // system.solve(solver);
    // // std::cout << system.u << '\n';

    // ConjugateGradient<Eigen::MatrixXd, Eigen::VectorXd> CG;
    // system.solve(CG);
    // std::cout << system.u << '\n';

    return 0;
}