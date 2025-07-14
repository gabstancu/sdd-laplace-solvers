#include <iostream>
#include <Eigen/Dense>
#include <cassert>
#include "solvers/Jacobi.hpp"
#include "basis/LinearSystem.hpp"
#include "Laplace2D.hpp"
#include <ginac/ginac.h>

void test_Jacobi_direct ()
{
    int N = 5;

    std::cout << "Testing Jacobi on a 5x5 system...\n";

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

    Jacobi<Eigen::MatrixXd, Eigen::VectorXd> solver;
    system.solve(solver);

    system.solve_directly();
    
    solver.log.print();
    std::cout << "Approximation accuracy: " << (system.u - system._u_).norm() << '\n';

}


void test_Jacobi_Laplace2D ()
{
    int GRID_SIZE = 10;

    std::cout << "Testing Jacobi on a 10x10 grid...\n";

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

    Laplace2D<Eigen::MatrixXd, Eigen::VectorXd> laplace2D(GRID_SIZE, variables, domain);
    laplace2D.bc.top.expr    = top;
    laplace2D.bc.bottom.expr = bottom;
    laplace2D.bc.left.expr   = left;
    laplace2D.bc.right.expr  = right;
    laplace2D.initialise_grid();

    auto system = laplace2D.construct_system();

    Jacobi<Eigen::MatrixXd, Eigen::VectorXd> solver;
    system.solve(solver);

    system.solve_directly();
    
    solver.log.print();
    std::cout << "Approximation accuracy: " << (system.u - system._u_).norm() << '\n';
}


int main ()
{
    test_Jacobi_direct      ();
    test_Jacobi_Laplace2D   ();
    
    return 0;
}