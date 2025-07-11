#include "Eigen/Eigen"
#include <string>
#include <ginac/ginac.h>
#include <gnuplot-iostream.h>
#include "basis/basis.hpp"
#include "solvers/solvers.hpp"
#include "Laplace2D.hpp"
#include "evaluate.hpp"
#include <iostream>
 

int main ()
{   
    std::vector<GiNaC::symbol> variables;
    GiNaC::symbol x_("x"), y_("y"); 
    variables.push_back(x_); variables.push_back(y_);

    GiNaC::ex top    = GiNaC::sin(GiNaC::Pi * x_);
    GiNaC::ex bottom = GiNaC::pow(x_, 2) + GiNaC::pow(y_, 2);
    double left      = 4;
    double right     = 5;

    for (int dim = START_GRID_DIMENSION; dim <= START_GRID_DIMENSION; dim+=STEP_SIZE)
    {
        Laplace2D<Eigen::MatrixXd, Eigen::VectorXd> laplace(dim, variables);
        laplace.bc.top.expr    = top;
        laplace.bc.bottom.expr = bottom;
        laplace.bc.left.expr   = left;
        laplace.bc.right.expr  = right;
        laplace.initialise_grid();

        LinearSystem<Eigen::MatrixXd, Eigen::VectorXd> system = laplace.construct_system();

        SSORPreconditioner<Eigen::MatrixXd, Eigen::VectorXd>    precon(system.A, 1.2);
        PCG<Eigen::MatrixXd, Eigen::VectorXd, decltype(precon)> PCG(precon);
        evaluate(system, PCG, dim);
        laplace.fill_grid(system.u);
        std::string filename = PCG.name + "/grid_" + std::to_string(dim) + ".txt";
        laplace.save_grid(filename);
        
        system = laplace.construct_system();
        ConjugateGradient<Eigen::MatrixXd, Eigen::VectorXd> CG;
        evaluate(system, CG, dim);
        laplace.fill_grid(system.u);
        std::string filename = CG.name + "/grid_" + std::to_string(dim) + ".txt";
        laplace.save_grid(filename);

        system = laplace.construct_system();
        GaussSeidel<Eigen::MatrixXd, Eigen::VectorXd> GS;
        evaluate(system, GS, dim);
        laplace.fill_grid(system.u);
        std::string filename = GS.name + "/grid_" + std::to_string(dim) + ".txt";
        laplace.save_grid(filename);

        system = laplace.construct_system();
        Jacobi<Eigen::MatrixXd, Eigen::VectorXd> Jacobi;
        evaluate(system, Jacobi, dim);
        laplace.fill_grid(system.u);
        std::string filename = Jacobi.name + "/grid_" + std::to_string(dim) + ".txt";
        laplace.save_grid(filename);

        system = laplace.construct_system();
        // TODO: write code to calculate optimal omega_
        SOR<Eigen::MatrixXd, Eigen::VectorXd> SOR(1.2);
        evaluate(system, SOR, dim);
        laplace.fill_grid(system.u);
        std::string filename = SOR.name + "/grid_" + std::to_string(dim) + ".txt";
        laplace.save_grid(filename);
    }


    /* ---------------------- test ---------------------- */
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