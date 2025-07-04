#ifndef SOLVERS_HPP
#define SOLVERS_HPP

#include <Eigen/Eigen>
#include <iostream>

#define TOL 1e-5
#define MAXITERS 400

Eigen::MatrixXd Conjugate_Gradient (Eigen::MatrixXd A, Eigen::VectorXd b, Eigen::VectorXd u_0);

Eigen::MatrixXd Preconditioned_Conjugate_Gradient (Eigen::MatrixXd A, Eigen::VectorXd b, Eigen::VectorXd u_0, Eigen::MatrixXd precon);

Eigen::MatrixXd Jacobi (Eigen::MatrixXd A, Eigen::VectorXd b, Eigen::VectorXd u_0);

Eigen::MatrixXd Gauss_Seidel (Eigen::MatrixXd A, Eigen::VectorXd b, Eigen::VectorXd u_0);

Eigen::MatrixXd SOR (Eigen::MatrixXd A, Eigen::VectorXd b, Eigen::VectorXd u_0, double omega); // Successive Over-Relaxation (SOR)

Eigen::MatrixXd Multigrid (Eigen::MatrixXd A, Eigen::VectorXd b, Eigen::VectorXd u_0);

Eigen::MatrixXd GMRM ();

#endif // SOLVERS_HPP