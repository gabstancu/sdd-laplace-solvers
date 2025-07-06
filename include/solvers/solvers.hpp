#ifndef SOLVERS_HPP
#define SOLVERS_HPP

#include <Eigen/Eigen>
#include "basis/LinearSystem.hpp"

#define DEFAULT_TOL 1e-4
#define MAX_ITERS 1500 

#include "solvers/ConjugateGradient.hpp"
#include "solvers/PCG.hpp"
#include "solvers/Jacobi.hpp"
#include "solvers/Multigrid.hpp"
#include "solvers/GaussSeidel.hpp"
#include "solvers/SOR.hpp"

#endif // SOLVERS_HPP