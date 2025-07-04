#include "PDE.hpp"
#include <fstream>

void set_boundary_value ( Eigen::MatrixXd& grid, 
                          double h,
                          BoundaryConditions bc,
                          std::pair<int, int> coordinates, 
                          std::optional<BoundaryConditions::BoundaryType> boundary, 
                          std::vector<GiNaC::symbol> vars )
{   
    int i = coordinates.first, j = coordinates.second;

    if (bc.is_set(boundary))
    {
        if (bc.is_constant(boundary))
        {
            grid(i, j) = std::get<double>(*boundary); 
        }
        else if (bc.is_expression(boundary))
        {
            GiNaC::exmap m; 
            m[vars[0]] = j * h; m[vars[1]] = i * h;
            GiNaC::ex expr = std::get<GiNaC::ex>(*boundary);
            GiNaC::ex evaluated_boundary = expr.subs(m).evalf();
            double val = GiNaC::ex_to<GiNaC::numeric>(evaluated_boundary).to_double();
            grid(i, j) = val;
        }  
    }
}


Eigen::MatrixXd initialise_grid ( int GRID_SIZE,  
                                  BoundaryConditions bc, 
                                  std::vector<GiNaC::symbol> vars )
{   
    int N = GRID_SIZE - 2; // inner grid dimension (A.dim)
    double h = 1.0 / (N + 1.0); // grid spacing
    // double x_, y_; // (cartesian) coordinates of grid points

    std::cout << "Initialising grid:\n";
    std::cout << "h = " << h << '\n';
    std::cout << "GRID SIZE = " << GRID_SIZE << '\n';

    Eigen::MatrixXd grid(GRID_SIZE, GRID_SIZE);
    GiNaC::exmap m;
    std::optional<BoundaryConditions::BoundaryType> boundary;

    // fill grid (set interior points to 0 - no solution)
    for (int i = 0; i < GRID_SIZE; i++)
    {
        for (int j = 0; j < GRID_SIZE; j++)
        {   
            grid(i, j) = 0.0;

            // top boundary
            if (i == 0)
            {   
                boundary = bc.top;
                set_boundary_value(grid, h, bc, std::make_pair(i, j), boundary, vars);
            }

            // bottom boundary
            if (i == GRID_SIZE - 1)
            {   
                boundary = bc.bottom;
                set_boundary_value(grid, h, bc, std::make_pair(i, j), boundary, vars);
            } 

            // left boundary
            if (j == 0 && (i > 0 && i < GRID_SIZE - 1)) // keep top's value
            {   
                boundary = bc.left;
                set_boundary_value(grid, h, bc, std::make_pair(i, j), boundary, vars);
            }

            // right boundary
            if (j == GRID_SIZE - 1 && (i > 0 && i < GRID_SIZE - 1))
            {   
                boundary = bc.right;
                set_boundary_value(grid, h, bc, std::make_pair(i, j), boundary, vars);
            }
        }
    }

    return grid;
}


void construct_system ( Eigen::MatrixXd grid, int GRID_SIZE, 
                        Eigen::MatrixXd& A, 
                        Eigen::VectorXd& b, 
                        Eigen::VectorXd& u )
{   
    int N = GRID_SIZE - 2; // inner grid dimension
    int dim = int(std::pow(N, 2));
    
    A = Eigen::MatrixXd::Zero(dim, dim);
    b = Eigen::VectorXd::Zero(dim);
    u = Eigen::VectorXd::Zero(dim);

    /* -4 u_{i, j} + u_{i+1, j} + u_{i-1, j} + u_{i, j+1} + u_{i, j-1} = b_{k} */
    int r = 0, k;
    for (int i = 1; i <= N; i++)
    {
        for (int j = 1; j <= N; j++)
        {   
            k = (i - 1) * N + (j - 1);
            // std::cout << k << "\n";
            // std::cout << "("<< i << " " << j << ")" << "\n";
            
            A(k, k) = -4;

            // bottom neighbor
            if (i == N)
            {
                b(k) -= grid(i + 1, j);
            }
            else {
                A(k, k + N) = 1.0;
            }

            // top neighbor
            if (i == 1)
            {
                b(k) -= grid(i - 1, j);
            }
            else {
                A(k, k - N) = 1.0;
            }

            // right neighbor
            if (j == N)
            {
                b(k) -= grid(i, j + 1);
            }
            else {
                A(k, k + 1) = 1.0;
            }

            // left neighbor
            if (j - 1 == 0)
            {
                b(k) -= grid(i, j - 1);
            }
            else {
                A(k, k - 1) = 1.0;
            }
        }
    }
}


void fill_grid (Eigen::MatrixXd& grid, Eigen::VectorXd u)
{   
    int ROWS = grid.rows();
    int COLS = grid.cols();

    int K = u.size();
    int i, j;

    for (int k = 0; k < K; k++)
    {
        i = k / (ROWS-2) + 1;
        j = k % (COLS-2) + 1;
        grid(i, j) = u(k);
    }
}


void save (Eigen::MatrixXd grid, std::string filename)
{
    int ROWS = grid.rows();
    int COLS = grid.cols(); 

    std::ofstream file(filename);

    for (int i = 0; i < ROWS; i++)
    {
        for (int j = 0; j < COLS; j++)
        {
            file << grid(i, j) << " ";
        }
        file << '\n';
    }
}