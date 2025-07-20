#ifndef LAPLACE_2D_HPP
#define LAPLACE_2D_HPP

#include "basis/BoundaryCondition.hpp"
#include "basis/LinearSystem.hpp"

template<typename Matrix, typename Vector>
struct Laplace2D
{   
    int    GRID_SIZE;
    Matrix grid;
    double h;
    std::pair<std::pair<double, double>, std::pair<double, double>> domain;
    BoundaryConditions bc;
    std::vector<GiNaC::symbol> vars;

    Laplace2D (int grid_size, 
               std::vector<GiNaC::symbol> variables, 
               std::pair<std::pair<double, double>, std::pair<double, double>> domain) 
        {
            this->GRID_SIZE = grid_size;
            this->vars      = variables;
            this->domain    = domain;   
            this->h         = 1.0 / (grid_size - 2 + 1.0);
        }


    void initialise_grid ()
    {
        int N = GRID_SIZE - 2; // inner grid dimension

        std::cout << "Initialising grid:\n";
        std::cout << "h = " << h << '\n';
        std::cout << "GRID SIZE = " << GRID_SIZE << '\n';

        this->grid.resize(GRID_SIZE, GRID_SIZE);
        double boundary_val;

        for (int i = 0; i < GRID_SIZE; i++)
        {
            for (int j = 0; j < GRID_SIZE; j++)
            {
                this->grid(i, j) = 0.0;

                // bottom boundary
                if (i == 0)
                {
                    this->grid(i, j) = bc.bottom.evaluate(std::make_pair(i, j), vars, h);
                }

                // top boundary
                if (i == GRID_SIZE - 1)
                {
                    this->grid(i, j) = bc.top.evaluate(std::make_pair(i, j), vars, h);
                }

                // left boundary
                if (j == 0 && (i > 0 && i < GRID_SIZE - 1)) // keep top value
                {
                    this->grid(i, j) = bc.left.evaluate(std::make_pair(i, j), vars, h);
                }

                // right boundary
                if (j == GRID_SIZE - 1 && (i > 0 && i < GRID_SIZE - 1))
                {
                    this->grid(i, j) = bc.right.evaluate(std::make_pair(i, j), vars, h);
                }
            }
        }
    }

    LinearSystem<Matrix, Vector> construct_system ()
    {   
        std::cout << "Constructing system...\n";
        LinearSystem<Matrix, Vector> system;
        int N    = GRID_SIZE - 2; // inner grid dimension
        int dim  = int(std::pow(N, 2));

        system.A = Eigen::MatrixXd::Zero(dim, dim);
        system.b = Eigen::VectorXd::Zero(dim);
        system.u = Eigen::VectorXd::Zero(dim);
        system.N = GRID_SIZE;

        /* ( 4 u_{i, j} - u_{i+1, j} - u_{i-1, j} - u_{i, j+1} - u_{i, j-1} ) / h^2 = 0 */
        int k;
        for (int i = 1; i <= N; i++)
        {
            for (int j = 1; j <= N; j++)
            {
                k = (i - 1) * N + (j - 1);
                // printf("k: %d i: %d j: %d\n", k, i, j);

                system.A(k, k) = 4.0;

                if (i == N) // top neighbor
                    system.b(k) += grid(i + 1, j);
                else
                    system.A(k, k + N) = -1.0;

                if (i == 1) // bottom neighbor
                    system.b(k) += grid(i - 1, j);
                else
                    system.A(k, k - N) = -1.0;


                if (j == N) // right neighbor
                    system.b(k) += grid(i, j + 1) ;
                else
                    system.A(k, k + 1) = -1.0;


                if (j == 1) // left neighbor
                    system.b(k) += grid(i, j - 1);
                else
                    system.A(k, k - 1) = -1.0;
            }
        }
        return system;
    }

    void fill_grid (Vector u)
    {
        int N = this->grid.rows();
        int K = u.size();
        int i, j;

        for (int k = 0; k < K; k++)
        {
            i = k / (N-2) + 1;
            j = k % (N-2) + 1;
            this->grid(i, j) = u(k);
        }
    }

    void save_grid (std::string folder, std::string filename)
    {   
        int N = this->grid.rows();
        std::filesystem::path full_path = std::filesystem::current_path() / folder / filename;
        std::filesystem::create_directories(full_path.parent_path());

        std::ofstream file(full_path);

        file << "C h:"         << this->h             << '\n';
        file << "C GRID SIZE:" << this->GRID_SIZE     << '\n';
        file << "C N:"         << this->GRID_SIZE - 2 << '\n';

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                file << this->grid(i, j) << " "; 
            }
            file << '\n';
        }
        std::cout << "Grid saved to: " << full_path << '\n';
    }

    void print ()
    {
        
    }
};


#endif // LAPLACE_2D_HPP