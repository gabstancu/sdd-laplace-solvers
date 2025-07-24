#ifndef LAPLACE_2D_HPP
#define LAPLACE_2D_HPP

#include "basis/BoundaryCondition.hpp"
#include "basis/LinearSystem.hpp"

template<typename Matrix, typename Vector>
struct Laplace2D
{   
    int    GRID_SIZE;
    Matrix grid;
    Matrix analytical_solution;
    double h;
    
    std::pair<std::pair<double, double>, std::pair<double, double>> domain;

    BoundaryConditions         bc;
    std::vector<GiNaC::symbol> vars;
    GiNaC::ex                  analytical_expression;

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

                // top boundary
                if (i == 0)
                {
                    this->grid(i, j) = bc.top.evaluate(std::make_pair(j, GRID_SIZE - 1 - i), vars, h);
                }

                // bottom boundary
                if (i == GRID_SIZE - 1)
                {
                    this->grid(i, j) = bc.bottom.evaluate(std::make_pair(j, GRID_SIZE - 1 - i), vars, h);
                }

                // left boundary
                if (j == 0 && (i > 0 && i < GRID_SIZE - 1)) // keep top value
                {
                    this->grid(i, j) = bc.left.evaluate(std::make_pair(j, GRID_SIZE - 1 - i), vars, h);
                }

                // right boundary
                if (j == GRID_SIZE - 1 && (i > 0 && i < GRID_SIZE - 1))
                {
                    this->grid(i, j) = bc.right.evaluate(std::make_pair(j, GRID_SIZE - 1 - i), vars, h);
                }
            }
        }
    }


    LinearSystem<Matrix, Vector> construct_system ()
    {   
        // std::cout << "Constructing system...\n";
        LinearSystem<Matrix, Vector> system;
        int N    = GRID_SIZE - 2; // inner grid dimension
        int dim  = int(std::pow(N, 2));

        system.A = Eigen::MatrixXd::Zero(dim, dim);
        system.b = Eigen::VectorXd::Zero(dim);
        system.u = Eigen::VectorXd::Zero(dim);
        system.N = N;

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
                // if (i == 1) 
                //     system.b(k) += grid(i + 1, j);
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


    void save_grid (std::string filename, bool analytical, std::string folder = "")
    {    
        std::string folder_name = folder;
        if (folder.empty())
        {
            std::string folder_name = "grids/";
        }
        int N = this->grid.rows();

        /*  log approximation  */
        std::filesystem::path full_path = std::filesystem::current_path() / folder_name / filename;
        std::filesystem::create_directories(full_path.parent_path());
        // std::cout << "FULL PATH: " << full_path << "\n\n";

        std::ofstream file(full_path);

        file << "C h:"         << this->h             << '\n';
        file << "C GRID SIZE:" << this->GRID_SIZE     << '\n';
        file << "C N:"         << this->GRID_SIZE - 2 << '\n';

        if (!analytical)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    file << this->grid(i, j) << " "; 
                }
                file << '\n';
            }
            std::cout << "Approximation saved to: " << full_path << '\n';
            return;
        }
        else if (analytical)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    file << this->analytical_solution(i, j) << " "; 
                }
                file << '\n';
            }
            std::cout << "Analytical solution saved to: " << full_path << '\n';
            return;
        }
    }


    void evaluate_analytical_solution ()
    {   
        this->analytical_solution.resize(GRID_SIZE, GRID_SIZE);

        double x_ = 0.0, y_ = 0.0;

        for (int i = 0; i < GRID_SIZE; i++)
        {   
            y_ = h * (GRID_SIZE - 1 - i); 
            for (int j = 0; j < GRID_SIZE; j++)
            {
                x_ = h * j;
                this->analytical_solution(i, j) = this->grid(i, j);
                // printf("(%d, %d), (x_, y_) -> (%f, %f)\n", i, j, x_, y_);
                if ((i>=1 && i <= GRID_SIZE - 2) && (j>=1 && j <= GRID_SIZE - 2))
                {
                    GiNaC::exmap m;
                    m[vars[0]] = x_; m[vars[1]] = y_;
                    GiNaC::ex evaluated_expression  = this->analytical_expression.subs(m).evalf();
                    double    value                 = GiNaC::ex_to<GiNaC::numeric>(evaluated_expression).to_double();
                    this->analytical_solution(i, j) = value;
                }
            }
        }
    }


    void print ()
    {
        
    }
};


#endif // LAPLACE_2D_HPP