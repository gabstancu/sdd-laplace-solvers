#ifndef LAPLACE_2D_HPP
#define LAPLACE_2D_HPP

#include "basis/BoundaryCondition.hpp"
#include "basis/LinearSystem.hpp"

template<typename Matrix, typename Vector>
struct Laplace2D
{   
    int GRID_SIZE;
    Matrix grid;
    double h;
    std::pair<std::pair<double, double>, std::pair<double, double>> domain;
    BoundaryConditions bc;
    std::vector<GiNaC::symbol> vars;

    Laplace2D (int grid_size, std::vector<GiNaC::symbol> variables) :
        GRID_SIZE(grid_size), vars(variables), h(1.0 / (grid_size - 1.0)) {}


    void initialise_grid ()
    {
        int N = GRID_SIZE - 2; // inner grid dimension

        std::cout << "Initialising grid:\n";
        std::cout << "h = " << h << '\n';
        std::cout << "GRID SIZE = " << GRID_SIZE << '\n';

        this->grid.resize(GRID_SIZE, GRID_SIZE);
        double boundary_val;

        std::optional<BoundaryType> boundary;

        for (int i = 0; i < GRID_SIZE; i++)
        {
            for (int j = 0; j < GRID_SIZE; j++)
            {
                this->grid(i, j) = 0.0;

                // top boundary
                if (i == 0)
                {
                    this->grid(i, j) = bc.top.evaluate(std::make_pair(i, j), vars, h);
                }

                // bottom boundary
                if (i == GRID_SIZE - 1)
                {
                    this->grid(i, j) = bc.bottom.evaluate(std::make_pair(i, j), vars, h);
                }

                // left boundary
                if (j == 0 && (i > 0 && i < GRID_SIZE - 1)) // keep top's value
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
        LinearSystem<Matrix, Vector> system;
        int N = GRID_SIZE - 2; // inner grid dimension
        int dim = int(std::pow(N, 2));
    }

    void fill_grid ()
    {

    }

    void save_grid (std::string label)
    {

    }
};


#endif // LAPLACE_2D_HPP