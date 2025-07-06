// #ifndef PDE_HPP
// #define PDE_HPP

// #include <vector>
// #include <map>
// #include <Eigen/Eigen>
// #include <ginac/ginac.h>

// enum class Side { Top, Bottom, Left, Right };

// struct BoundaryCondition
// {   
//     using BoundaryType = std::variant<GiNaC::ex, double>;

//     std::optional<BoundaryType> expr = std::nullopt;
//     Side side;

//     bool is_set()
//     {
//         return expr.has_value();
//     }

//     bool is_expression()
//     {
//         return expr.has_value() && std::holds_alternative<GiNaC::ex>(expr.value());
//     }

//     bool is_constant()
//     {
//         return expr.has_value() && std::holds_alternative<double>(expr.value());   
//     }
// };


// struct BoundaryConditions
// {   
//     using BoundaryType = std::variant<GiNaC::ex, double>;

//     std::optional<BoundaryType>  top   =  std::nullopt;
//     std::optional<BoundaryType> bottom =  std::nullopt;
//     std::optional<BoundaryType>  left  =  std::nullopt;
//     std::optional<BoundaryType> right  =  std::nullopt;


//     bool is_set(std::optional<BoundaryType> boundary)
//     {
//         return boundary.has_value();
//     }

//     bool is_expression(std::optional<BoundaryType> boundary)
//     {
//         return boundary.has_value() && std::holds_alternative<GiNaC::ex>(boundary.value());
//     }

//     bool is_constant(std::optional<BoundaryType> boundary)
//     {
//         return boundary.has_value() && std::holds_alternative<double>(boundary.value());   
//     }
// };


// Eigen::MatrixXd initialise_grid ( int GRID_SIZE, 
//                                   BoundaryConditions bc, 
//                                   std::vector<GiNaC::symbol> vars);


// void set_boundary_value ( Eigen::MatrixXd& grid, double h, 
//                           BoundaryConditions bc,
//                           std::pair<int, int> coordinates, 
//                           std::optional<BoundaryConditions::BoundaryType> boundary, 
//                           std::vector<GiNaC::symbol> vars);


// void construct_system ( Eigen::MatrixXd grid, int GRID_SIZE, 
//                         Eigen::MatrixXd& A, Eigen::VectorXd& b, Eigen::VectorXd& u);


// void fill_grid (Eigen::MatrixXd& grid, Eigen::VectorXd u);


// void save (Eigen::MatrixXd grid, std::string filename);

// #endif // PDE_HPP