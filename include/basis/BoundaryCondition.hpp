#ifndef BOUNDARY_CONDITION_HPP
#define BOUNDARY_CONDITION_HPP

#include <ginac/ginac.h>
#include <optional>
#include <variant>
#include <iostream>

enum class Side { Top, Bottom, Left, Right };

inline std::string to_string(Side side) {
    switch (side) {
        case Side::Top:    return "Top";
        case Side::Bottom: return "Bottom";
        case Side::Left:   return "Left";
        case Side::Right:  return "Right";
        default:           return "Unknown side.";
    }
}

enum class Type { Dirichlet, Neumann, Robin };

inline std::string to_string(Type type) {
    switch (type) {
        case Type::Dirichlet: return "Dirichlet";
        case Type::Neumann:   return "Neumann";
        case Type::Robin:     return "Robin";
        default:              return "Unknown type of boundary condition.";
    }
}


using BoundaryType = std::variant<GiNaC::ex, double>;

struct BoundaryCondition
{
    std::optional<BoundaryType> expr = std::nullopt;
    Side side;
    Type type;

    bool is_set()
    {
        return expr.has_value();
    }

    bool is_expression()
    {
        return expr.has_value() && std::holds_alternative<GiNaC::ex>(expr.value());
    }

    bool is_constant()
    {
        return expr.has_value() && std::holds_alternative<double>(expr.value());   
    }

    double evaluate (std::pair<int, int> coordinates, 
                     std::vector<GiNaC::symbol> vars, 
                     double h)
    {
        double val = 0.0;
        int i = coordinates.first, j = coordinates.second;

        if (this->is_set())
        {
            if (this->is_constant())
            {   
                val = std::get<double>(*expr);
                return val;
            }
            else if (this->is_expression())
            {
                GiNaC::exmap m; 
                m[vars[0]] = j * h; m[vars[1]] = i * h;
                GiNaC::ex b_expr = std::get<GiNaC::ex>(*expr);
                GiNaC::ex evaluated_boundary = b_expr.subs(m).evalf();
                val = GiNaC::ex_to<GiNaC::numeric>(evaluated_boundary).to_double();
                return val;
            }
        }
        return val;
    }

    void print ()
    {
        if (!is_set())
        {
            std::cout << "Side " << to_string(side) << " was not initialised\n.";
            return;
        }
        else
        {
            std::cout << "Side: " << to_string(side) << '\n';
            if (is_constant())
            {
                std::cout << "Value: " << std::get<double>(*expr) << '\n';
            }
            else if (is_expression())
            {
                std::cout << "Expression: " << std::get<GiNaC::ex>(*expr) << '\n';
            }
        }
    }
};

struct BoundaryConditions
{
    BoundaryCondition top    {std::nullopt, Side::Top,    Type::Dirichlet};
    BoundaryCondition bottom {std::nullopt, Side::Bottom, Type::Dirichlet};
    BoundaryCondition left   {std::nullopt, Side::Left,   Type::Dirichlet}; 
    BoundaryCondition right  {std::nullopt, Side::Right,  Type::Dirichlet};
};

#endif // BOUNDARY_CONDITION_HPP