#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP

#include <ginac/ginac.h>
#include <optional>
#include <variant>
#include <iostream>

enum class Side { Top, Bottom, Left, Right };
enum class Type { Dirichlet, Neumann, Robin };

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
};


struct BoundaryConditions
{
    BoundaryCondition top{std::nullopt, Side::Top}, 
                      bottom{std::nullopt, Side::Bottom}, 
                      left{std::nullopt, Side::Left}, 
                      right{std::nullopt, Side::Right} ;

};


#endif // BOUNDARY_CONDITIONS_HPP