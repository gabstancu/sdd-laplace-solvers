#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP

#include <ginac/ginac.h>
#include <optional>
#include <variant>

enum class Side { Top, Bottom, Left, Right };

using BoundaryType = std::variant<GiNaC::ex, double>;

struct BoundaryCondition
{
    std::optional<BoundaryType> expr = std::nullopt;
    Side side;

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
};


struct BoundaryConditions
{
    BoundaryCondition top{std::nullopt, Side::Top}, 
                      bottom{std::nullopt, Side::Bottom}, 
                      left{std::nullopt, Side::Left}, 
                      right{std::nullopt, Side::Right} ;
};


#endif // BOUNDARY_CONDITIONS_HPP