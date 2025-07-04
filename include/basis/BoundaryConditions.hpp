#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP

#include <ginac/ginac.h>
#include <optional>
#include <variant>

enum class Side { Top, Bottom, Left, Right };

struct BoundaryCondition
{
    using BoundaryType = std::variant<GiNaC::ex, double>;

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
   
};


#endif // BOUNDARY_CONDITIONS_HPP