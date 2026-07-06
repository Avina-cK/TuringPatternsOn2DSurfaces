cd(@__DIR__)
using RuntimeGeneratedFunctions
using NaNMath
include("../../include/func_DzuikSurface.jl")

# compiled function for the r.h.s of the PDE
rhs_func = include("rhs_f_0.jl")

function u_chosen(x,y,z)
    return x
end

@info "Loaded u_chosen(x,y,z)=x and rhs_func=-Δₛu_chosen"