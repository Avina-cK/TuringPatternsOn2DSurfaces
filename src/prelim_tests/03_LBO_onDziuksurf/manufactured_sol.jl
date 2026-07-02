cd(@__DIR__)
using RuntimeGeneratedFunctions
using NaNMath
include("../../include/func_DzuikSurface.jl")

# compiled function for the r.h.s of the PDE
rhs_func = include("../../include/cfunc_negLBO.jl")

function u_chosen(x,y,z)
    return x*y
end

@info "Loaded u_chosen(x,y,z)=xy and rhs_func=-Δₛu_chosen"