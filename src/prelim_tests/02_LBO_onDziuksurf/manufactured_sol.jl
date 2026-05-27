cd(@__DIR__)
using RuntimeGeneratedFunctions
using NaNMath
include("../../include/func_DzuikSurface.jl")
rhs_func = include("../../include/cfunc_negLBO.jl")
function u_chosen(x,y,z)
    return x*y
end
@info "Loaded u_chosen(x,y,z)=xy and rhs_func=-Δₛu_chosen"