cd(@__DIR__)
using RuntimeGeneratedFunctions
using NaNMath
include("../../include/func_DzuikSurface.jl")

function u₁(x,y,z)
    return x*y
end

function u₂(t)
    return sin(t + π/2)
end

u_chosen(x,y,z,t) = u₁(x,y,z) * u₂(t)

negΔₛu1 = include("cfunc_negLBO_u1.jl")

#dt_u2 = include("cfunc_dt_u2.jl")
dt_u2(t) = cos(t + π/2)

rhs_func(x,y,z,t,k) = ( u₁(x,y,z) * dt_u2(t) ) + (u₂(t) * k * negΔₛu1(x,y,z))

@info "Loaded u_chosen(x,y,z)=xy*sin(t + π/2) and rhs_func=-Δₛu_chosen"
