include("../../include/funcs_sym.jl")
include("../../include/func_sym_dt.jl")
include("../../include/func_DzuikSurface.jl")

function u₁(x,y,z)
    return x*y
end

function u₂(t)
    return sin(t + π/2)
end

cd(@__DIR__)
negΔₛu1 = negΔₛ(u₁, Dziuk_surface, "cfunc_negLBO_u1.jl")
do_t_u2 = dₜu(u₂, "cfunc_dt_u2.jl")
