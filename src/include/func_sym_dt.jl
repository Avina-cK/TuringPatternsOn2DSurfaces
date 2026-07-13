using LinearAlgebra
using Symbolics

# Define the variables (do not use commas!)
@variables sym_t
Dt = Symbolics.Differential(sym_t)

function dₜu(u::Function, opfilename::String="cfunc_dt_u.jl")
    u_sym = u(sym_t)
    dt_u_sym = expand_derivatives(Dt(u_sym))
    
    dt_u_expr = build_function(dt_u_sym, sym_t;
                                 expression = Val{true})
    
    code = "dt_u_generated = $(dt_u_expr)\n"
    
    write(opfilename, code)

    # Return function for immediate use
    return build_function(dt_u_sym, sym_t;
                         expression = Val{false})
end