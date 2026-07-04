using LinearAlgebra
using Symbolics

# Define the variables (do not use commas!)
@variables sym_x₁ sym_x₂ sym_x₃

"""
normal_atpt(d::Function)
    - d(x,y,z) is the level set function that defines the surface (d=0 on the surface)
    The function returns a callable that computes the normal at any point (x,y,z).
        n = ∇d / |∇d|
Example:
    d(x,y,z) = z^2 + y^2 + (x - z^2)^2 - 1.0
    n_d = normal_atpt(d)
    n_d(0.0, -0.50, -1.0)
Then you can evaluate n_d at any point (x,y,z) to get the value of the normal at that point.
"""
function normal_atpt(d::Function)
    d_sym = d(sym_x₁, sym_x₂, sym_x₃)
    gradd = Symbolics.gradient(d_sym, [sym_x₁, sym_x₂, sym_x₃])
    norm_n = sqrt(sum(gradd.^2))
    normalfuncsym = gradd ./ norm_n
    normalfunc_out,_ = build_function(normalfuncsym, sym_x₁, sym_x₂, sym_x₃;
                                 expression = Val{false}
                                 , simplify=Val{true}
                                 , target=Symbolics.JuliaTarget())
    return (x,y,z) -> normalfunc_out(x,y,z)
end

"""
∇ₛ(u::Function, d::Function)
    - u(x,y,z) is the function for which we want to compute the tangential gradient
    - d(x,y,z) is the level set function that defines the surface (d=0 on the surface)
    The function returns a callable that computes the tangential gradient at any point (x,y,z).
The tangential gradient, on a surface,d=0 , is given by:
    ∇ₛ(u) = ∇u - (∇u ⋅ n) n, where n = ∇d / |∇d| is the unit normal vector to the surface.
Example:
    d(x,y,z) = z^2 + y^2 + (x - z^2)^2 - 1.0
    u(x,y,z) = x*y
    ∇ₛ_u_given = ∇ₛ(u,d)
    ∇ₛ_u_given(0.0, -0.50, -1.0)
Then you can evaluate ∇ₛ_u_given at any point (x,y,z) to get the value of the tangential gradient at that point.
"""
function ∇ₛ(u::Function, d::Function; outputfile::Bool=false, opfilename::String="func_tangrad.jl")
    u_sym = u(sym_x₁, sym_x₂, sym_x₃);
    d_sym = d(sym_x₁, sym_x₂, sym_x₃);
    gradu = Symbolics.gradient(u_sym, [sym_x₁, sym_x₂, sym_x₃]);
    gradd = Symbolics.gradient(d_sym, [sym_x₁, sym_x₂, sym_x₃]);
    norm_n = sqrt(sum(gradd.^2));
    n = gradd ./ norm_n;
    TG_sym = gradu .- (LinearAlgebra.dot(gradu, n) * n);
    TG_func, _ = build_function(TG_sym, sym_x₁, sym_x₂, sym_x₃;
                                 expression = Val{true}
                                 , simplify=Val{true}
                                 , target=Symbolics.JuliaTarget())
    if outputfile
        write(opfilename, string(TG_func))
    end
    return (x,y,z) -> TG_func(x,y,z)
end

"""
Δₛ(u::Function, d::Function)
    - u(x,y,z) is the function for which we want to compute the tangential gradient
    - d(x,y,z) is the level set function that defines the surface (d=0 on the surface)
    The function returns a callable that computes the Laplace-Beltrami operator on u at any point (x,y,z).
The Laplace-Beltrami operator on u (on a surface, d=0) is given by:
    Δₛu = ∇ₛ ⋅ ∇ₛ(u)
Example:
    d(x,y,z) = z^2 + y^2 + (x - z^2)^2 - 1.0
    u(x,y,z) = x*y
    Δₛ_u_given = Δₛ(u, d)
    Δₛ_u_given(0.0, -0.50, -1.0)
Then you can evaluate ∇ₛ_u_given at any point (x,y,z) to get the value of the tangential gradient at that point.
"""
function Δₛ(u::Function, d::Function; outputfile::Bool=false, opfilename::String="func_tangrad.jl")
    u_sym = u(sym_x₁, sym_x₂, sym_x₃)
    d_sym = d(sym_x₁, sym_x₂, sym_x₃)
    gradu = Symbolics.gradient(u_sym, [sym_x₁, sym_x₂, sym_x₃])
    gradd = Symbolics.gradient(d_sym, [sym_x₁, sym_x₂, sym_x₃])
    norm_n = sqrt(sum(gradd.^2))
    n = gradd ./ norm_n
    TG_sym = gradu .- (sum(gradu .* n) .* n)

    # Compute jacobian of tangential gradient
    JTG_sym = Symbolics.jacobian(TG_sym, [sym_x₁, sym_x₂, sym_x₃])
    Δₛ_sym = JTG_sym[1,1] + JTG_sym[2,2] + JTG_sym[3,3]  # trace = divergence

    Δₛ_func,_ = build_function(Δₛ_sym, sym_x₁, sym_x₂, sym_x₃;
                                 expression = Val{false}
                                 , simplify=Val{true}
                                 , target=Symbolics.JuliaTarget())
    if outputfile
        write(opfilename, string(Δₛ_func))
    end
   return (x,y,z) -> Δₛ_func(x,y,z)
end

"""
negΔₛ(u::Function, d::Function)
    - u(x,y,z) is the function for which we want to compute the tangential gradient
    - d(x,y,z) is the level set function that defines the surface (d=0 on the surface)
    The function returns a callable that computes the Laplace-Beltrami operator on u at any point (x,y,z).
The Laplace-Beltrami operator on u (on a surface, d=0) is given by:
    -Δₛu = -∇ₛ ⋅ ∇ₛ(u)
"""
function negΔₛ(u::Function, d::Function, opfilename::String="cfunc_negLBO.jl")
    u_sym = u(sym_x₁, sym_x₂, sym_x₃)
    d_sym = d(sym_x₁, sym_x₂, sym_x₃)
    gradu = Symbolics.gradient(u_sym, [sym_x₁, sym_x₂, sym_x₃])
    gradd = Symbolics.gradient(d_sym, [sym_x₁, sym_x₂, sym_x₃])
    norm_n = sqrt(sum(gradd.^2))
    n = gradd ./ norm_n
    TG_sym = gradu .- (sum(gradu .* n) .* n)

    # Compute jacobian of tangential gradient
    JTG_sym = Symbolics.jacobian(TG_sym, [sym_x₁, sym_x₂, sym_x₃])
    nΔₛ_sym = -(JTG_sym[1,1] + JTG_sym[2,2] + JTG_sym[3,3])  # trace = divergence

    nΔₛ_expr = build_function(nΔₛ_sym, sym_x₁, sym_x₂, sym_x₃;
                                 expression = Val{true},
                                 simplify = Val{true})
    
    code = """
        $(nΔₛ_expr)
    """
    
    write(opfilename, code)

    # Return function for immediate use
    return build_function(nΔₛ_sym, sym_x₁, sym_x₂, sym_x₃;
                         expression = Val{false},
                         simplify = Val{true})
end
@info "Loaded symbolically defined functions: normal_atpt(d), ∇ₛ(u,d) and Δₛ(u,d)."

#= De-comment the following if you want to (re-)generate functions for the Dzuik surface and given u
include("func_DzuikSurface.jl")

function u_chosen(x,y,z)
    return x*y
end
cd(@__DIR__)
#∇ₛu_chosen = ∇ₛ(u_chosen, Dziuk_surface; outputfile=true, opfilename="cfunc_tangrad.jl")
#Δₛu_chosen = Δₛ(u_chosen, Dziuk_surface; outputfile=true, opfilename="cfunc_LBO.jl")
negΔₛu_chosen = negΔₛ(u_chosen, Dziuk_surface, "cfunc_negLBO.jl")
=#