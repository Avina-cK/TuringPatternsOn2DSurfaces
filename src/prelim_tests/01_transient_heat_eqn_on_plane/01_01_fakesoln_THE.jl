#=
This file generates a manufactured solution for the transient heat equation on a plane (Ω=[-L₁,L₁] x [-L₂,L₂]), i.e., 
    ∂u/∂t - ∇ · (k ∇u) = f on Ω × (0, Tf]
    u = 0 on ∂Ω × (0, Tf]
where u is the temperature, k is the thermal conductivity, and f is a source term.
    The solution is constructed using a bump function (steady state solution) and the eigenfunctions of the Laplacian with zero Dirichlet boundary conditions. 
    The source term f is derived from the manufactured solution by substituting it into the PDE and computing the necessary derivatives, via the Symbolics.jl package.
=#
using Symbolics

@variables x₁ x₂ tᵥ

const Dx₁ = Differential(x₁)
const Dx₂ = Differential(x₂)
const Dt  = Differential(tᵥ)

"""
bump function on a square domain
    b_f(a, b) = exp(-1.0 / ((1.0 - (a/L₁)^2) * (1.0 - (b/L₂)^2)))
    - a, b are the spatial coordinates in x₁ and x₂ directions respectively
    - L₁, L₂ are the half-lengths of the domain in x₁ and x₂ directions respectively
"""
b_f(a, b) = exp(-1.0 / ((1.0 - (a/L₁)^2) * (1.0 - (b/L₂)^2)))

"""
Eigenvalues of the Laplacian on a square domain with zero Dirichlet BCs
    λ(n::Int, m::Int) = (π^2/4) * (n^2/L₁^2 + m^2/L₂^2)
    - n, m are the mode numbers in x₁ and x₂ directions respectively
    - L₁, L₂ are the half-lengths of the domain in x₁ and x₂ directions respectively
"""
λ(n::Int, m::Int) = (π^2/4) * (n^2/L₁^2 + m^2/L₂^2)

"""
Eigenvectors of the Laplacian on a square domain with zero Dirichlet BCs
    v_nm(n::Int, m::Int, a, b)  = sin(n*π*(a+L₁)/(2L₁)) * sin(m*π*(b+L₂)/(2L₂))
    - n, m are the mode numbers in x₁ and x₂ directions respectively
    - a, b are the spatial coordinates in x₁ and x₂ directions respectively
    - L₁, L₂ are the half-lengths of the domain in x₁ and x₂ directions respectively
"""
v_nm(n::Int, m::Int, a, b)  = sin(n*π*(a+L₁)/(2L₁)) * sin(m*π*(b+L₂)/(2L₂))

ne, me = 1, 1

g₁(tv) = -λ(ne,me) * k * exp(-λ(ne,me) * k * tv)
g₂(a, b) =  λ(ne,me) * v_nm(ne, me, a, b)
u_true(a, b, tv) = g₁(tv) * g₂(a, b) + b_f(a, b)

# manufactured solution in symbolic form for Symbolics.jl to compute derivatives on
u_sym   = u_true(x₁, x₂, tᵥ)
# source term in symbolic form, derived by substituting the manufactured solution into the PDE 
rhs_sym = expand_derivatives(Dt(u_sym) - k*(Dx₁(Dx₁(u_sym)) + Dx₂(Dx₂(u_sym)))) 

# build callable functions from symbolic expressions
u_true_fn = build_function(u_sym,   x₁, x₂, tᵥ, expression=Val{false})
rhs_f_fn = build_function(rhs_sym, x₁, x₂, tᵥ, expression=Val{false})

@info "Generated manufactured solution, u_true_fn, and source term, rhs_f_fn, for the transient heat equation"


#= Visualise the functions
using Plots
L₁, L₂ = 1.0, 1.0
b_f(a, b) = exp(-1.0 / ((1.0 - (a/L₁)^2) * (1.0 - (b/L₂)^2)))
λ(n::Int, m::Int) = (π^2/4) * (n^2/L₁^2 + m^2/L₂^2)
v_nm(n::Int, m::Int, a, b)  = sin(n*π*(a+L₁)/(2L₁)) * sin(m*π*(b+L₂)/(2L₂))
ne, me = 1, 1
k = 1.0e-3
g₁(tv) = -λ(ne,me) * k * exp(-λ(ne,me) * k * tv)
g₂(a, b) =  λ(ne,me) * v_nm(ne, me, a, b)
u_true(a, b, tv) = g₁(tv) * g₂(a, b) + b_f(a, b)
gm7_r = cgrad([colorant"rgb(61, 26, 120)"
                    , colorant"rgb(80, 73, 204)"
                    , colorant"rgb(123, 173, 226)"
                    , colorant"rgb(255,255,255)"
                    , colorant"rgb(152, 232, 193)"
                    , colorant"rgb(38, 206, 170)"
                    ,colorant"rgb(7, 141, 112)"]
                );
les7 = cgrad([colorant"rgb(213, 45, 0)"
                    , colorant"rgb(239, 118, 39)"
                    , colorant"rgb(255, 154, 86)"
                    , colorant"rgb(255,255,255)"
                    , colorant"rgb(209, 98, 164)"
                    , colorant"rgb(181, 86, 144)"
                    , colorant"rgb(163, 2, 98)"]
                );
x1, x2 = -L₁:0.01:L₁, -L₂:0.01:L₂
p1 = plot(x1, x2, b_f, st=:surface
    , title="b_f(x₁,x₂)=exp(-1/((1-x₁²)(1-x₂²)))", titlefontsize=10
    , xlabel="x₁", ylabel="x₂", zlabel="b_f(x₁,x₂)"
    , xlims=(-L₁, L₁), ylims=(-L₂, L₂)
    , aspect_ratio=:equal
    , color=gm7_r
    )
p2 = plot(x1, x2, (a,b) -> g₂(a, b), st=:surface
    , title="g₂(x₁,x₂)=λ₁₁⋅v₁₁(x₁,x₂)", titlefontsize=10
    , xlabel="x₁", ylabel="x₂", zlabel="g₂(x₁,x₂)"
    , xlims=(-L₁, L₁), ylims=(-L₂, L₂)
    , aspect_ratio=:equal
    , color=gm7_r
    )
p3 = plot(0.0:10:1000.0, (t) -> g₁(t), st=:line
    , title="g₁(t)=-λ₁₁⋅k⋅exp(-λ₁₁⋅k⋅t)"
    , xlabel="t", ylabel="g₁(t)"
    , xlims=(0.0, 1000.0)
    , color=:blue, legend=false
    )
p4 = plot(x1, x2, (a,b) -> u_true(a, b, 1.5), st=:surface
    , title="u_true = g₁(t) ⋅ g₂(x₁,x₂) + b_f(x₁,x₂)"
    , xlabel="x₁", ylabel="x₂", zlabel="u_true(x₁,x₂,1.5)"
    , xlims=(-L₁, L₁), ylims=(-L₂, L₂)
    , color=les7
    )
plot(p1, p2, p3,p4, layout=4, tickfontsize=6, guidefont=(8), titlefont=(9))
=#