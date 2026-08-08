#=
PDE:
    ∂ₜu - Dᵤ⋅Δₛu = α(1-u) - uv²
    ∂ₜv - Dᵥ⋅Δₛv = -(α+β)v + uv²
=#

using Ferrite, FerriteGmsh
using BlockArrays, LDLFactorizations, SparseArrays, LinearAlgebra
using Random
using WriteVTK

include("../../include/surface_cell_values.jl")
include("../../include/gen_read_meshes/func_load_tri_surfmesh.jl")

struct GrayScottMaterial{T}
    Dᵤ::T
    Dᵥ::T
    α::T
    β::T
end

## Parameters
given_m = GrayScottMaterial{Float64}(0.00016, 0.00008, 0.06, 0.062)
T = 3200.0
Δt = 10.0
nsteps = ceil(Int, T/Δt)
save_every = Int(floor(nsteps/50))
α, β = given_m.α, given_m.β

## Initialise mesh
refinement = 2
filename = "Dzuik_mesh_$(refinement).msh"
cd(@__DIR__)
cd("../../Dziuk_surf_meshes/")
Ωₕ, h = load_tri_surfmesh(filename)

## Cell values
ip = Lagrange{RefTriangle, 1}()
qr = QuadratureRule{RefTriangle}(1)
cellvalues_Ω = SurfaceCellValues(qr, ip)

# Two scalar fields
dh = DofHandler(Ωₕ)
add!(dh, :u, ip)
add!(dh, :v, ip)
close!(dh)

ndof   = ndofs(dh)

## initial condition
include("func_assemble_initial_cond.jl")
uv_0 = zeros(ndof)
uv_0 = setup_initial_conditions(uv_0, dh, 314)
## Assemble matrices
include("func_assemble_matrices.jl")

#Mₕ = allocate_matrix(dh; coupling = [true false; false true])
#Kₕ = allocate_matrix(dh; coupling = [true false; false true])
Mₕ = spzeros(ndof, ndof)
Kₕ = spzeros(ndof, ndof)

assemble_globalKM!(Mₕ, Kₕ, cellvalues_Ω, dh, given_m)
Aₕ = Mₕ + Δt .* Kₕ
#Aₕ = sparse(Aₕ)
## Zero-mean constraints, for u and v

eᵤ = zeros(ndof); eᵤ[u_dofs] .= 1.0
eᵥ = zeros(ndof); eᵥ[v_dofs] .= 1.0

Cᵤ = Mₕ * eᵤ    # nonzero only on u-dofs (M is block-diagonal in u,v)
Cᵥ = Mₕ * eᵥ    # nonzero only on v-dofs

C  = hcat(Cᵤ, Cᵥ)   # ndof × 2 constraint matrix

## Augmented system
Aug = [Aₕ sparse(C);
       sparse(C') spzeros(2, 2)]
Aₕ_F = ldl(Aug)

## Time-stepping: semi-implicit Euler
# w:=[u,v]
# (M + Δt K) wₙ₊₁ = M wₙ + Δt M R(wₙ)  with zero-mean

rhs_vecₕ = zeros(ndof+2) # +2 ∵ augmented system
Rₕ = zeros(ndof)
wₙ = copy(uv_0)

## Save solutions
cd(@__DIR__)
mkpath("results2")
cd("results2/")
pvd = paraview_collection("GS_Dzuiksurf_$(refinement)")
VTKGridFile("GS_Dzuiksurf_$(refinement)-0", dh) do vtk
    write_solution(vtk, dh, uv_0, "_h")
    pvd[0.0] = vtk
end

##
for (n, t) in enumerate(Δt:Δt:T)
    u = copy(wₙ[u_dofs])
    v = copy(wₙ[v_dofs])

    Rₕ[u_dofs] = (α .* (1.0 .- u)) .- (u .* v .^2)
    Rₕ[v_dofs] = -((α .+ β) .* v) .+ (u .* v .^2)

    rhs_vecₕ[1:ndof] .= (Mₕ * wₙ) .+ Δt .* (Mₕ * Rₕ)
    rhs_vecₕ[ndof+1:end] .= C' * uv_0 #zero-mean
    @assert !any(isnan, rhs_vecₕ) "NaN values in rhs"

    wₙ₊₁ = Aₕ_F \ rhs_vecₕ

    wₙ .= @view wₙ₊₁[1:ndof]
    @assert !any(isnan, wₙ) "NaN values in solution"
    
    if n % save_every == 0
        VTKGridFile("GS_Dzuiksurf_$(refinement)-$(n)", dh) do vtk
            write_solution(vtk, dh, wₙ, "_h")
            pvd[t] = vtk
        end
        @info "Step $(t)/$(T)"
    end

end
cd(@__DIR__)