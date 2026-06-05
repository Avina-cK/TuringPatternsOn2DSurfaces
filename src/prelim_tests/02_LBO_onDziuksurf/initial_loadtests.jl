using Ferrite, Tensors, SparseArrays
include("setup_surface.jl")

## Load mesh
include("../../include/func_loadmesh.jl")
refinement = 5
filename = "lowresmesh_$(refinement).msh"
cd(@__DIR__)
cd("../../Dziuk_surf_meshes/")
Ωₕ, h = load_mesh(filename);
# flip to counter-clockwise orientation
flipped_cells = [Triangle((c.nodes[2], c.nodes[1], c.nodes[3])) for c in Ωₕ.cells]
Ωₕ = Grid(flipped_cells, Ωₕ.nodes)

@info "Grid imported for refinement level $refinement"
cd(@__DIR__)
#= Ωₕ.nodes; Ωₕ.cells =#

##
qr = QuadratureRule{RefTriangle}(1)
ip = Lagrange{RefTriangle, 1}()
cellvalues_Ω = setup_surface_cellvalues(Ωₕ)

dh = DofHandler(Ωₕ)
add!(dh, :u, ip)
close!(dh)

## Check dofs
# for P1 scalar field: 1 DoF per node
@assert ndofs(dh) == getnnodes(Ωₕ) "Error: no. DoFs ≠ no. nodes"

# Verify all cells have 3 dofs (P1 triangle)
for i in 1:getncells(Ωₕ)
    @assert length(celldofs(dh, i)) == 3 "Error: Cell $i does not have 3 dofs"
end

##
n_dofs = ndofs(dh)
include("func_assembleKF.jl")
F = zeros(n_dofs)
include("manufactured_sol.jl")

Kₕ, F = assemble_globalKF(F, dh, cellvalues_Ω, rhs_func)

## Testing Kₕ
@assert Kₕ == Kₕ' "Kₕ is not symmetric"   #symmetry
@assert norm(Kₕ * ones(n_dofs))<1.0e-12 "1-vec is not in null-space of Kₕ"

## TODO

include("func_assembleconstraintC.jl")
C = zeros(n_dofs)
assemble_constraint_C!(C, dh, qr, ip)
for i in 1:n_dofs; @assert C[i]>0 "C has zero-valued elements"; end 

include("func_buildlagrangesys.jl")
A, b = build_lagrange_system(Kₕ, F, C)
uμ = A\b
u_h = uμ[1:n_dofs]

include("func_assemble_u_true.jl")
u_true_vec = assemble_true_solution(u_chosen, dh,ip)

mean_u_true = dot(C, u_true_vec) / sum(C)
u_true_centered = u_true_vec .- mean_u_true

@assert abs(dot(C, u_true_centered)) < 1e-10 "u_true not zero-mean after centering"

include("func_L2_norm.jl")
l2_err  = L2_error(u_true_centered, u_h, dh, qr, ip)
