function create_embedded_sphere(refinements)
    gmsh.initialize()

    # Add a unit sphere in 3D space
    gmsh.model.occ.addSphere(0.0, 0.0, 0.0, 1.0)
    gmsh.model.occ.synchronize()

    # Generate nodes and surface elements only, hence we need to pass 2 into generate
    gmsh.model.mesh.generate(2)

    # To get good solution quality refine the elements several times
    for _ in 1:refinements
        gmsh.model.mesh.refine()
    end

    # Now we create a Ferrite grid out of it. Note that we also call toelements
    # with our surface element dimension to obtain these.
    nodes = tonodes()
    elements, _ = toelements(2)
    gmsh.finalize()
    return Grid(elements, nodes)
end

rhs_f(x,y,z) = y

Ωₕ = create_embedded_sphere(2)

@info "Grid imported for refinement level $refinement"
cd(@__DIR__)
#= Ωₕ.nodes; Ωₕ.cells =#

## Setup FE spaces and cell values
qr = QuadratureRule{RefTriangle}(1)
ip = Lagrange{RefTriangle, 1}()
cellvalues_Ω = setup_surface_cellvalues(Ωₕ)

dh = Ferrite.DofHandler(Ωₕ)
Ferrite.add!(dh, :u, ip)
Ferrite.close!(dh)

## Check dofs 
# for P1 scalar field: 1 DoF per node
@assert Ferrite.ndofs(dh) == Ferrite.getnnodes(Ωₕ) "Error: no. DoFs ≠ no. nodes"

# Verify all cells have 3 dofs (P1 triangle)
for i in 1:getncells(Ωₕ)
    @assert length(celldofs(dh, i)) == 3 "Error: Cell $i does not have 3 dofs"
end

## Assemble global stiffness matrix Kₕ and load vector F
n_dofs = ndofs(dh)
include("func_assembleKF_manual.jl")
F = zeros(n_dofs)
include("manufactured_sol.jl")

Kₕ, F = assemble_globalKF(dh, cellvalues_Ω, rhs_f)

## Testing Kₕ
@assert Kₕ == Kₕ' "Kₕ is not symmetric"   #symmetry
@assert norm(Kₕ * ones(n_dofs))<1.0e-12 "1-vec is not in null-space of Kₕ"

## Assemble constraint vector C

include("func_assembleconstraintC.jl")
C = zeros(n_dofs)
C = assemble_globalC(dh, cellvalues_Ω)
for i in 1:n_dofs; @assert C[i]>0 "C has zero-valued elements"; end 

## Build Lagrange system and solve for u_h
include("func_buildlagrangesys.jl")
A, b = build_lagrange_system(Kₕ, F, C)
uμ = A\b
u_h = uμ[1:n_dofs]

## Compute L2 error against true solution
include("func_assemble_u_true.jl")

u_true_vec = assemble_true_solution(rhs_f, dh,ip)

mean_u_true = dot(C, u_true_vec) / sum(C)
u_true_centered = u_true_vec .- mean_u_true

@assert abs(dot(C, u_true_centered)) < 1e-10 "u_true not zero-mean after centering"

include("func_L2_norm.jl")
l2_err  = L2_error(u_true_centered, u_h, dh, qr, ip)
