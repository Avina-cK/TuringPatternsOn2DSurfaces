using Ferrite, Tensors, SparseArrays
include("surface_cell_values.jl")
include("../../include/read_Dzuik_mesh.jl")
include("../../include/funcs_error_analysis.jl")
include("func_assembleKF.jl")
include("manufactured_sol.jl")

## Load mesh

refinement = 2
filename = "lowresmesh_$(refinement).msh"
cd(@__DIR__)
cd("../../Dziuk_surf_meshes/")
Ωₕ, h = load_Dzuik_mesh(filename)

@info "Grid imported for refinement level $refinement"
cd(@__DIR__)
#= Ωₕ.nodes; Ωₕ.cells =#

## Setup FE spaces and cell values
qr = QuadratureRule{RefTriangle}(1)
ip = Lagrange{RefTriangle, 1}()

dh = DofHandler(Ωₕ)
add!(dh, :u, ip)
close!(dh)

cellvalues_Ω = SurfaceCellValues(qr, ip)

## Check dofs 
# for P1 scalar field: 1 DoF per node
@assert ndofs(dh) == getnnodes(Ωₕ) "Error: no. DoFs ≠ no. nodes"

# Verify all cells have 3 dofs (P1 triangle)
for i in 1:getncells(Ωₕ)
    @assert length(celldofs(dh, i)) == 3 "Error: Cell $i does not have 3 dofs"
end

## Assemble global stiffness matrix Kₕ and load vector F
n_dofs = ndofs(dh)
Kₕ, F = assemble_globalKF(dh, cellvalues_Ω, (x, y, z) -> rhs_func(x, y, z))

## Testing Kₕ
@assert Kₕ == Kₕ' "Kₕ is not symmetric"   #symmetry
@assert norm(Kₕ * ones(n_dofs))<1.0e-12 "1-vec is not in null-space of Kₕ"

## Assemble constraint vector C

_, c = assemble_globalKF(dh, cellvalues_Ω, (x, y, z) -> 1.0)  # c_i = ∫φᵢ dΩ
for i in 1:n_dofs; @assert c[i]>0 "c has zero-valued elements"; end 

## Build Lagrange system and solve for u_h
A = [Kₕ c; c' 0.0]
A = sparse(A)
b = [F; 0.0]

sol = A \ b
u_h = sol[1:n_dofs]
u_h = u_h .- mean(u_h)
λ = sol[end]  # Lagrange multiplier
abs(λ)>1e-3 && @warn "Langrange multiplier is large: $(λ)"   # should be very small

## Compute the true solution
u_ex = zeros(length(u_h))
u_ex = apply_analytical!(u_ex, dh, :u, x -> u_chosen(x[1],x[2],x[3]))

## Compute L2 error against true solution
u_ex = u_ex .- mean(u_ex)
err = abs.(u_h .- u_ex)
L2error = compute_L2_error(u_h, u_ex, cellvalues_Ω, dh)

VTKGridFile("surface_laplace_Dzuik_surf_$(refinement)", dh) do vtk
    write_solution(vtk, dh, u_h, "_h")
    write_solution(vtk, dh, err, "_err")
    write_solution(vtk, dh, u_ex, "_true")
end

println("max nodal error    : ", maximum(abs, err))
println("L2 err: ", L2error)
