using Ferrite

include("../../include/func_loadmesh.jl")
refinement = 3
filename = "lowresmesh_$(refinement).msh"
cd(@__DIR__)
cd("../../Dziuk_surf_meshes/")
Ωₕ, h = load_mesh(filename);
@info "Grid imported for refinement level $refinement"
cd(@__DIR__)
#= Ωₕ.nodes; Ωₕ.cells =#
include("custom_interpolation.jl")
include("custom_geomapping.jl")
qr = QuadratureRule{RefTriangle}(2)
ip = LTI()
#cellvalues_Ω = CellValues(qr, ip, ip^3)

dh = DofHandler(Ωₕ)
add!(dh, :u, ip)
close!(dh)
n_dofs = ndofs(dh)
include("func_assembleK.jl")
Kₕ = assemble_globalK_surface(dh, qr, ip)

F = zeros(n_dofs)
include("manufactured_sol.jl")
include("func_assembleRHS.jl")
assemble_rhs_forces!(F, dh, qr, ip, rhs_func)

include("func_assembleconstraintC.jl")
C = zeros(n_dofs)
assemble_constraint_C!(C, dh, qr, ip)

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
l2_err  = L₂_error(u_true_centered, u_h, dh, qr, ip)
#l2_true = L2_norm(u_true_centered, dh, qr, ip)