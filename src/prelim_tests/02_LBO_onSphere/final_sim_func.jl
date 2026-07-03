using Ferrite, FerriteGmsh, LinearAlgebra

include("../../include/gen_sphere.jl")
include("../../include/funcs_error_analysis.jl")
include("surface_cell_values.jl")
include("func_assembleKF_manual.jl")

function solve_laplace_eq_on_sphere(refinement, u_exact::Function, f_rhs::Function, tosavesol::Bool=false; dir_to_sphere_meshes="../../sphere_meshes")
    # --- mesh, quadrature rule, interpolation ---------#
    Ωₕ, h = load_sphere(refinement, dir_to_sphere_meshes)
    qr = QuadratureRule{RefTriangle}(1)
    ip = Lagrange{RefTriangle, 1}()

    # --- dof handler ---------------------------------#
    dh = DofHandler(Ωₕ)
    add!(dh, :u, ip)
    close!(dh)

    cv = SurfaceCellValues(qr, ip)

    # --- assemble stiffness/load, and the constraint (mass) vector -------#
    K, f = assemble_globalKF(dh, cv, (x, y, z) -> f_rhs(x, y, z))
    _, c = assemble_globalKF(dh, cv, (x, y, z) -> 1.0)  # c_i = ∫φᵢ dΩ

    # --- border the singular system with the zero-mean constraint cᵀu = 0 --#
    n = length(f)
    A = [K c; c' 0.0]
    b = [f; 0.0]

    sol = A \ b
    u_h = sol[1:n]
    λ = sol[end]  # Lagrange multiplier

    # --- check against the manufactured solution --------------------------#
    u_ex = zeros(length(u_h))
    u_ex = apply_analytical!(u_ex, dh, :u, x -> u_exact(x[1], x[2], x[3]))
    err = abs.(u_h .- u_ex)
    L2error = compute_L2_error(u_h, u_ex, cv, dh)

    if tosavesol
        # --- export for visualization (e.g. ParaView) ----------------------------#
        VTKGridFile("surface_laplace_sphere_$(refinement)", dh) do vtk
            write_solution(vtk, dh, u_h, "_h")
            write_solution(vtk, dh, err, "_err")
            write_solution(vtk, dh, u_ex, "_true")
        end
    end
    return h, L2error
end

## Example
#=
u_exact(x, y, z) = z
f_rhs(x, y, z) = 2 * z

cd(@__DIR__)
mkpath("results")
cd("results")
solve_laplace_eq_on_sphere(0, u_exact, f_rhs, true)
=#