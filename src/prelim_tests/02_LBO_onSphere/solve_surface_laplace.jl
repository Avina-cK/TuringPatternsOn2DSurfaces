using Ferrite, FerriteGmsh, LinearAlgebra

include("../../include/gen_sphere.jl")
include("surface_cell_values.jl")
include("func_assembleKF_manual.jl")

# --- mesh, quadrature rule, interpolation ---------#
#TODO: store meshes in files 
Ωₕ = create_embedded_sphere(0)
qr = QuadratureRule{RefTriangle}(1)
ip = Lagrange{RefTriangle, 1}()

# --- dof handler ---------------------------------#
dh = DofHandler(Ωₕ)
add!(dh, :u, ip)
close!(dh)

cv = SurfaceCellValues(qr, ip)

# --- problem data --------------------------------#
u_exact(x, y, z) = z
f_rhs(x, y, z) = 2 * z

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

U = assemble_integral_global_uh(dh, cv, u_h)

# --- check against the manufactured solution --------------------------#
#u_ex = [node.x[3] for node in Ωₕ.nodes]
# TODO: why is it not in nodal order?
# https://ferrite-fem.github.io/Ferrite.jl/stable/topics/degrees_of_freedom/#DoF-numbering-is-decoupled-from-the-node-numbering-f75c9a3d32993a20

u_ex = zeros(length(u_h))
u_ex = apply_analytical!(u_ex, dh, :u, x -> x[3])
err = u_h .- u_ex
println("nnodes            : ", n)
println("max nodal error    : ", maximum(abs, err))
println("relative L2-ish err: ", norm(err) / norm(u_ex))

# --- export for visualization (e.g. ParaView) ----------------------------#
VTKGridFile("surface_laplace_sphere", dh) do vtk
    write_solution(vtk, dh, u_h)
end
VTKGridFile("surface_laplace_sphere_err", dh) do vtk
    write_solution(vtk, dh, err)
end

VTKGridFile("surface_laplace_sphere_ex", dh) do vtk
    write_solution(vtk, dh, u_ex)
end
