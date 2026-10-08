using Ferrite, FerriteGmsh
using BlockArrays, SparseArrays, LinearAlgebra
using Random
using WriteVTK

include("../../include/struct_Schankenberg.jl")

const a = 0.2
const b = 1.0

# steady state
u0 = steadystate(a,b)[1]
v0 = steadystate(a,b)[2]
dc = d_critical(a,b)


include("../../include/surface_cell_values.jl")
include("../../include/gen_read_meshes/gen_sphere.jl")
include("../../include/fem_assembly_funcs/func_assemble_matrices_RDE.jl")
include("func_assemble_initial_cond.jl")

##
#= From paper: 
    Spatio-temporal pattern formation on spherical surfaces: 
    numerical simulation and application to solid tumour growth
    M.A.J. Chaplain · M. Ganesh · I.G. Graham.
    J. Math. Biol. 42, 387–423 (2001)
    DOI: 10.1007/s002850000067
=#

refinement::Int=3
Δt::Float64=0.01
T::Float64=10.0
savesol::Bool=true
no_frames_to_save::Int=100

## Experiment 1.1
k_sq = 6
γ = γ_critical(a,b, sqrt(k_sq))

d = 15
Dᵤ = 1.0
Dᵥ = Dᵤ*d

given_m = SchnakenbergMaterial{Float64}(Dᵤ, Dᵥ, a, b)
dir_to_sphere_meshes = "../../sphere_meshes/"

nsteps = ceil(Int, T/Δt) 
save_every = Int(floor(nsteps/no_frames_to_save))
Du, Dv = given_m.Dᵤ, given_m.Dᵥ
α, β = given_m.a, given_m.b

cd(@__DIR__)
filename="$(dir_to_sphere_meshes)/embedded_sphere_$(refinement).msh"
if isfile(filename)
    Ωₕ, h = load_sphere(refinement, dir_to_sphere_meshes)
else
    Ωₕ, h = create_embedded_sphere(refinement)
end

## Cell values
ip = Lagrange{RefTriangle, 1}()
qr = QuadratureRule{RefTriangle}(1)
cellvalues_Ω = SurfaceCellValues(qr, ip)

dh = DofHandler(Ωₕ)
add!(dh, :uv, ip^2)
close!(dh)

## assemble matrices
M = allocate_matrix(dh; coupling=[true false; false true])
D = allocate_matrix(dh; coupling=[true false; false true])

assemble_matrices_RDE!(M, D, cellvalues_Ω, dh, given_m)

A = M + Δt .* D
cholA = cholesky(A)

## initial condition
uₜ = zeros(ndofs(dh))
uₜ₋₁ = zeros(ndofs(dh))
setup_random_initcond!(uₜ₋₁, cellvalues_Ω, dh, u0, v0;noise_amp=0.0001)

cd(@__DIR__)
fts(x) = replace(string(x), "." => "-") #float to string
basefolder_name = "Eg1-1_Schnakenberg_Sphere_$(refinement)_$(fts(Du))_$(fts(Dv))_$(fts(γ))_$(fts(α))_$(fts(β))"
mkpath("results/$(basefolder_name)")
cd("results/$(basefolder_name)/")
basefile_name = "Schnakenberg_Sphere"
pvd = paraview_collection("$(basefile_name)")
VTKGridFile("$(basefile_name)-0", dh) do vtk
    write_solution(vtk, dh, uₜ₋₁, "_h")
    pvd[0.0] = vtk
end

## main loop

for (iₜ, t) in enumerate(Δt:Δt:T)

    # solve main heat problem ()
    uₜ .= cholA \ (M * uₜ₋₁)

    uv_dofs = ndofs(dh) ÷ 2
    rvₜ = reshape(uₜ, (2, uv_dofs))
    for i in 1:uv_dofs
        u = rvₜ[1, i]
        v = rvₜ[2, i]
        rvₜ[1, i] += Δt * reaction1(u, v, α, β, γ)
        rvₜ[2, i] += Δt * reaction2(u, v, α, β, γ)
    end
    if iₜ % save_every == 0
        VTKGridFile("$(basefile_name)-$(iₜ)", dh) do vtk
            write_solution(vtk, dh, uₜ, "_h")
            pvd[t] = vtk
        end
        @info "Step $(t)/$(T)"
    end
    if any(isnan.(uₜ))
        @error "NaN detected in solution at step $(iₜ), t = $(t)"
        break
    end
    if isapprox(uₜ, uₜ₋₁; rtol=1e-8, atol=1e-10)
        @info "Converged at step $(iₜ), t = $(t)"
        break
    end
    uₜ₋₁ .= uₜ
end
close(pvd)

cd(@__DIR__)