#=
PDE:
    ∂ₜu - Dᵤ⋅Δₛu = α(1-u) - uv²
    ∂ₜv - Dᵥ⋅Δₛv = -(α+β)v + uv²
=#

using Ferrite, FerriteGmsh
using BlockArrays, SparseArrays, LinearAlgebra
using Random
using WriteVTK

include("../../include/surface_cell_values.jl")
include("../../include/gen_read_meshes/gen_sphere.jl")
include("../../include/struct_GrayScott.jl")

include("func_assemble_matrices.jl")
include("func_assemble_initial_cond.jl")


given_m_default = GrayScottMaterial{Float64}(0.00016, 0.00008, 0.06, 0.062)
dir_to_sphere_meshes_default = "../../sphere_meshes/"

function sim_gs_on_sphere(given_m::GrayScottMaterial{Float64}=given_m_default,
                        refinement::Int=3, Δt::Float64=1.0, T::Float64=50000.0,
                        savesol::Bool=true, no_frames_to_save::Int=100,
                        dir_to_sphere_meshes::String=dir_to_sphere_meshes_default
    )
        
    nsteps = ceil(Int, T/Δt) 
    save_every = Int(floor(nsteps/no_frames_to_save))
    Du, Dv = given_m.Dᵤ, given_m.Dᵥ
    α, β = given_m.α, given_m.β

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

    assemble_matrices!(M, D, cellvalues_Ω, dh, given_m)

    A = M + Δt .* D
    cholA = cholesky(A)

    ## initial condition
    uₜ = zeros(ndofs(dh))
    uₜ₋₁ = zeros(ndofs(dh))
    setup_initial_conditions!(uₜ₋₁, cellvalues_Ω, dh)

    ## Save solutions
    if savesol
        cd(@__DIR__)
        fts(x) = replace(string(x), "." => "-") #float to string
        basefolder_name = "GS_Sphere_$(refinement)_$(fts(Du))_$(fts(Dv))_$(fts(α))_$(fts(β))"
        mkpath("results/$(basefolder_name)")
        cd("results/$(basefolder_name)/")
        basefile_name = "GS_Sphere"
        pvd = paraview_collection("$(basefile_name)")
        VTKGridFile("$(basefile_name)-0", dh) do vtk
            write_solution(vtk, dh, uₜ₋₁, "_h")
            pvd[0.0] = vtk
        end
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
            rvₜ[1, i] += Δt * reaction1(u, v, α, β)
            rvₜ[2, i] += Δt * reaction2(u, v, α, β)
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
    @info "Gray-Scott simulation completed on sphere."
    @info "Refinement: $(refinement)"
    @info "Parameters: Du=$(Du), Dv=$(Dv), α=$(α), β=$(β)"
    cd(@__DIR__)
end


given_m_2 = GrayScottMaterial{Float64}(0.00016, 0.00008, 0.064, 0.062)

sim_gs_on_sphere(given_m_2)