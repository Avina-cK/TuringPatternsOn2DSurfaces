#=
PDE:
    ∂ₜu - Dᵤ⋅Δₛu = α(1-u) - uv²
    ∂ₜv - Dᵥ⋅Δₛv = -(α+β)v + uv²
=#

using Ferrite, FerriteGmsh
using BlockArrays, SparseArrays, LinearAlgebra
using Random
#using WriteVTK
using VTKHDF

include("../../include/surface_cell_values.jl")
include("../../include/gen_read_meshes/func_load_tri_surfmesh.jl")
include("../../include/struct_GrayScott.jl")

include("../../include/fem_assembly_funcs/func_assemble_matrices_RDE.jl")
include("../../include/fem_assembly_funcs/func_assemble_polar_initial_cond.jl")
include("../../include/fem_assembly_funcs/func_assemble_const_initial_cond.jl")

include("save_load_state.jl")

given_m_default = GrayScottMaterial{Float64}(0.00016, 0.00008, 0.06, 0.062)
dir_to_meshes_default = "../../Dziuk_surf_meshes"

function fts(x, n)
    # Small values: format as *-**e-*
    if x != 0 && abs(x) < 1.0e-4
        e = floor(Int, log10(abs(x)))
        m = round(x / 10.0^e, digits=n)
        if abs(m) >= 10 # rounding pushed mantissa to 10.0
            m /= 10
            e += 1
        end
        s = string(m)
        i = findfirst('.', s)
        s = rpad(s, i + n, '0')
        return replace(s, "." => "-") * "e" * string(e)
    end

    # Old version (if breaks)
    s = string(round(x, digits=n))
    i = findfirst('.', s)
    s = rpad(s, i + n, '0')
    return replace(s, "." => "-")
end

function sim_gs_on_Dziukmesh(given_m::GrayScottMaterial{Float64}=given_m_default,
                        refinement::Int=4, Δt::Float64=1.0, T::Float64=50000.0,
                        savesol::Bool=true, no_frames_to_save::Int=100,
                        dir_to_meshes::String=dir_to_meshes_default;
                        restart_from::Union{Nothing,String}=nothing,
                        init_const=[0.0; 0.0], noise_amp_given::Float64=0.01
    )
        
    nsteps = ceil(Int, T/Δt) 
    save_every = Int(floor(nsteps/no_frames_to_save))
    Du, Dv = given_m.Dᵤ, given_m.Dᵥ
    α, β = given_m.α, given_m.β
    cd(@__DIR__)
    cd(dir_to_meshes)
    filename="Dzuik_mesh_$(refinement).msh"
    Ωₕ, h = load_tri_surfmesh(filename)
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
    restarted = restart_from !== nothing
    if restarted
        cd(@__DIR__)
        uₜ₋₁, t_start, iter_start = load_state(restart_from)
        @assert length(uₜ₋₁) == ndofs(dh) "Loaded state does not match the number of degrees of freedom in the mesh."
        @info "Restarting from $(restart_from) at t=$(t_start), iter=$(iter_start)"
    else
        t_start=0.0
        iter_start=0
        if init_const==[0.0; 0.0]
            setup_polar_initial_conditions!(uₜ₋₁, cellvalues_Ω, dh; noise_amp = noise_amp_given)
        else
            setup_const_init_conditions!(uₜ₋₁, cellvalues_Ω, dh, init_const)
        end
    end

    ## Save solutions
    if savesol
        cd(@__DIR__)
        basefolder_name = "GS_Dziuk_surf_$(refinement)_$(fts(Du,5))_$(fts(Dv,5))_$(fts(α,3))_$(fts(β,6))"
        mkpath("results/$(basefolder_name)")
        cd("results/$(basefolder_name)/")
        basefile_name = "GS_Dziuk_surf"
        open(basefolder_name*".txt", "a") do io
            write(io, "Δt = $(Δt), D_u = $(fts(Du,10)), D_v = $(fts(Dv,10)), alpha = $(fts(α,3)), beta = $(fts(β,6))")
        end
        vtkhdf = VTKHDFGridFile("$(basefile_name).vtkhdf", dh; temporal = true)
        #pvd = paraview_collection("$(basefile_name)"; append=restarted)

        if !restarted
            write_timestep(vtkhdf, 0.0) do vtk
                write_solution(vtk, dh, uₜ₋₁, "_h")
            end
            #VTKGridFile("$(basefile_name)-0", dh) do vtk
            #    write_solution(vtk, dh, uₜ₋₁, "_h")
            #    pvd[0.0] = vtk
            #end
        end
    end

    ## main loop
    iₜ_final = iter_start
    t_final = t_start

    for (iₜ_r, t_r) in enumerate(Δt:Δt:T)
        iₜ = iter_start + iₜ_r
        t = t_start + t_r
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
        iₜ_final, t_final = iₜ, t
        if iₜ_r % save_every == 0
            write_timestep(vtkhdf, t) do vtk
                write_solution(vtk, dh, uₜ, "_h")
            end
            #VTKGridFile("$(basefile_name)-$(iₜ)", dh) do vtk
            #    write_solution(vtk, dh, uₜ, "_h")
            #    pvd[t] = vtk
            #end
            @info "Step $(t)/$(t_start + T)"
        end
        if any(isnan.(uₜ))
            @error "NaN detected in solution at step $(iₜ), t = $(t)"
            break
        end
        if isapprox(uₜ, uₜ₋₁; rtol=1e-8, atol=1e-10)
            @info "Converged at step $(iₜ), t = $(t)"
            write_timestep(vtkhdf, t) do vtk
                write_solution(vtk, dh, uₜ, "_h")
            end
            #VTKGridFile("$(basefile_name)-$(iₜ)", dh) do vtk
            #    write_solution(vtk, dh, uₜ, "_h")
            #    pvd[t] = vtk
            #end
            break
        end
        uₜ₋₁ .= uₜ
    end
    if savesol
        close(vtkhdf)
        #close(pvd)
        # always persist final state for a future restart, regardless of exit reason
        save_state(uₜ₋₁, t_final, iₜ_final, "$(basefile_name)_state_$(iₜ_final).jld2")
    end
    @info "Gray-Scott simulation completed on Dziuk surface."
    @info "Refinement: $(refinement)"
    @info "Parameters: Du=$(Du), Dv=$(Dv), α=$(α), β=$(β)"
    cd(@__DIR__)
end