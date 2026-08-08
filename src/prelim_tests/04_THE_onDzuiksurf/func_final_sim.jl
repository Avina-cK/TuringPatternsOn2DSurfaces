# main library
using Ferrite
# storage and solvers
using LinearAlgebra, LDLFactorizations, SparseArrays, Tensors

# to load and read mesh
using FerriteGmsh
using Gmsh.gmsh

using WriteVTK # for output

include("../../include/surface_cell_values.jl")
include("../../include/gen_read_meshes/func_load_tri_surfmesh.jl")
include("../../include/funcs_error_analysis.jl")

include("../../include/fem_assembly_funcs/func_assembleF.jl")
include("../../include/fem_assembly_funcs/func_assembleK.jl")
include("../../include/fem_assembly_funcs/func_assembleM.jl")

include("../../include/fem_assembly_funcs/func_assembleU.jl")

include("manufactured_sol.jl")

#TODO: docstring
"""
    sim_THE_onDziuksurf(refinement::Int, T::Float64, Δt::Float64, k_given::Float64=1.0; tosavesol::Bool=false)
    
Simulate the following PDE on a Dzuik surface mesh:
    ∂ₜu - Δₛu = f on Ω
"""
function sim_THE_onDziuksurf(refinement::Int, T::Float64, Δt::Float64, k_given::Float64=1.0;
                             tosavesol::Bool=false)
    filename = "Dzuik_mesh_$(refinement).msh"
    cd(@__DIR__)
    cd("../../Dziuk_surf_meshes/")
    Ωₕ, h = load_tri_surfmesh(filename)

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

    ## Assemble global stiffness and mass matrices, Kₕ and Mₕ
    n_dofs = ndofs(dh)
    Kₕ = assemble_globalK(cellvalues_Ω, dh; k=k_given)
    Mₕ = assemble_globalM(cellvalues_Ω, dh)
    #@assert norm(Kₕ * ones(n_dofs))≈0 "Kₕ 𝕀 = $(norm(Kₕ * ones(n_dofs))) ≠ 0"

    # initial load vector
    rhs_func_0(x,y,z) = rhs_func(x,y,z,0.0,k_given)
    Fₕ = assemble_globalF(cellvalues_Ω, dh, rhs_func_0)
    C = assemble_globalF(cellvalues_Ω, dh, (x,y,z) -> 1.0)
    
    # initial conditions
    uₙ = zeros(length(Fₕ));
    apply_analytical!(uₙ, dh, :u, x -> u_chosen(x[1], x[2], x[3], 0.0));
    @assert abs.(mean(uₙ))<1e-15

    ## Compute the true solution at t = 0 for the initial error field
    u_ex0 = assemble_U(cellvalues_Ω, dh, (x,y,z) -> u_chosen(x,y,z, 0.0))
    u_err0 = uₙ .- u_ex0
    @assert norm(u_err0)==0
    
    ##
    if tosavesol
        cd("results/")
        pvd = paraview_collection("THE_Dzuiksurf_$(refinement)")
        VTKGridFile("THE_Dzuiksurf_$(refinement)-0", dh) do vtk
            write_solution(vtk, dh, uₙ, "_h")
            write_solution(vtk, dh, u_err0, "_err")
            pvd[0.0] = vtk
        end
    end
    Aₕ = Mₕ + Δt .* Kₕ
    Aug = [ Aₕ  C;
            C'  0.0]
    Aug = sparse(Aug)
    Aₕ_F = ldl(Aug)

    for (step, t) in enumerate(Δt:Δt:T)
        
        Fₕ = assemble_globalF(cellvalues_Ω, dh, (x,y,z) -> rhs_func(x,y,z,t,k_given))
        b = Δt .* Fₕ + Mₕ * uₙ
        rhs = [b; 0.0]

        sol = Aₕ_F \ rhs
        uₙ₊₁ = sol[1:n_dofs]
        #λ    = sol[end]
        uₙ .= uₙ₊₁

        # true solution and error at this timestep
        u_ex_t = assemble_U(cellvalues_Ω, dh, (x,y,z) -> u_chosen(x,y,z, t))
        @assert abs.(mean(u_ex_t))<1e-15
        u_err_t = uₙ .- u_ex_t

        if tosavesol
            VTKGridFile("THE_Dzuiksurf_$(refinement)-$(step)", dh) do vtk
                write_solution(vtk, dh, uₙ, "_h")
                write_solution(vtk, dh, u_err_t, "_err")
                pvd[t] = vtk
            end
        end
    end
    
    if tosavesol
        vtk_save(pvd)
    end

    ## Final solution, true solution, and error
    u_h = uₙ
    u_ex = assemble_U(cellvalues_Ω, dh, (x,y,z) -> u_chosen(x,y,z, T))
    u_err = u_h .- u_ex

    ## Compute L2 error against true solution
    L2error = compute_L2_error(u_h, u_ex, cellvalues_Ω, dh)

    if tosavesol
        ## Write final output (separate file, includes all three fields)
        VTKGridFile("THE_Dzuiksurf_$(refinement)_final", dh) do vtk
            write_solution(vtk, dh, u_h, "_h")
            write_solution(vtk, dh, u_err, "_err")
            write_solution(vtk, dh, u_ex, "_true")
        end
        cd(@__DIR__)
    end

    return refinement, h, L2error
end

## Example: refinement=2, Time to run simulation=2π, Δt=0.04, k=1.0, tosavesol=true
#sim_THE_onDziuksurf(2, 2.0*π, 0.04, 1.0; tosavesol=true)