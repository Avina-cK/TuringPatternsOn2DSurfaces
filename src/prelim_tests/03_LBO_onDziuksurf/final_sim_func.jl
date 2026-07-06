using LinearAlgebra
using Ferrite, Tensors, SparseArrays
using WriteVTK

include("surface_cell_values.jl")
include("../../include/read_Dzuik_mesh.jl")
include("../../include/funcs_error_analysis.jl")
include("../../include_wip/check_mesh.jl")

include("func_assembleKF.jl")
include("func_assemble_u_true.jl")
include("manufactured_sol.jl")
include("func_buildlagrangesys.jl")


function solve_LB_eq_on_Dzuiksurf(refinement)

    filename = "Dzuik_mesh_$(refinement).msh"
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

    @assert check_normal_orientation(dh, cellvalues_Ω) == [] "Some cells are not oriented correctly. Check your mesh!"

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
    A, b = build_lagrange_system(Kₕ, F, c)

    sol = A \ b
    u_h = sol[1:n_dofs]
    #u_h = u_h .- mean(u_h)
    λ = sol[end]  # Lagrange multiplier
    abs(λ)>1e-3 && @warn "Langrange multiplier is large: $(λ)"
    # should be very small, according to Dzuik (2013)

    ## Compute the true solution
    u_ex = assemble_trueU_nodal(dh, cellvalues_Ω, u_chosen)

    ## Compute L2 error against true solution
    u_err = u_h .- u_ex
    #u_n_err = normalize(u_h) .- normalize(u_ex)

    L2error = compute_L2_error(u_h, u_ex, cellvalues_Ω, dh)

    ## Write output
    VTKGridFile("surface_laplace_Dzuik_surf_$(refinement)", dh) do vtk
        write_solution(vtk, dh, u_h, "_h")
        write_solution(vtk, dh, u_err, "_err")
        #write_solution(vtk, dh, u_n_err, "_n_err")
        write_solution(vtk, dh, u_ex, "_true")
    end
        #print([refinement, h, L2error])

    return refinement, h, L2error
end
