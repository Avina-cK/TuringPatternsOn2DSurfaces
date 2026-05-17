using Ferrite, FerriteGmsh, Gmsh.gmsh
using SparseArrays
using WriteVTK

#= Simulating the Dziuk surface example from the main paper 
    PDE:    ∇⋅(∇u) = f on Ω,
    where Ω is the Dziuk surface.
    Assumptions: 
        - uₜᵣᵤₑ(x)=x₁⋅x₂
        - f(x) is computed from the PDE using uₜᵣᵤₑ
        - zero mean constraint: ∫udΩ = 0
=#

include("../../include/funcs_error_analysis.jl")
include("funcs_Dziuk.jl")   #functions in the example of the Dziuk paper
include("funcs_assembly.jl") #functions to assemble stiffness matrix and force vector
include("../../include/func_loadmesh.jl") #functions to import mesh into Ferrite.Grid

function sim_heateqn_fromfile(refinement::Int)
    filename = "lowresmesh_$(refinement).msh"
    Ωₕ, h = load_mesh(filename);
    @info "Grid imported for refinement level $refinement"
    
    #-- set up interpolation and quadrature ---#
    ip = Lagrange{RefTriangle, 1}()
    qr = QuadratureRule{RefTriangle}(2)
    cellvalues_Ω = CellValues(qr, ip, ip^3)   #geometric interpolation for 3D space

    #--Degrees of freedom-#
    dh = DofHandler(Ωₕ)
    add!(dh, :u, ip)
    close!(dh);

    #--- Creating zero mean constraint ---------------#
    # ∫udΩ = 0 => Σuᵢ∫φᵢdΩ = 0
    ch = ConstraintHandler(dh);
    # integrate shape functions over the surface to get the sum of dofs
    sum_dofs = zeros(ndofs(dh));
    for cell in CellIterator(Ωₕ)
        reinit!(cellvalues_Ω, cell)
        celldofs = dof_range(dh, :u)
        for q_point in 1:getnquadpoints(cellvalues_Ω)
            dΩ = getdetJdV(cellvalues_Ω, q_point)
            for i in 1:getnbasefunctions(cellvalues_Ω)
                φᵢ = shape_value(cellvalues_Ω, q_point, i)
                sum_dofs[celldofs[i]] += φᵢ .* dΩ
            end
        end
    end
    #--- Creating Affine Constraint by assuming one DoF depends on all the rest--#
    # u₁ = -(Σⱼ₌₂ⁿuⱼ∫φⱼdΩ) / ∫φ₁dΩ
    slave_dof = 1
    # let wₘ[j] = -∫φⱼdΩ / ∫φ₁dΩ
    wₘ = [Pair(j, -sum_dofs[j] / sum_dofs[slave_dof]) for j in 2:ndofs(dh)]

    zero_constraint = AffineConstraint(slave_dof, wₘ, 0.0)
    add!(ch, zero_constraint)
    close!(ch)
    update!(ch, 0.0)

    #--Initiating Stiffness matrices, and RHS vector--#
    K = allocate_matrix(dh, ch)
    f = zeros(ndofs(dh))
    apply!(K, f, ch);
    #--Final system: u = K/F ---#
    K, f = assemble_globalK(cellvalues_Ω, K, dh), assemble_rhs_forces!(f, cellvalues_Ω, dh, rhs_f)
    apply!(K, f, ch);

    #-- Solving system --#
    uₕ = K \ f

    #-- True solution ---------------------------#
    uₜᵣᵤₑ = zeros(ndofs(dh))
    apply_analytical!(uₜᵣᵤₑ, dh, :u, x -> x[1]*x[2])

    #-- Error analysis -------------------------#
    L2error= compute_L2_error(uₕ, uₜᵣᵤₑ, cellvalues_Ω, dh)
    rL2error = relative_L2_error(uₕ, uₜᵣᵤₑ, cellvalues_Ω, dh)
    return h, L2error, rL2error
end

cd(@__DIR__)
cd("../../Dziuk_surf_meshes/")
rs = 0:2
nrs = length(rs)
hs = zeros(nrs)
L2errors = zeros(nrs)
rL2errors = zeros(nrs)

for i in 1:nrs
    r = rs[i]
    hᵣ, L2errorᵣ, rL2errorᵣ = sim_heateqn_fromfile(r);
    #print(hᵣ, L2errorᵣ, rL2errorᵣ)
    hs[i] = hᵣ
    L2errors[i] = L2errorᵣ
    rL2errors[i] = rL2errorᵣ
end

EOCs = zeros(nrs)
for i in 1:nrs-1
    EOCs[i] = EOC(hs[i], hs[i+1], L2errors[i], L2errors[i+1])
end
EOCs[nrs] = NaN

data = [hs L2errors rL2errors EOCs]

using PrettyTables
pretty_table(data;
            column_labels = ["h", "L₂ error", "Relative L₂ error", "EOC"]
            )

using Plots
plot(log.(data[:,1]), log.(data[:,2])
    , xaxis = "log(h)"
    , yaxis="log(L₂-error)"
    , marker=true
    , legend=false)
