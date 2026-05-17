# Transient heat equation: ∂u/∂t - ∇⋅(k∇u) = f in domain
# with u(x,t) = 0 on the boundary

function initial_setup(refine)
        
    #--- Creating domain ----------------------#
    # left most point
    left_pt = Vec{2, Float64}([-L₁, -L₂])
    # right most point
    right_pt = Vec{2, Float64}([L₁, L₂]) 
    # generate a rectangular grid, with triangules, 'refine' divisions in either directions, between left_pt and right_pt
    Ωₕ = generate_grid(Triangle, (refine, refine), left_pt, right_pt)
    # @info "Generated grid, refinement=$(refine)"
    # boundaries
    ∂Ω = union(getfacetset.((Ωₕ,),["top", "bottom", "left", "right"])...)

    #--- Test and Trial functions ------------#
    ip = Lagrange{RefTriangle, 1}()
    qr = QuadratureRule{RefTriangle}(2)
    cellvalues = CellValues(qr, ip)

    #--- Degrees of Freedom -----------------#
    dh = DofHandler(Ωₕ)
    add!(dh, :u, ip)
    close!(dh)

    #--- Initiating Mass and Stiffness matrices, and RHS -------#
    K = allocate_matrix(dh)
    M = allocate_matrix(dh)
    F = zeros(ndofs(dh))
    # @info "Initialised empty mass and stiffness matrices, and RHS vector"

    #--- Boundary conditions -----------------------------------#
    # Condition on ∂Ω: u=0
    dbc = Dirichlet(:u, ∂Ω, (x,t) -> 0);

    # Create ConstraintHandler and add boundary constraints
    ch = ConstraintHandler(dh);
    add!(ch, dbc);
    close!(ch);
    # @info "Set up dirichlet boundary conditions"
    return Ωₕ, ∂Ω, cellvalues, dh, K, M, F, ch
end

#--- Create initial force vector, from true solution ------------#
function assemble_rhs_forces!(F::Vector, cellvalues::CellValues, dh::DofHandler, func_rhs_assemble::Function)
    fill!(F, 0.0);
    n_basefuncs = getnbasefunctions(cellvalues)
    # Create a temporary array for the cell's local contributions to the external force vector
    fe = zeros(n_basefuncs)
    for cell in CellIterator(dh)
        # Update the cellvalues to the correct cell number
        reinit!(cellvalues, cell)
        # Reset the temporary array for the next cell
        fill!(fe, 0.0)
        # Access the cell's coordinates
        cell_coordinates = getcoordinates(cell)
        for q_point in 1:getnquadpoints(cellvalues)
            # Calculate the global coordinate of the quadrature point.
            x = spatial_coordinate(cellvalues, q_point, cell_coordinates)
            fᵥ = func_rhs_assemble(x)
            # Get the integration weight for the current quadrature point.
            dΩᵣ = getdetJdV(cellvalues, q_point)
            for i in 1:n_basefuncs
                φᵢ = shape_value(cellvalues, q_point, i)
                fe[i] += fᵥ ⋅ φᵢ * dΩᵣ
            end
        end
        # Add the local contributions to the correct indices in the global external force vector
        assemble!(F, celldofs(cell), fe)
    end
    return F
end

#--- Assembly functions for K and M -----------------------------#
# Assembling local K, M contributions
function assemble_KeMe!(Kₑ::Matrix, Mₑ::Matrix, cellvalues::CellValues)
    n_basefuncs = getnbasefunctions(cellvalues)
    # reset to 0
    fill!(Kₑ, 0)
    fill!(Mₑ, 0)
    # loop over quadrature points
    for q_point in 1:getnquadpoints(cellvalues)
        # quadrature weight:
        dΩᵣ = getdetJdV(cellvalues, q_point)
        # loop over test functions
        for i in 1:n_basefuncs
            φᵢ = shape_value(cellvalues, q_point, i)
            ∇φᵢ = shape_gradient(cellvalues, q_point, i)
            # loop over trial functions
            for j in 1:n_basefuncs
                φⱼ = shape_value(cellvalues, q_point, j)
                ∇φⱼ = shape_gradient(cellvalues, q_point, j)
                Kₑ[i,j] += k * (∇φᵢ ⋅ ∇φⱼ) * dΩᵣ
                Mₑ[i,j] += (φᵢ * φⱼ) *dΩᵣ
            end
        end
    end
    return Kₑ, Mₑ
end

#--- Assemble global matrices --------------------------#
"""
assemble_KM!(K::SparseMatrixCSC, M::SparseMatrixCSC, cellvalues::CellValues, dh::DofHandler)
Assemble global stiffness and mass matrices by
    1. looping over all cells,
    2. computing local contributions, and
    3. adding them to the global matrices.
"""
function assemble_KM!(K::SparseMatrixCSC, M::SparseMatrixCSC, cellvalues::CellValues, dh::DofHandler)
    n_basefuncs = getnbasefunctions(cellvalues)
    Ke = zeros(n_basefuncs, n_basefuncs)
    Me = zeros(n_basefuncs, n_basefuncs)
    assembler_K = start_assemble(K)
    assembler_M = start_assemble(M)

    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        assemble_KeMe!(Ke, Me, cellvalues)
        assemble!(assembler_K, celldofs(cell), Ke)
        assemble!(assembler_M, celldofs(cell), Me)
    end
    return K, M
end

#--- Final system -------------------------------------------#

function sim_for_error(refine, Tf, u_true::Function, func_rhs2::Function)
    Ωₕ, ∂Ω, cellvalues, dh, K, M, F, ch = initial_setup(refine);
    # Assembling K, M, F
    K, M = assemble_KM!(K, M, cellvalues, dh)
    A = (Δt .* K) + M
    apply!(A, ch) 
    RHS_data = get_rhs_data(ch, A)
    
    # initiate initial state
    # Initial conditions: u₀ = u_true(x, 0)
    uₙ = zeros(ndofs(dh))
    apply_analytical!(uₙ, dh, :u, x -> u_true(x, 0.0))

    for t in Δt:Δt:Tf
        # update Dirichlet boundary values
        update!(ch, t)
        # Compute RHS: b = Δt*F + Muₙ
        func_rhs_t(x) = func_rhs2(x,t)
        F = assemble_rhs_forces!(F, cellvalues, dh, func_rhs_t)
        b = Δt .* F .+ M * uₙ
        # Apply boundary conditions for this time step
        apply_rhs!(RHS_data, b, ch)
        apply!(A, ch)
        # Solve the time test
        uₙ₊₁ = A \ b
        # L2 error from previous time step
        change_in_u = compute_L2_error(uₙ, uₙ₊₁, cellvalues, dh)
        if change_in_u < 1e-6
            @info "Solution converged at time t=$(t) w.r.t. L-2 change in u = $(change_in_u). Stopping simulation."
            Tf = t;
            break
        end
        # Update uₙ for the next iteration
        uₙ .= uₙ₊₁
    end

    #--- Error Analysis ----------------------------------------------------------#
    # true solution
    uₜᵣᵤₑ = zeros(ndofs(dh))
    apply_analytical!(uₜᵣᵤₑ, dh, :u, x -> u_true(x, Tf))

    # Vectors: uₜᵣᵤₑ and uₙ
    return compute_L2_error(uₙ, uₜᵣᵤₑ, cellvalues, dh), relative_L2_error(uₙ, uₜᵣᵤₑ, cellvalues, dh), Tf
end


function error_table(refinements, Tf, u_true::Function, func_rhs::Function)
    rₛ = length(refinements)
    L2_errors = zeros(rₛ)
    rel_L2_errors = zeros(rₛ)
    T_finals = zeros(rₛ)

    for i in 1:rₛ
        L2_errors[i], rel_L2_errors[i], T_finals[i] = sim_for_error(refinements[i], Tf, u_true, func_rhs)
    end

    h_j = 2.0 ./refinements
    EOC_succesive = zeros(rₛ)
    for i in 2:rₛ
        EOC_succesive[i] = EOC(h_j[i-1], h_j[i], L2_errors[i-1], L2_errors[i])
    end
    EOC_succesive[1] = NaN;
    return [h_j L2_errors rel_L2_errors EOC_succesive T_finals]
end

@info "Loaded setup functions for the transient heat equation with Dirichlet boundary conditions:\n ∂u/∂t - ∇⋅(k∇u) = f \n on Ω = [-$(L₁), $(L₁)] × [-$(L₂), $(L₂)], \n with k = $(k), Δt = $(Δt)."
