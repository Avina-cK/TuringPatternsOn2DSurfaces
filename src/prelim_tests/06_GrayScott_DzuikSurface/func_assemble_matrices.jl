function assemble_element_mass!(Me::Matrix, cellvalues)
    n_basefuncs = getnbasefunctions(cellvalues)
    num_reactants = 2
    r₁range = 1:num_reactants:(num_reactants * n_basefuncs)
    r₂range = 2:num_reactants:(num_reactants * n_basefuncs)
    Me₁ = @view Me[r₁range, r₁range]
    Me₂ = @view Me[r₂range, r₂range]
    # Loop over quadrature points
    for q_point in 1:getnquadpoints(cellvalues)
        # Get the quadrature weight
        dΩ = getdetJdV(cellvalues, q_point)
        # Loop over test shape functions
        for i in 1:n_basefuncs
            δuᵢ = shape_value(cellvalues, q_point, i)
            # Loop over trial shape functions
            for j in 1:n_basefuncs
                δuⱼ = shape_value(cellvalues, q_point, j)
                # Add contribution to Ke
                Me₁[i, j] += (δuᵢ * δuⱼ) * dΩ
                Me₂[i, j] += (δuᵢ * δuⱼ) * dΩ
            end
        end
    end
    return nothing
end
function assemble_element_diffusion!(De::Matrix, cellvalues, material::GrayScottMaterial)
    n_basefuncs = getnbasefunctions(cellvalues)
    Dᵤ = material.Dᵤ
    Dᵥ = material.Dᵥ
    # The diffusion between the reactions is not coupled, so we get a blocked-strided matrix.
    num_reactants = 2
    r₁range = 1:num_reactants:(num_reactants * n_basefuncs)
    r₂range = 2:num_reactants:(num_reactants * n_basefuncs)
    De₁ = @view De[r₁range, r₁range]
    De₂ = @view De[r₂range, r₂range]
    # Loop over quadrature points
    for q_point in 1:getnquadpoints(cellvalues)
        # Get the quadrature weight
        dΩ = getdetJdV(cellvalues, q_point)
        # Loop over test shape functions
        for i in 1:n_basefuncs
            ∇δuᵢ = shape_gradient(cellvalues, q_point, i)
            # Loop over trial shape functions
            for j in 1:n_basefuncs
                ∇δuⱼ = shape_gradient(cellvalues, q_point, j)
                # Add contribution to Ke
                De₁[i, j] += Dᵤ * (∇δuᵢ ⋅ ∇δuⱼ) * dΩ
                De₂[i, j] += Dᵥ * (∇δuᵢ ⋅ ∇δuⱼ) * dΩ
            end
        end
    end
    return nothing
end
function assemble_matrices!(M::SparseMatrixCSC, D::SparseMatrixCSC, cellvalues, dh::DofHandler, material::GrayScottMaterial)
    n_basefuncs = getnbasefunctions(cellvalues)
    # Allocate the element stiffness matrix and element force vector
    Me = zeros(2 * n_basefuncs, 2 * n_basefuncs)
    De = zeros(2 * n_basefuncs, 2 * n_basefuncs)
    # Create an assembler
    M_assembler = start_assemble(M)
    D_assembler = start_assemble(D)
    # Loop over all cells
    for cell in CellIterator(dh)
        # Reinitialize cellvalues for this cell
        reinit!(cellvalues, cell)
        fill!(Me, 0)
        fill!(De, 0)
        # Compute element contribution
        assemble_element_mass!(Me, cellvalues)
        assemble!(M_assembler, celldofs(cell), Me)
        assemble_element_diffusion!(De, cellvalues, material)
        assemble!(D_assembler, celldofs(cell), De)
    end
    return nothing
end;