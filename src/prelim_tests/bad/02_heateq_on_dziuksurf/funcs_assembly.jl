#-- Assembly of RHS ---------------------------------------------------------------------------------------#
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

# Assembling contributions of stiffness matrix for finite elements
function assemble_elementK!(Ke::Matrix, cellvalues::CellValues)
    n_base_funcs = getnbasefunctions(cellvalues);
    # Reset matrix and vector
    fill!(Ke, 0)
    no_quadpoints = getnquadpoints(cellvalues);
    #Loop over quadrature points
    for q_point in 1:no_quadpoints
        # get quadrature weight × det(J)
        dΩᵣ = getdetJdV(cellvalues, q_point);
        #Loop over test shape functions
        for i in 1:n_base_funcs
            φᵢ = shape_value(cellvalues, q_point, i)
            ∇φᵢ = shape_gradient(cellvalues, q_point, i)
            #Loop over trial shape functions
            for j in 1:n_base_funcs
                ∇φⱼ = shape_gradient(cellvalues, q_point, j)
                #Add to Ke
                Ke[i,j] += (∇φᵢ ⋅ ∇φⱼ) * dΩᵣ
            end
        end
    end
    return Ke
end

function assemble_globalK(Kₕ::SparseMatrixCSC, cellvalues::CellValues, dh::DofHandler)
    n_base_funcs = getnbasefunctions(cellvalues);
    # Initiate element stiffness matrix and element force vector
    Ke = zeros(n_base_funcs, n_base_funcs);
    # Create assembler: allocates workspace necessary for efficient matrix assembly
    assembler = start_assemble(Kₕ);
    # Loop over all cells
    for cell in CellIterator(dh)
        # reinitalise cellvalues for cell, to update derivatives of shape functions etc
        reinit!(cellvalues, cell);
        # Compute element contribution
        assemble_elementK!(Ke, cellvalues);
        # Assemble constributions of Ke into Kₕ
        #   equivalent to K[dofs, dofs] += Ke, 
        #   where K is the global stiffness matrix 
        assemble!(assembler, celldofs(cell), Ke)
    end
    return Kₕ
end
@info "Loaded assembly functions for stiffness matrix and force vector"