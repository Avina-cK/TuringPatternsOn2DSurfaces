"""
    assemble_elementK_surface!(Ke, cell_nodes, qr, ip, mapping)

Element stiffness using the custom SurfaceMapping directly.

Pipeline per quadrature point:
  1. compute_surface_mapping  ->  J, G, G⁻¹, √det(G)
  2. reference_shape_gradient -> ∇_hat_φ ∈ ℝ²
  3. map_shape_gradient        ->  ∇_x φ = J G⁻¹ ∇_hat_φ ∈ ℝ³
  4. (∇_x φᵢ ⋅ ∇_x φⱼ) * √det(G) * w  ->  Kₑ[i,j]
"""
function assemble_elementK_surface!(
    Ke::AbstractMatrix{T},
    cell_nodes::AbstractVector{<:Vec{3, Float64}},
    qr::QuadratureRule{RefTriangle},
    ip::LTI,
    mapping::SurfaceMapping
) where {T}
    n_basefuncs = getnbasefunctions(ip)
    fill!(Ke, zero(T))

    for (ξ, w) in zip(Ferrite.getpoints(qr), Ferrite.getweights(qr))

        # -- Step 1: compute all geometry at this quadrature point ---------#
        mv = compute_surface_mapping(cell_nodes, ξ)
        # area-weighted integration measure
        dΩ = get_detJ(mapping, mv) * w

        for i in 1:n_basefuncs

            # -- Step 2: reference gradient  ∇_hat_φᵢ ∈ ℝ² -----------------#
            ∇_hat_φᵢ = Ferrite.reference_shape_gradient(ip, ξ, i)

            # -- Step 3: physical surface gradient  ∇_x φᵢ ∈ ℝ³ ------------#
            # map_shape_gradient applies:  J G⁻¹ ∇_hat_φᵢ
            ∇φᵢ = map_shape_gradient(mapping, ∇_hat_φᵢ, mv)

            for j in 1:n_basefuncs
                ∇_hat_φⱼ = Ferrite.reference_shape_gradient(ip, ξ, j)
                ∇φⱼ = map_shape_gradient(mapping, ∇_hat_φⱼ, mv)

                # -- Step 4: weak Laplace–Beltrami contribution ------------#
                # ∇φᵢ ⋅ ∇φⱼ is the ℝ³ dot product of surface gradients
                Ke[i, j] += (∇φᵢ ⋅ ∇φⱼ) * dΩ
            end
        end
    end
    return Ke
end

"""
    assemble_globalK_surface(dh, qr, ip) → SparseMatrixCSC

Global stiffness assembly using SurfaceMapping directly.
Node coordinates are extracted as Vec{3} from the grid via CellIterator.
"""
function assemble_globalK_surface(
    dh::DofHandler,
    qr::QuadratureRule{RefTriangle},
    ip::LTI
)
    sp = init_sparsity_pattern(dh)  # uses DOF connectivity from mesh topology
    add_sparsity_entries!(sp, dh)   # fills entries for all shared-cell DOF pairs
    Kₕ = allocate_matrix(sp)        # -> SparseMatrixCSC{Float64,Int}
    
    mapping = SurfaceMapping()
    
    assembler = start_assemble(Kₕ)

    T = eltype(Ferrite.matrix_handle(assembler))

    n_basefuncs = getnbasefunctions(ip)
    Ke = zeros(T, n_basefuncs, n_basefuncs)
    

    for cell in CellIterator(dh)
        # getcoordinates returns Vector{Vec{3, Float64}} for a 3D grid
        # - the only geometric information we need
        cell_nodes =collect(getcoordinates(cell))

        assemble_elementK_surface!(Ke, cell_nodes, qr, ip, mapping)
        assemble!(assembler, celldofs(cell), Ke)
    end
    return Kₕ
end