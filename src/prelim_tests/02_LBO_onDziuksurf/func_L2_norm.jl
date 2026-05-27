function L2_norm(
    u_given :: AbstractVector{Float64},
    dh :: DofHandler,
    qr :: QuadratureRule{RefTriangle},
    ip :: LTI
)
    mapping = SurfaceMapping()
    norm_sq = 0.0

    for cell in CellIterator(dh)
        cell_nodes = collect(getcoordinates(cell))     # Vector{Vec{3,Float64}}
        dofs = celldofs(cell)
        u_local = u_given[dofs]                 # local DOF values

        for (ξ, w) in zip(Ferrite.getpoints(qr), Ferrite.getweights(qr))

            # -- Geometry ---------------------------------------#
            mv = compute_surface_mapping(cell_nodes, ξ)
            dΩ = get_detJ(mapping, mv) * w  # √det(G)*wq

            # -- Interpolated value at quadrature point ---------------#
            # u_h(ξ) = ∑ᵢ uᵢ · φᵢ_hat  (scalar, identity mapping)
            uh_q = sum(
                u_local[i] * map_shape_value(mapping, Ferrite.reference_shape_value(ip, ξ, i), mv)
                for i in 1:getnbasefunctions(ip)
            )

            norm_sq += uh_q^2 * dΩ
        end
    end

    return sqrt(norm_sq)
end

function L₂_error(
    u_true :: AbstractVector{Float64},
    u_h ::AbstractVector{Float64},
    dh :: DofHandler,
    qr :: QuadratureRule{RefTriangle},
    ip :: LTI
)
    mapping = SurfaceMapping()
    norm_sq = 0.0

    for cell in CellIterator(dh)
        cell_nodes =collect(getcoordinates(cell))     # Vector{Vec{3,Float64}}
        dofs = celldofs(cell)
        u_local_true = u_true[dofs] # local DOF values
        u_local_h = u_h[dofs]
        for (ξ, w) in zip(Ferrite.getpoints(qr), Ferrite.getweights(qr))

            # -- Geometry ---------------------------------------#
            mv = compute_surface_mapping(cell_nodes, ξ)
            dΩ = get_detJ(mapping, mv) * w  # √det(G)*wq

            # -- Interpolated value at quadrature point ---------------#
            # u_h(ξ) = ∑ᵢ uᵢ · φᵢ_hat  (scalar, identity mapping)
            uth_q = sum(
                u_local_true[i] * map_shape_value(mapping, Ferrite.reference_shape_value(ip, ξ, i), mv)
                for i in 1:getnbasefunctions(ip)
            )
            uhh_q = sum(
                u_local_h[i] * map_shape_value(mapping, Ferrite.reference_shape_value(ip, ξ, i), mv)
                for i in 1:getnbasefunctions(ip)
            )

            norm_sq += (uth_q-uhh_q)^2 * dΩ
        end
    end

    return sqrt(norm_sq)
end