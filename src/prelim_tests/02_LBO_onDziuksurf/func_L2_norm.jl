function L2_norm(
        u_given::AbstractVector{Float64},
        dh::DofHandler,
        qr::QuadratureRule{RefTriangle},
        ip::Lagrange{RefTriangle}
    )
    norm_sq = 0.0
    n_basefuncs = getnbasefunctions(ip)
    for cell in CellIterator(dh)
        cell_nodes = collect(getcoordinates(cell))
        dofs = celldofs(cell)
        u_local = u_given[dofs]
        for (ξ, w) in zip(Ferrite.getpoints(qr), Ferrite.getweights(qr))
            mv = compute_surface_mapping(cell_nodes, ξ)
            dΩ = mv.detJ * w
            uh_q = sum(u_local[i] * Ferrite.reference_shape_value(ip, ξ, i) for i in 1:n_basefuncs)
            norm_sq += uh_q^2 * dΩ
        end
    end
    return sqrt(norm_sq)
end

function L2_error(
        u_true::AbstractVector{Float64},
        u_h::AbstractVector{Float64},
        dh::DofHandler,
        qr::QuadratureRule{RefTriangle},
        ip::Lagrange{RefTriangle}
    )
    norm_sq = 0.0
    n_basefuncs = getnbasefunctions(ip)
    for cell in CellIterator(dh)
        cell_nodes = collect(getcoordinates(cell))
        dofs = celldofs(cell)
        u_local_true = u_true[dofs]
        u_local_h = u_h[dofs]
        for (ξ, w) in zip(Ferrite.getpoints(qr), Ferrite.getweights(qr))
            mv = compute_surface_mapping(cell_nodes, ξ)
            dΩ = mv.detJ * w
            diff_q = sum((u_local_true[i] - u_local_h[i]) * Ferrite.reference_shape_value(ip, ξ, i) for i in 1:n_basefuncs)
            norm_sq += diff_q^2 * dΩ
        end
    end
    return sqrt(norm_sq)
end