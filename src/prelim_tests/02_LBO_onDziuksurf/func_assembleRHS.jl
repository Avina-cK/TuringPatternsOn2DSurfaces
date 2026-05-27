function assemble_rhs_forces!(
    F:: AbstractVector,
    dh:: DofHandler,
    qr:: QuadratureRule{RefTriangle},
    ip:: LTI,
    func_rhs_assemble:: Function
)
    fill!(F, 0.0)

    mapping= SurfaceMapping()
    n_basefuncs = getnbasefunctions(ip)
    fe = zeros(n_basefuncs)

    for cell in CellIterator(dh)
        cell_nodes =collect(getcoordinates(cell)) # ∈ℝ³
        fill!(fe, 0.0)

        for (ξ, w) in zip(Ferrite.getpoints(qr), Ferrite.getweights(qr))

            # -- Geometry at this quadrature point -------------------------
            mv = compute_surface_mapping(cell_nodes, ξ)  # J, G, invG, detJ
            dΩ = get_detJ(mapping, mv) * w  # √det(G) · wq

            # -- Physical coordinate x = Φₜ(ξ) ∈ ℝ³ ----------------------#
            x = Φₜ(cell_nodes, ξ)

            # -- Evaluate RHS forcing function at the physical point ------#
            fᵥ = func_rhs_assemble(x[1], x[2],x[3])

            for i in 1:n_basefuncs
                # Scalar shape value = identity mapping, same in ref and physical
                φᵢ_hat = Ferrite.reference_shape_value(ip, ξ, i)
                φᵢ = map_shape_value(mapping, φᵢ_hat, mv)

                fe[i] += fᵥ * φᵢ * dΩ
            end
        end

        # Scatter local fe into global F at the correct DOF indices
        # assemble!(g, dofs, ge)->vector form from assembly devdocs
        assemble!(F, celldofs(cell), fe)
    end

    return F
end