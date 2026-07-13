function assemble_globalF(cv::SurfaceCellValues, dh::Ferrite.DofHandler, func_rhs::Function)
    ndofs = Ferrite.ndofs(dh)
    f = zeros(ndofs)

    no_bfs = Ferrite.getnbasefunctions(cv)
    fₑ = zeros(no_bfs)

    for cell in Ferrite.CellIterator(dh)
        Ferrite.reinit!(cv, cell)

        # compute local element vector
        # ∫ (f(x) ⋅ φᵢ) * dΩ_q
        fill!(fₑ, 0.0)
        for q in 1:Ferrite.getnquadpoints(cv)
            dΩ_q = getdetJdV(cv, q) #√det(G)*w_q
            x = cv.mapping.x[q]
            fᵥ = func_rhs(x[1], x[2], x[3])
            for i in 1:no_bfs
                φᵢ = Ferrite.shape_value(cv, q, i)
                fₑ[i] += fᵥ * φᵢ * dΩ_q
            end
        end

        # global dof indices for this cell
        dofs = Ferrite.celldofs(cell)
        for i in 1:no_bfs
            gi = dofs[i]
            f[gi] += fₑ[i]
        end
    end

    return f
end
