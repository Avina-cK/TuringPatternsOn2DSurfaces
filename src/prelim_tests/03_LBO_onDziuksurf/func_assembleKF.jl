function assemble_elementKF!(
    Kₑ::Matrix, fₑ::Vector,
    cv::SurfaceCellValues,
    func_rhs::Function
    )

    no_bfs = Ferrite.getnbasefunctions(cv)
    fill!(Kₑ, 0.0)
    fill!(fₑ, 0.0)
    for q in 1:Ferrite.getnquadpoints(cv)
        dΩ_q = getdetJdV(cv, q) #√det(G)*w_q
        x = cv.mapping.x[q]
        fᵥ = func_rhs(x[1], x[2], x[3])
        for i in 1:no_bfs
            ∇φᵢ = Ferrite.shape_gradient(cv, q, i)
            φᵢ = Ferrite.shape_value(cv, q, i)
            fₑ[i] += fᵥ * φᵢ * dΩ_q
            for j in 1:no_bfs
                ∇φⱼ = Ferrite.shape_gradient(cv, q, j)
                Kₑ[i,j] += (∇φᵢ ⋅ ∇φⱼ) * dΩ_q
            end
        end
    end
    return Kₑ, fₑ
end

function assemble_globalKF(dh::Ferrite.DofHandler, cv::SurfaceCellValues, func_rhs::Function)
    ndofs = Ferrite.ndofs(dh)
    K = spzeros(ndofs, ndofs)
    f = zeros(ndofs)

    no_bfs = Ferrite.getnbasefunctions(cv)
    Kₑ = zeros(no_bfs, no_bfs)
    fₑ = zeros(no_bfs)

    for cell in Ferrite.CellIterator(dh)
        Ferrite.reinit!(cv, cell)

        # compute local element matrix/vector
        assemble_elementKF!(Kₑ, fₑ, cv, func_rhs)

        # global dof indices for this cell
        dofs = Ferrite.celldofs(cell)

        for i in 1:no_bfs
            gi = dofs[i]
            f[gi] += fₑ[i]
            for j in 1:no_bfs
                gj = dofs[j]
                K[gi, gj] += Kₑ[i, j]
            end
        end
    end

    return K, f
end
