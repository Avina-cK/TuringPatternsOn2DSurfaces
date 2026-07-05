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
    K = zeros(ndofs, ndofs)
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


function assemble_integral_uh(
    cv::SurfaceCellValues,
    uₕ::Vector
    )
    no_bfs = Ferrite.getnbasefunctions(cv)
    uₑ = zeros(no_bfs)
    fill!(uₑ, 0.0)
    for q in 1:Ferrite.getnquadpoints(cv)
        dΩ_q = getdetJdV(cv, q) # √det(G) * w_q
        x = cv.mapping.x[q]

        u_q = 0.0
        for j in 1:no_bfs
            φⱼ = Ferrite.shape_value(cv, q, j)
            u_q += uₕ[j] * φⱼ
        end

        for i in 1:no_bfs
            φᵢ = Ferrite.shape_value(cv, q, i)
            uₑ[i] += u_q * φᵢ * dΩ_q
        end
    end
    return uₑ
end

function assemble_integral_global_uh(dh::Ferrite.DofHandler, cv::SurfaceCellValues, uₕ::Vector)
    ndofs = Ferrite.ndofs(dh)
    U = zeros(ndofs)

    no_bfs = Ferrite.getnbasefunctions(cv)

    uₑ = zeros(no_bfs)

    for cell in Ferrite.CellIterator(dh)
        Ferrite.reinit!(cv, cell)

        # compute local element matrix/vector
        uₑ = assemble_integral_uh(cv, uₕ[celldofs(cell)])

        # global dof indices for this cell
        dofs = Ferrite.celldofs(cell)

        for i in 1:no_bfs
            gi = dofs[i]
            U[gi] += uₑ[i]
        end
    end

    return U
end
