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



## To check
function assemble_globalKF_sparse(dh::Ferrite.DofHandler, cv::SurfaceCellValues, func_rhs::Function)
    ndofs = Ferrite.ndofs(dh)
    f = zeros(ndofs)

    no_bfs = Ferrite.getnbasefunctions(cv)
    Kₑ = zeros(no_bfs, no_bfs)
    fₑ = zeros(no_bfs)

    I = Int[]; J = Int[]; V = Float64[]
    sizehint!(I, Ferrite.getncells(dh) * no_bfs^2)
    sizehint!(J, Ferrite.getncells(dh) * no_bfs^2)
    sizehint!(V, Ferrite.getncells(dh) * no_bfs^2)

    for cell in Ferrite.CellIterator(dh)
        Ferrite.reinit!(cv, cell)
        assemble_elementKF!(Kₑ, fₑ, cv, func_rhs)
        dofs = Ferrite.celldofs(cell)

        for i in 1:no_bfs
            f[dofs[i]] += fₑ[i]
            for j in 1:no_bfs
                push!(I, dofs[i])
                push!(J, dofs[j])
                push!(V, Kₑ[i, j])
            end
        end
    end

    K = sparse(I, J, V, ndofs, ndofs)
    return K, f
end