"""
Assemble global mass matrix M
"""
function assemble_globalM(cv::SurfaceCellValues, dh::Ferrite.DofHandler)
    ndofs = Ferrite.ndofs(dh)
    M = spzeros(ndofs, ndofs)
    no_bfs = Ferrite.getnbasefunctions(cv)
    Mₑ = zeros(no_bfs, no_bfs)
    for cell in Ferrite.CellIterator(dh)
        Ferrite.reinit!(cv, cell)  
        
        # Assemble elemental mass matrix Mₑ
        # ∫ (φᵢ ⋅ φⱼ) * dΩ_q
        fill!(Mₑ, 0.0)
        for q in 1:Ferrite.getnquadpoints(cv)
            dΩ_q = getdetJdV(cv, q)
            for i in 1:no_bfs
                φᵢ = Ferrite.shape_value(cv, q, i)
                for j in 1:no_bfs
                    φⱼ = Ferrite.shape_value(cv, q, j)
                    Mₑ[i,j] += (φᵢ ⋅ φⱼ) * dΩ_q
                end
            end
        end

        dofs = Ferrite.celldofs(cell)
        for i in 1:no_bfs
            gi = dofs[i]
            for j in 1:no_bfs
                gj = dofs[j]
                M[gi, gj] += Mₑ[i, j]
            end
        end
    end
    return M
end