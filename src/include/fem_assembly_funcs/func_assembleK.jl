"""
assemble_globalK(cv::SurfaceCellValues, dh::Ferrite.DofHandler; k::Float64=1.0)
Assemble global stiffness matrix K that corresponds to the weak form
    ∫ k(∇φᵢ ⋅ ∇φⱼ) dΩₕ

with element-contribution
    Kₑ[i,j] += k * (∇φᵢ ⋅ ∇φⱼ) * dΩ_q
where ∇φᵢ, ∇φⱼ are the shape gradients 
"""
function assemble_globalK(cv::SurfaceCellValues, dh::Ferrite.DofHandler; k::Float64=1.0)
    ndofs = Ferrite.ndofs(dh)
    K = spzeros(ndofs, ndofs)
    no_bfs = Ferrite.getnbasefunctions(cv)
    Kₑ = zeros(no_bfs, no_bfs)
    for cell in Ferrite.CellIterator(dh)
        Ferrite.reinit!(cv, cell)  
        
        # Assemble elemental stiffness matrix Kₑ
        # ∫ k(∇φᵢ ⋅ ∇φⱼ) * dΩ_q
        fill!(Kₑ, 0.0)
        for q in 1:Ferrite.getnquadpoints(cv)
            dΩ_q = getdetJdV(cv, q)
            for i in 1:no_bfs
                ∇φᵢ = Ferrite.shape_gradient(cv, q, i)
                for j in 1:no_bfs
                    ∇φⱼ = Ferrite.shape_gradient(cv, q, j)
                    Kₑ[i,j] += k * (∇φᵢ ⋅ ∇φⱼ) * dΩ_q
                end
            end
        end

        dofs = Ferrite.celldofs(cell)
        for i in 1:no_bfs
            gi = dofs[i]
            for j in 1:no_bfs
                gj = dofs[j]
                K[gi, gj] += Kₑ[i, j]
            end
        end
    end
    return K
end