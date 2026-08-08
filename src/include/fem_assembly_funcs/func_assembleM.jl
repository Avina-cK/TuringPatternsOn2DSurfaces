include("func_assembleMe.jl")
"""
assemble_globalM(cv::SurfaceCellValues, dh::Ferrite.DofHandler; k::Float64=1.0)
Assemble global mass matrix M that corresponds to the weak form
    ∫ (φᵢ ⋅ φⱼ) dΩₕ

with element-contribution
    Mₑ[i,j] += (φᵢ ⋅ φⱼ) * dΩ_q
where φᵢ, φⱼ are the shape values 
"""
function assemble_globalM(cv::SurfaceCellValues, dh::Ferrite.DofHandler)
    ndofs = Ferrite.ndofs(dh)
    M = spzeros(ndofs, ndofs)
    no_bfs = Ferrite.getnbasefunctions(cv)
    Mₑ = zeros(no_bfs, no_bfs)
    for cell in Ferrite.CellIterator(dh)
        Ferrite.reinit!(cv, cell)
        # Assemble elemental mass matrix Mₑ
        Mₑ = assemble_elementM(Mₑ, cv)
        # Add Mₑ to the global matrix
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