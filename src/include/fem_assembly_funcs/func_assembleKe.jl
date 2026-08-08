"""
assemble_elementK(Kₑ, cv::SurfaceCellValues; k::Float64=1.0)
Assemble stiffness matrix K's element-contribution
    Kₑ[i,j] += k * (∇φᵢ ⋅ ∇φⱼ) * dΩ_q
where ∇φᵢ, ∇φⱼ are the shape gradients 
"""
function assemble_elementK(Kₑ, cv::SurfaceCellValues; k::Float64=1.0)
    no_bfs = Ferrite.getnbasefunctions(cv)
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
    return Kₑ
end