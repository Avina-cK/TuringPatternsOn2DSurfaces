"""
assemble_elementM(Mₑ, cv::SurfaceCellValues)
Assemble mass matrix M's element-contribution
    Mₑ[i,j] += (φᵢ ⋅ φⱼ) * dΩ_q
where φᵢ, φⱼ are the shape values 
"""
function assemble_elementM(Mₑ, cv::SurfaceCellValues)
    no_bfs = Ferrite.getnbasefunctions(cv)
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
    return Mₑ
end