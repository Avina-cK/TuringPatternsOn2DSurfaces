function assemble_elementKF!(
    Kₑ::Matrix, fₑ::Vector,
    cv::SurfaceCellValues,
    cell_coords, func_rhs::Function
    )

    no_bfs = getnbasefunctions(cv)
    fill!(Kₑ, 0.0)
    fill!(fₑ, 0.0)
    for q in 1:Ferrite.getnquadpoints(cv)
        dΩ_q = getdetJdV(cv, q) #√det(G)*w_q
        x = cv.mapping.x[q]
        fᵥ = func_rhs(x[1], x[2], x[3])
        for i in 1:no_bfs
            ∇φᵢ = shape_gradient(cv, q, i)
            φᵢ = shape_value(cv, q, i)
            fₑ[i] += fᵥ * φᵢ * dΩ_q
            for j in 1:no_bfs
                ∇φⱼ = shape_gradient(cv, q, j)
                Kₑ[i,j] += (∇φᵢ ⋅ ∇φⱼ) * dΩ_q
            end
        end
    end
    return Kₑ, fₑ
end

function assemble_globalKF(F::AbstractVector, dh::DofHandler, cv::SurfaceCellValues, func_rhs_g::Function)
    sp = init_sparsity_pattern(dh)
    add_sparsity_entries!(sp, dh)
    Kₕ = allocate_matrix(sp)
    fill!(F, 0.0)
    
    assembler = start_assemble(Kₕ, F)
    T = eltype(Ferrite.matrix_handle(assembler))
    no_bfs = getnbasefunctions(cv)

    Kₑ = zeros(T, no_bfs, no_bfs)
    fₑ = zeros(no_bfs)

    for cell in CellIterator(dh)
        Ferrite.reinit!(cv, cell)
        cell_coords = getcoordinates(cell)
        assemble_elementKF!(Kₑ, fₑ, cv, cell_coords, func_rhs_g)
        assemble!(assembler, celldofs(cell), Kₑ, fₑ)
    end
    return Kₕ, F
    
end