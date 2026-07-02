"""
    assemble_constraint_C!(C, dh, qr, ip) -> C

Assemble the zero-mean constraint row vector C ∈ ℝ^(1 × n_dofs):

    C[j] = ∫_Ω φ_j dΩ  for each global DOF j

Built by looping over cells and quadrature points, exactly mirroring
the RHS assembly, but integrating shape values only (no forcing function).
"""
function assemble_constraint_C!(
    C:: AbstractVector,     # length n_dofs: stores C as a dense row
    dh:: DofHandler,
    qr:: QuadratureRule{RefTriangle},
    ip:: Lagrange{RefTriangle, 1}
)
    fill!(C, 0.0)

    n_basefuncs = getnbasefunctions(ip)
    ce = zeros(n_basefuncs)   # local constraint contributions

    for cell in CellIterator(dh)
        cell_nodes =collect(getcoordinates(cell))  #Vector{Vec{3,Float64}}
        fill!(ce, 0.0)

        for (ξ, w) in zip(Ferrite.getpoints(qr), Ferrite.getweights(qr))
            mv = compute_surface_mapping(cell_nodes, ξ)
            dΩ = mv.detJ * w
            for i in 1:n_basefuncs
                φᵢ = Ferrite.reference_shape_value(ip, ξ, i)    #φᵢ = φᵢ_hat
                ce[i] += φᵢ * dΩ
            end
        end

        # Scatter local ce into global C at the correct DOF positions
        dofs = celldofs(cell)
        for (i, dof) in enumerate(dofs)
            C[dof] += ce[i]
        end
    end

    return C
end


function assemble_ele_constraint_C!(
    cₑ::Vector,
    cv::SurfaceCellValues
    )
    no_bfs = Ferrite.getnbasefunctions(cv)
    fill!(cₑ, 0.0)
    for q in 1:Ferrite.getnquadpoints(cv)
        dΩ_q = getdetJdV(cv, q) #√det(G)*w_q
        for i in 1:no_bfs
            φᵢ = Ferrite.shape_value(cv, q, i)
            cₑ[i] += φᵢ * dΩ_q
        end
    end
    return cₑ
end

function assemble_globalC(dh::Ferrite.DofHandler, cv::SurfaceCellValues)
    ndofs = Ferrite.ndofs(dh)
    C = zeros(ndofs)

    no_bfs = Ferrite.getnbasefunctions(cv)
    cₑ = zeros(no_bfs)

    for cell in Ferrite.CellIterator(dh)
        Ferrite.reinit!(cv, cell)
        assemble_ele_constraint_C!(cₑ, cv)
        dofs = Ferrite.celldofs(cell)
        for (i, dof) in enumerate(dofs)
            C[dof] += cₑ[i]
        end
    end

    return C
end