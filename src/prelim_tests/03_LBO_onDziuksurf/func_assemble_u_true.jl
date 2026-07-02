function assemble_true_solution(
    u_true::Function,
    dh::DofHandler,
    ip::Lagrange{RefTriangle, 1}
)
    u_true_v = zeros(ndofs(dh))
    n_basefuncs = Ferrite.getnbasefunctions(ip)
    ue = zeros(n_basefuncs)

    for cell in Ferrite.CellIterator(dh)
        cell_nodes = collect(Ferrite.getcoordinates(cell))
        fill!(ue, 0.0)
        for (i, ξ_node) in enumerate(Ferrite.reference_coordinates(ip))
            x = Φₜ(cell_nodes, ξ_node)
            ue[i] = u_true(x[1], x[2], x[3])
        end
        dofs = Ferrite.celldofs(cell)
        for (i, dof) in enumerate(dofs)
            u_true_v[dof] = ue[i]
        end
    end
    return u_true_v
end


function assemble_elementtrue!(
    uₑ::Vector,
    cv::SurfaceCellValues,
    u_true::Function
    )

    no_bfs = Ferrite.getnbasefunctions(cv)
    fill!(uₑ, 0.0)
    for q in 1:Ferrite.getnquadpoints(cv)
        dΩ_q = getdetJdV(cv, q) #√det(G)*w_q
        x = cv.mapping.x[q]
        uᵥ = u_true(x[1], x[2], x[3])
        for i in 1:no_bfs

            uₑ[i] += uᵥ
        end
    end
    return uₑ
end

function assemble_globalU(dh::Ferrite.DofHandler, cv::SurfaceCellValues, u_true::Function)
    ndofs = Ferrite.ndofs(dh)
    U = zeros(ndofs)

    no_bfs = Ferrite.getnbasefunctions(cv)
    uₑ = zeros(no_bfs)

    for cell in Ferrite.CellIterator(dh)
        Ferrite.reinit!(cv, cell)
        assemble_elementtrue!(uₑ, cv, u_true)
        dofs = Ferrite.celldofs(cell)
        for (i, dof) in enumerate(dofs)
            U[dof] += uₑ[i]
        end
    end

    return U
end
