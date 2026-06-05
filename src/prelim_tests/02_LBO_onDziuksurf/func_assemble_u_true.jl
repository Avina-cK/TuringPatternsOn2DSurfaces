function assemble_true_solution(
    u_true::Function,
    dh::DofHandler,
    ip::Lagrange{RefTriangle, 1}
)
    u_true_v = zeros(ndofs(dh))
    n_basefuncs = getnbasefunctions(ip)
    ue = zeros(n_basefuncs)

    for cell in CellIterator(dh)
        cell_nodes = collect(getcoordinates(cell))
        fill!(ue, 0.0)
        for (i, ξ_node) in enumerate(Ferrite.reference_coordinates(ip))
            x = Φₜ(cell_nodes, ξ_node)
            ue[i] = u_true(x[1], x[2], x[3])
        end
        dofs = celldofs(cell)
        for (i, dof) in enumerate(dofs)
            u_true_v[dof] = ue[i]
        end
    end
    return u_true_v
end