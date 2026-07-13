function assemble_U(cv::SurfaceCellValues, dh::Ferrite.DofHandler, u_true::Function)
    ndofs = Ferrite.ndofs(dh)
    u = fill(NaN, ndofs)   # NaN sentinel: "not yet filled"
    no_bfs = Ferrite.getnbasefunctions(cv)

    for cell in Ferrite.CellIterator(dh)
        Ferrite.reinit!(cv, cell)

        cell_nodes = Ferrite.getcoordinates(cell)
        dofs = Ferrite.celldofs(cell)

        for i in 1:no_bfs
            dof = dofs[i]
            isnan(u[dof]) || continue

            x = cell_nodes[i]
            u[dof] = u_true(x[1], x[2], x[3])
        end
    end

    return u
end

surface_mean(u::Vector, c::Vector) = (c ⋅ u) / sum(c)