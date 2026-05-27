function assemble_true_solution(
    u_true::Function,
    dh::DofHandler,
    ip::LTI
)
    u_true_v = zeros(ndofs(dh))
    for cell in CellIterator(dh)
        cell_nodes =collect(getcoordinates(cell)) # ∈ ℝ³
        dofs = celldofs(cell) #global DoF indices
        for (local_dof, ξ_node) in enumerate(Ferrite.reference_coordinates(ip))
            x = Φₜ(cell_nodes, ξ_node)  # map reference coordinate to surface coordinate
            global_dof = dofs[local_dof]
            u_true_v[global_dof] = u_true(x[1], x[2],x[3])
        end
    end
    return u_true_v
end