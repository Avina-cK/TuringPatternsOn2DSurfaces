using LinearAlgebra

"""
    check_normal_orientation(dh, cv::SurfaceCellValues) -> Vector{Int}

For every cell, reinit `cv`, take the (constant, since the triangle map is
affine) surface normal, and check that it points outward — i.e. that
`n̂ · x̂_centroid > 0`, where `x̂_centroid` is the unit vector from the origin
to the cell's centroid. 

Returns the list of cell indices whose normal fails the check (candidates
for flipped orientation / inconsistent winding).
"""
function check_normal_orientation(dh::DofHandler, cv::SurfaceCellValues)
    grid = dh.grid
    flipped = Int[]
    min_dot = Inf
    max_dot = -Inf

    for (cellid, cell) in enumerate(CellIterator(dh))
        reinit!(cv, cell)

        # Physical centroid: average of the 3 node coordinates.
        cell_nodes = Ferrite.getcoordinates(cell)
        centroid = (cell_nodes[1] + cell_nodes[2] + cell_nodes[3]) / 3
        x̂_centroid = centroid / norm(centroid)

        # For a linear (affine) triangle map, J and hence the normal is
        # constant over the cell, so any quadrature point works; use q=1.
        n = getnormal(cv, 1)

        d = n ⋅ x̂_centroid
        min_dot = min(min_dot, d)
        max_dot = max(max_dot, d)

        if d <= 0
            push!(flipped, cellid)
        end
    end
    return flipped
end

"""
    flip_cells!(grid::Ferrite.Grid, cellids::Vector{Int})

Reverse the local node ordering of the given cells to flip their outward
normal. Swapping any two local nodes reverses the winding direction
(here nodes 2 and 3), which is enough to fix the cross-product sign in
`surface_normal(J_Φ)` without changing which physical triangle is defined.
"""
function flip_cells!(grid::Ferrite.Grid, cellids::Vector{Int})
    for cellid in cellids
        cell = grid.cells[cellid]
        n1, n2, n3 = cell.nodes
        grid.cells[cellid] = typeof(cell)((n1, n3, n2))  # swap 2 ↔ 3
    end
    return grid
end


"""
    fix_orientation!(grid::Ferrite.Grid) -> Vector{Int}

Make triangle orientation globally consistent via edge-adjacency flood-fill,
then fix the overall sign using the divergence-theorem volume check.
Returns the indices of cells that were flipped.
"""
function fix_orientation!(grid::Ferrite.Grid)
    ncells = length(grid.cells)

    edge_map = Dict{Tuple{Int,Int}, Vector{Tuple{Int,Int}}}()
    for (cellid, cell) in enumerate(grid.cells)
        n1, n2, n3 = cell.nodes
        for (a, b) in ((n1, n2), (n2, n3), (n3, n1))
            key = a < b ? (a, b) : (b, a)
            dir = a < b ? 1 : -1
            push!(get!(edge_map, key, Tuple{Int,Int}[]), (cellid, dir))
        end
    end

    adjacency = [Tuple{Int,Bool}[] for _ in 1:ncells]  # (neighbor, same_dir)
    for (_, lst) in edge_map
        if length(lst) == 2
            (c1, d1), (c2, d2) = lst
            same_dir = d1 == d2
            push!(adjacency[c1], (c2, same_dir))
            push!(adjacency[c2], (c1, same_dir))
        else
            @warn "Edge shared by $(length(lst)) cells (expected 2) — not a closed manifold?"
        end
    end

    visited = falses(ncells)
    to_flip = falses(ncells)
    for start in 1:ncells
        visited[start] && continue
        visited[start] = true
        queue = [start]
        while !isempty(queue)
            c = popfirst!(queue)
            for (nb, same_dir) in adjacency[c]
                want_flip = same_dir ? !to_flip[c] : to_flip[c]
                if !visited[nb]
                    visited[nb] = true
                    to_flip[nb] = want_flip
                    push!(queue, nb)
                elseif to_flip[nb] != want_flip
                    @warn "Inconsistent (non-orientable?) mesh at cells $c, $nb"
                end
            end
        end
    end

    flipped_ids = findall(to_flip)
    flip_cells!(grid, flipped_ids)

    V6 = sum(grid.cells) do cell
        n1, n2, n3 = cell.nodes
        q1, q2, q3 = grid.nodes[n1].x, grid.nodes[n2].x, grid.nodes[n3].x
        q1 ⋅ (q2 × q3)
    end
    if V6 < 0
        all_ids = collect(1:ncells)
        flip_cells!(grid, all_ids)
        flipped_ids = symdiff(flipped_ids, all_ids)
        @info "Entire mesh was inside-out — flipped all $ncells cells"
    end

    return flipped_ids
end