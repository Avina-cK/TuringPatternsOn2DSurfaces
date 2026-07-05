using Gmsh.gmsh, Ferrite, FerriteGmsh

function extract_triangles(etypes, enodes)
    for (i, t) in enumerate(etypes)
        if t == 2  # 3-node triangle
            return reshape(enodes[i], 3, :)
        end
    end
    error("No triangle elements found")
end

function edge_length_stats2(node_tags, coords, etypes, enodes)
    n = length(coords) ÷ 3
    V = reshape(coords, 3, n)
    tag2idx = Dict(tag => i for (i, tag) in enumerate(node_tags))
    T = extract_triangles(etypes, enodes)
    min_len = Inf
    max_len = 0.0
    for t in eachcol(T)
        a, b, c = t
        pa = V[:, tag2idx[a]]
        pb = V[:, tag2idx[b]]
        pc = V[:, tag2idx[c]]
        max_len = max(max_len,
                      norm(pb - pa),
                      norm(pc - pb),
                      norm(pa - pc))
        min_len = min(min_len,
                      norm(pb - pa),
                      norm(pc - pb),
                      norm(pa - pc))
    end
    return min_len, max_len
end

function signed_area(pa, pb, pc)
    x1, y1 = pa[1], pa[2]
    x2, y2 = pb[1], pb[2]
    x3, y3 = pc[1], pc[2]

    return 0.5 * (
        (x2 - x1)*(y3 - y1) -
        (x3 - x1)*(y2 - y1)
    )
end

function load_mesh(filename::String)
    gmsh.isInitialized()!=0 && gmsh.finalize()
        # -- Open the mesh ----#
        gmsh.initialize()
    try
        gmsh.open(filename)

        # --- Extract node coordinates and tags---#
        node_tags, node_coords, _ = gmsh.model.mesh.getNodes(2)

        n_nodes = length(node_tags)

        # Gmsh node tag -> contiguous Ferrite index
        id_map = Dict(tag => i for (i, tag) in enumerate(node_tags))

        nodes = [
            Ferrite.Node(Vec{3,Float64}((
                node_coords[3i-2],
                node_coords[3i-1],
                node_coords[3i]
            )))
            for i in 1:n_nodes
        ]


        # --- Extract element coordinates and tags---#
        etypes, _, enodes = gmsh.model.mesh.getElements(2)
        min_edge_length, max_edge_length = edge_length_stats2(node_tags, node_coords, etypes, enodes)

        tri_idx = findfirst(==(2), etypes)
        tri_idx === nothing && error("No triangle elements found.")
        tri_ntags = enodes[tri_idx]
        n_cells   = length(tri_ntags) ÷ 3
        @assert length(tri_ntags) % 3 == 0 "Triangle node tag count not divisible by 3"
        
        cells = Ferrite.Triangle[]
        cells = [
            Ferrite.Triangle((
                id_map[tri_ntags[3i-2]],
                id_map[tri_ntags[3i-1]],
                id_map[tri_ntags[3i]]
            ))
            for i in 1:n_cells
        ]

        # Validate all cell node indices are in range
        n_nodes = length(nodes)
        for (ci, cell) in enumerate(cells)
            for nidx in cell.nodes
                1 ≤ nidx ≤ n_nodes ||
                    error("Cell $ci references out-of-range node index $nidx")
            end
        end
        # -- Build the Ferrite grid ---#
        Ωₕ = Ferrite.Grid(cells, nodes)
        return Ωₕ, max_edge_length
    finally
        gmsh.isInitialized()!=0 && gmsh.finalize()
    end
end