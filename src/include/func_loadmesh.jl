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

function load_mesh(filename::String)
        # -- Open the mesh ----#
        gmsh.initialize()
        gmsh.open(filename)

        # --- Extract node coordinates and tags---#
        node_tags, node_coords, _ = gmsh.model.mesh.getNodes()

        # Map gmsh node-tag -> contiguous 1-based index
        id_map = Dict(id => i for (i, id) in enumerate(node_tags))

        nodes = [Ferrite.Node(Vec{3,Float64}((node_coords[3i-2], node_coords[3i-1], node_coords[3i])))
            for i in 1:length(node_tags)]

        # --- Extract element coordinates and tags---#
        etypes, _, ntags = gmsh.model.mesh.getElements(2)
        min_edge_length, max_edge_length = edge_length_stats2(node_tags, node_coords, etypes, ntags)

        tri_idx = findfirst(==(2), etypes)
        tri_idx === nothing && error("No triangle elements found.")

        tri_ntags = ntags[tri_idx]
        n_cells = length(tri_ntags) ÷ 3

        cells = [
                Ferrite.Triangle((
                id_map[tri_ntags[3i-2]],
                id_map[tri_ntags[3i-1]],
                id_map[tri_ntags[3i]]
                ))
                for i in 1:n_cells
        ]
        gmsh.finalize()

        # -- Build the Ferrite grid ---#
        Ωₕ = Ferrite.Grid(cells, nodes)
        return Ωₕ, max_edge_length
end