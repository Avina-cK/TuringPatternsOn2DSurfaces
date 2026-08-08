# To use only inside an open Gmsh instance
using LinearAlgebra

"""
    mesh_size() -> max_len::Float64

    measures the maximum edge length of the mesh.
    To use only inside an open Gmsh instance.
"""
function mesh_size()
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    etypes, _, enodes = gmsh.model.mesh.getElements(2)

    tri_nodes = nothing
    for (i, t) in enumerate(etypes)
        if t == 2
            tri_nodes = reshape(enodes[i], 3, :)
            break
        end
    end
    tri_nodes === nothing && error("No triangle elements found")

    n = length(node_coords) ÷ 3
    V = reshape(node_coords, 3, n)
    tag2idx = Dict(tag => i for (i, tag) in enumerate(node_tags))

    max_len = 0.0
    for t in eachcol(tri_nodes)
        a, b, c = t
        pa = V[:, tag2idx[a]]
        pb = V[:, tag2idx[b]]
        pc = V[:, tag2idx[c]]

        max_len = max(max_len,
                      norm(pb - pa),
                      norm(pc - pb),
                      norm(pa - pc))
    end
    return max_len
end
