# Source: ferrite-fem.github.io/Ferrite.jl/stable/tutorials/reactive_surface/
using Ferrite, FerriteGmsh
function create_embedded_sphere(refinements; tosave=false)
    gmsh.initialize()

    # Add a unit sphere in 3D space
    gmsh.model.occ.addSphere(0.0, 0.0, 0.0, 1.0)
    gmsh.model.occ.synchronize()

    # Generate nodes and surface elements only, hence we need to pass 2 into generate
    gmsh.model.mesh.generate(2)

    # To get good solution quality refine the elements several times
    for _ in 1:refinements
        gmsh.model.mesh.refine()
    end

    if tosave
        cd(@__DIR__)
        gmsh.write("../sphere_meshes/embedded_sphere_$(refinements).msh")
    end
    # Now we create a Ferrite grid out of it. Note that we also call toelements
    # with our surface element dimension to obtain these.
    nodes = tonodes()
    elements, _ = toelements(2)
    gmsh.finalize()
    return Grid(elements, nodes)
end

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


function load_sphere(refinements, dir_to_spheres)
    gmsh.initialize()
    gmsh.open("$(dir_to_spheres)/embedded_sphere_$(refinements).msh")
    nodes = tonodes()
    elements, _ = toelements(2)
    Ω = Grid(elements, nodes)
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    etypes, _, enodes = gmsh.model.mesh.getElements(2)
    _, max_edge_length = edge_length_stats2(node_tags, node_coords, etypes, enodes)

    gmsh.finalize()
    return Ω, max_edge_length
end