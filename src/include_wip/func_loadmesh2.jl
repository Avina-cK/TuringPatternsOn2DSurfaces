include("func_mesh_size.jl")
function load_mesh(refinements, dir_to_meshes="../Dziuk_surf_meshes")
    gmsh.initialize()
    gmsh.open("$(dir_to_meshes)/lowresmesh_$(refinements).msh")
    nodes = tonodes()
    elements, _ = toelements(2)
    Ω = Grid(elements, nodes)
    max_edge_length = mesh_size()
    gmsh.finalize()
    return Ω, max_edge_length
end