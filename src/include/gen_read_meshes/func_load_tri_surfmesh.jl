using Gmsh.gmsh
using FerriteGmsh

include("func_mesh_size.jl")

"""
    load_tri_surfmesh(mesh_file) -> Ω::Ferrite.Grid, max_edge_length::Float64

Read a triangular surface mesh into a Ferrite.Grid format and calculate its mesh size.
"""
function load_tri_surfmesh(mesh_file)

    gmsh.isInitialized()!=0 && gmsh.finalize()
    gmsh.initialize()
    
    gmsh.open(mesh_file)
    node_order, node_coords = gmsh.model.mesh.getNodes()
    node_order = Int.(node_order)
    nnodes = length(node_order)
    node_coords = reshape(node_coords,3, nnodes)

    node_coords_vec_list::Vector{Node{3, Float64}} = []
    for i in 1:nnodes
        node_i = Node((node_coords[1, i], node_coords[2, i], node_coords[3,i]))
        push!(node_coords_vec_list, node_i)
    end
    node_coords_vec_list = node_coords_vec_list[node_order]
    elements, _ = FerriteGmsh.toelements(2)
    
    Ω = Grid(elements, node_coords_vec_list)
    max_edge_length = mesh_size()
    gmsh.finalize()

    return Ω, max_edge_length
end
