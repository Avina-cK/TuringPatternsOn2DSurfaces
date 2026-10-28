##
using LinearAlgebra
using Gmsh.gmsh

include("../func_floattostring.jl")
include("func_BloodCellSurface.jl")

d_giv=4.0
b_giv = 0.2
h_giv = 1.0

given_surface(x,y,z) = blood_cell(x,y,z;d=d_giv,b=b_giv, h=h_giv)
given_grad_surface(x,y,z) = grad_blood_cell(x,y,z;d=d_giv,b=b_giv, h=h_giv)

include("func_project_to_surf.jl")
include("func_surface_lloyd.jl")

##

gmsh.initialize(append!(["gmsh"], ARGS))

cd(@__DIR__)
init_ref=0
final_ref = 7

mesh_name="Blood_cell_mesh_$(fts(d_giv, 2))_$(fts(b_giv, 2))_$(fts(h_giv, 2))"

filename="../../blood_cell_meshes/$(mesh_name)_$(init_ref).msh"
gmsh.open(filename)

for r in init_ref+1:final_ref
    gmsh.model.mesh.refine()
    node_tags, node_coords = gmsh.model.mesh.getNodes()
    no_nodes = length(node_tags)
    node_coords=reshape(node_coords, 3, no_nodes)

    newnodes=[]
    for i in 1:no_nodes
        node_i = [node_coords[1, i], node_coords[2, i], node_coords[3, i]]
        new_node_i = [0.,0.,0.]
        new_node_i .= project_to_surf(node_i[1], node_i[2], node_i[3])
        push!(newnodes, new_node_i)
    end

    for i in 1:no_nodes
        gmsh.model.mesh.setNode(node_tags[i], newnodes[i], [])
    end
   
    triangles = gmsh.model.mesh.getElements(2)[3][1]
    triangles = reshape(triangles, 3, :)
    vertices = gmsh.model.mesh.getNodes()[2]
    vertices = reshape(node_coords, 3, :)
    lloyded_vertices = surface_lloyd(vertices, triangles)
            
    if r<=5
        n_v = length(lloyded_vertices)
        for i in 1:n_v
            if abs(lloyded_vertices[i])<1e-18
                lloyded_vertices[i]=0.0
            end
        end
    end
    lloyded_vertices=reshape(lloyded_vertices, 3, no_nodes)

    for i in 1:no_nodes
        gmsh.model.mesh.setNode(node_tags[i], lloyded_vertices[1:3,i], [])
    end

    cd(@__DIR__)
    gmsh.write("../../blood_cell_meshes/$(mesh_name)_$(r).msh")

end

gmsh.finalize()
