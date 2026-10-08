##
using LinearAlgebra
include("func_DzuikSurface.jl")

given_surface = Dziuk_surface
given_grad_surface = gradDziuk_surface

include("func_project_to_surf.jl")
include("func_surface_lloyd.jl")

##
using Gmsh.gmsh

gmsh.initialize(append!(["gmsh"], ARGS))

cd(@__DIR__)
init_ref=0
final_ref = 8


filename="../../Dziuk_surf_meshes/Dzuik_mesh_$(init_ref).msh"
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
    gmsh.write("../../Dziuk_surf_meshes/Dzuik_mesh_$(r).msh")

end

gmsh.finalize()
