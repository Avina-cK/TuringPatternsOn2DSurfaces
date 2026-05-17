using Gmsh.gmsh
using LinearAlgebra
using FerriteGmsh
using Ferrite
#= 
Dziuk surface: z^2 + y^2 + (x - z^2)^2 - 1.0 = 0
This script manually creates a 6-node 8-triangle mesh of a surface, and then refines it `ref` times.
=#

include("funcs_gensurface.jl")

function gensurfmesh_from0(ref::Int=0)

    # initial surface has 6 points:
    vertices0 = [
            1.25 0.0 -sqrt(3)/2 #1 : bottom pole
            1.0 0.0 0.0         #2 : equator point 1
            1.25 0.0 sqrt(3)/2  #3 : top pole
            0.0 1.0 0.0         #4 : equator point 2
            -1.0 0.0 0.0        #5 : equator point 3
            0.0 -1.0 0.0        #6 : equator point 4
        ]

    # Initial surface has 8 triangles:
    triangles0 = [
            1 2 4
            1 4 5
            1 5 6
            1 6 2
            3 4 2
            3 5 4
            3 6 5
            3 2 6
        ]
    
    # initialise gmsh
    gmsh.initialize()
    gmsh.model.add("initial_surface")
    
    # create discrete surface entity of dim=2
    surf_tag = gmsh.model.addDiscreteEntity(2)

    # create vertex coordinates
    node_coords::Vector = vec(vertices0')
    node_tags = collect(1:size(vertices0, 1))

    # add nodes to the model
    gmsh.model.mesh.addNodes( 
        2,  # dim
        surf_tag,  # tag of the entity to which the nodes belong
        node_tags,  # node tags
        node_coords  # node coordinates
    )

    # create triangle elements
    element_type = 2
    element_tags = collect(1:size(triangles0, 1))
    element_connectivity = vec(triangles0')
    
    # add elements to the model
    gmsh.model.mesh.addElements(
        2,  # dim
        surf_tag,  # tag of the entity to which the elements belong
        [element_type],  # element type
        [element_tags],  # element tags
        [element_connectivity]  # connectivity of the elements
    )

    # Refine mesh `ref` times
    if ref>=1
        for j in 1:(ref)
            gmsh.model.mesh.refine()

            node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
            etypes, _, enodes = gmsh.model.mesh.getElements()
            tri = enodes[findfirst(==(2), etypes)]
            triangles = reshape(tri, 3, :)

            # Project all vertices to the implicit surface
            project_allvertices!(node_coords)

            # apply Lloyd's algorithm to smooth the mesh
            vertices = reshape(node_coords, 3, :)
            vertices, triangles = surface_lloyd(vertices, triangles; iterations = 50, repel=true)

            # update the node coordinates in the gmsh model 
            node_coords .= vec(vertices)
            for (i, tag) in enumerate(node_tags)
                idx = 3*(i-1)
                gmsh.model.mesh.setNode(
                    Int(tag),
                    node_coords[idx+1:idx+3],
                    Float64[]
                )
            end              
        end
    end

    # generate mesh in gmsh
    gmsh.model.mesh.generate(2)

    # calculate min and max edge length
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    etypes, _, enodes = gmsh.model.mesh.getElements()
    min_edge_length, max_edge_length = edge_length_stats2(node_tags, node_coords, etypes, enodes)

    # write mesh into file and close gmsh
    gmsh.write("lowresmesh_$(ref).msh")
    gmsh.finalize()

    @info "Generated $(ref)-refined mesh. h=$(max_edge_length)."
    return min_edge_length, max_edge_length
end
