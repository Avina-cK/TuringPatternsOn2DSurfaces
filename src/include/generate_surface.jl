using Gmsh.gmsh
using LinearAlgebra

include("funcs_gensurface.jl")

#constraint vertices
v1 = vec([1.25 0.0 sqrt(3)/2])
v2 = vec([1.25 0.0 -sqrt(3)/2])
v3 = vec([1.0 0.0 0.0])
v4 = vec([-1.0 0.0 0.0])
v5 = vec([0.0 1.0 0.0])
v6 = vec([0.0 -1.0 0.0])
constraint_vs = [v1, v2, v3, v4, v5, v6]

# to be consistent with machine precision
constraint_pts = [collect(project_to_surf(v...)) for v in constraint_vs]

# TODO: create unique name 

function generate_Dzuik_surface(refinement::Float64 = 0.5; to_constrain::Bool=true, to_lloyd::Bool=true, la_iter::Int=200,)
    
    #--- initialise gmsh model --------------------#
    gmsh.initialize()

    modelname="SurfaceMesh_$(refinement)_ref"
    gmsh.model.add(modelname)

    #--- Create sphere of radius 1, centre at origin ---#
    gmsh.model.occ.addSphere(0,0,0,1.0) # add CAD sphere "object"
    gmsh.model.occ.synchronize()    # synchronise sphere to model

    #--- Set mesh size limits --------------------#
    gmsh.option.setNumber("Mesh.MeshSizeMin", refinement)
    gmsh.option.setNumber("Mesh.MeshSizeMax", refinement)

    #--- Generate spherical mesh ----------------------#
    gmsh.model.mesh.generate(2) # generate mesh of dim=2

    node_tags, node_coords,  = gmsh.model.mesh.getNodes()
    vertices = reshape(node_coords, 3, :)

    #--- Project all vertices to implicit surface ----------------------#
    for i in axes(vertices, 2)
        px, py, pz = project_to_surf(vertices[1,i], vertices[2,i], vertices[3,i])
        vertices[1,i] = px
        vertices[2,i] = py
        vertices[3,i] = pz
        gmsh.model.mesh.setNode(node_tags[i], [px,py,pz], [] )
    end

    _, _, elem_node_tags = gmsh.model.mesh.getElements(2)
    triangles = reshape(elem_node_tags[1], 3,:)

    #fixing constrainted vertices
    fixed_indices = Int[]

    if to_constrain
        for pc in constraint_pts
            # euclidean distance from every current vertex to this constraint point
            dists = [norm(vertices[:, i] - pc) for i in axes(vertices, 2)]
            idx   = argmin(dists)
            # snap the nearest vertex to the exact constraint position
            vertices[:, idx[1]] .= pc
            push!(fixed_indices, idx[1])
        end
        unique!(fixed_indices)
        idx_v1, idx_v2 = fixed_indices[1], fixed_indices[2]
        vertices, triangles = break_direct_edge!(vertices, triangles, idx_v1, idx_v2)
    end

    min_edge_length_before, max_edge_length_before = edge_length_stats(vertices, triangles)
    println("Min edge length before Lloyd's algorithm: ", min_edge_length_before)
    println("Max edge length before Lloyd's algorithm: ", max_edge_length_before)

    if to_lloyd
        lloyded="lloyded"
        #--- Run Lloyd's algorithm to smoothen mesh ----------------------#
        vertices, triangles = surface_lloyd(vertices, triangles; iterations=la_iter,  fixed_indices = fixed_indices,)
        
        min_edge_length_after, max_edge_length_after = edge_length_stats(vertices, triangles)
        println("Min edge length after Lloyd's algorithm: ", min_edge_length_after)
        println("Max edge length after Lloyd's algorithm: ", max_edge_length_after)
    else
        lloyded="notlloyded"
    end
    
    #--- clear whatever is in gmsh -------------------------#
    gmsh.model.removeEntities(gmsh.model.getEntities(), true)
    gmsh.model.mesh.clear()

    #--- initialise final surface as a Discrete Entity of dim=2 -#
    final_surf = gmsh.model.addDiscreteEntity(2)

    #--- Add vertices ----------------------------------------#
    N_v = size(vertices, 2) # no. of vertices
    gmsh.model.mesh.addNodes(
        2, #dim of entity
        final_surf,
        collect(1:N_v), #unique id numbers for nodes (nodeTags)
        vec(vertices)   #coordinates of vertices
    )

    #--- Add triangles --------------------------------------#
    N_t = size(triangles, 2)
    gmsh.model.mesh.addElements(
        2,
        final_surf,
        [2], #msh type of element (2 -> 3-node triangles)
        [collect(1:N_t)], #unique id numbers for triangles
        [vec(triangles)] 
    )

    #--- Write mesh into file ------#
    gmsh.write("Constrained_"*modelname*"_"*lloyded*".msh")
    #--- Close gmsh -----#
    gmsh.finalize()
end