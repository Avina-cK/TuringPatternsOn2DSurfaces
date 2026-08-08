##
using LinearAlgebra
include("../include/func_DzuikSurface.jl")

"""
Project to surface
    project_to_surf(x,y,z; max_iters=200, tol=1e-15, 
                dist_func::Function=Dziuk_surface,
                grad_func::Function=gradDziuk_surface)
- `x,y,z` : coordinates of points to project onto surface
- `max_iters` : maximum number of iterations to run projection
- `tol` : how close to the surface is close enough
- `dist_func` : distance function describing implicit surface
- `grad_func` : gradient of dist_func
"""
function project_to_surf(
    x,y,z; max_iters::Int=200, tol=1e-15,
    dist_func::Function=Dziuk_surface,
    grad_func::Function=gradDziuk_surface
    )
    px,py,pz = x,y,z
    for _ in 1:max_iters
        d_ist = dist_func(px,py,pz)
        if abs(d_ist) < tol
            break 
        end
        gx,gy,gz = grad_func(px,py,pz)
        norm_gₛ = gx^2 + gy^2 + gz^2
        px -= d_ist/norm_gₛ * gx
        py -= d_ist/norm_gₛ * gy
        pz -= d_ist/norm_gₛ * gz
    end
    return px,py,pz
end

"""
Lloyd's algorithm
    surface_lloyd(coords, triangles; iterations, fixed_indices)
- `coords`: 3 x N matrix of vertex coordinates
- `triangles`: 3 x M matrix of triangle vertex indices
- `iterations`: number of Lloyd iterations to perform (default: 50)
- `fixed_indices`: vector of vertex indices to keep fixed during iterations (default: empty)
"""
function surface_lloyd(coords, triangles; iterations=50, fixed_indices=Int[])
    pts = [coords[:, i] for i in axes(coords, 2)]
    fixed_set = Set(fixed_indices)
    n = length(pts)
    adj = [Int[] for _ in 1:n]

    for t in eachcol(triangles)
        a,b,c = t
        push!(adj[a], b); push!(adj[a], c)
        push!(adj[b], a); push!(adj[b], c)
        push!(adj[c], a); push!(adj[c], b)
    end

    for _ in 1:iterations
        mass = zeros(n)
        centroid = [zeros(3) for _ in 1:n]

        for t in eachcol(triangles)
            a,b,c = t
            #coordinates of triangle's vertices
            pa, pb, pc = pts[a], pts[b], pts[c]
            #compute centroid of triangle
            tri_centroid = (pa + pb + pc) / 3
            #compute area of triangle = |(b-a) × (c-a)|/2
            area = 0.5 * norm(cross(pb-pa, pc-pa))
            #for each vertex v, 
            for v in (a,b,c)
                #accumulate 'area * tri_centroid' ∀ Δles touching v
                centroid[v] += area * tri_centroid
                #accumulate 'area' ∀ Δles touching v
                mass[v] += area
            end
        end

        for i in 1:n
            # ignore fixed vertices
            i in fixed_set && continue
            # ensure no division by 0
            mass[i] == 0 && continue
            # new vertex is computed as an area-weighted centroid
            newp = centroid[i] / mass[i]
            # project to implicit surface
            px,py,pz = project_to_surf(newp...)
            # updates pts vector
            pts[i] .= (px,py,pz)
        end
    end
    # reshape vertex list appropriately
    newcoords = Vector(hcat(pts...)[1:end])
    return newcoords
end


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
