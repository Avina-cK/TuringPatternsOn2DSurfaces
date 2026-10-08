"""
Lloyd's algorithm
    surface_lloyd(coords, triangles; iterations, fixed_indices)
- `coords`: 3 x N matrix of vertex coordinates (required)
- `triangles`: 3 x M matrix of triangle vertex indices (required)
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
