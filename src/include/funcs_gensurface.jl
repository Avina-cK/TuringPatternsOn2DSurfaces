#=
Source for surface function:
Gerhard Dziuk. Finite Elements for the Beltrami operator on arbitrary surfaces, 
pages 142–155. Springer Berlin Heidelberg, Berlin, Heidelberg, 1988.
=#

function Dziuk_surface(x,y,z)
    z^2 + y^2 + (x - z^2)^2 - 1.0
end

function gradDziuk_surface(x,y,z)
    gx = 2.0 *(x - z^2)
    gy = (2.0 *y)
    gz = 2.0 *z - (4.0 *z*(x - (z^2)))
    gx,gy,gz
end

function project_to_surf(x,y,z)
    px,py,pz = x,y,z
    for _ in 1:100
        d_ist = Dziuk_surface(px,py,pz)
        if abs(d_ist) < 1e-10 
            break 
        end
        gx,gy,gz = gradDziuk_surface(px,py,pz)
        norm_gₛ = gx^2 + gy^2 + gz^2
        px -= d_ist/norm_gₛ * gx
        py -= d_ist/norm_gₛ * gy
        pz -= d_ist/norm_gₛ * gz
    end
    return px,py,pz
end

function project_allvertices!(verts; max_itersg::Int=100, tolg::Float64=1e-10)
    n = length(verts) ÷ 3
    for i in 1:n
        idx = 3*(i-1)
        px = verts[idx + 1]
        py = verts[idx + 2]
        pz = verts[idx + 3]
        verts[idx + 1], verts[idx + 2], verts[idx + 3] = project_to_surf(px,py,pz; max_iters=max_itersg, tol=tolg)
    end
    return nothing
end

"""
    break_direct_edge!(vertices, triangles, i1, i2) -> (vertices, triangles)
If vertex indices `i1` and `i2` share one or more triangle edges, split each such
triangle by inserting a midpoint vertex (projected back onto the implicit surface).
Returns updated `vertices` (3×N) and `triangles` (3×M) matrices.
"""
function break_direct_edge!(vertices, triangles, i1, i2)
    shared_tris = [t for t in axes(triangles, 2)
                   if i1 ∈ triangles[:, t] && i2 ∈ triangles[:, t]] 
    isempty(shared_tris) && return vertices, triangles
    # Project the midpoint back onto the implicit surface
    mid_raw = (vertices[:, i1] .+ vertices[:, i2]) ./ 2
    mid = collect(project_to_surf(mid_raw[1], mid_raw[2], mid_raw[3]))
    new_idx = size(vertices, 2) + 1
    vertices = hcat(vertices, mid)
    # Keep non-shared triangles; split each shared triangle into two
    kept_mask = trues(size(triangles, 2))
    kept_mask[shared_tris] .= false
    new_tris = triangles[:, kept_mask]
    for t in shared_tris
        tri = triangles[:, t]
        other = setdiff(tri, [i1, i2])[1]   # third vertex of this triangle
        # Each half retains one of the original endpoints and the new midpoint
        new_tris = hcat(new_tris,
                        [i1; new_idx; other],
                        [new_idx; i2; other])
    end
    return vertices, new_tris
end

"""
    edge_length_stats(vertices::Matrix{Float64}, triangles::Matrix{<:Integer}) -> (min_len, max_len)
"""
function edge_length_stats(vertices::Matrix{Float64}, triangles::Matrix{<:Integer})
    max_len = 0.0
    min_len = Inf
    for t in eachcol(triangles)
        a, b, c = t
        pa, pb, pc = vertices[:, a], vertices[:, b], vertices[:, c]
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

"""
Lloyd's algorithm
    surface_lloyd(coords, triangles; iterations, fixed_indices)
- `coords`: 3 x N matrix of vertex coordinates
- `triangles`: 3 x M matrix of triangle vertex indices
- `iterations`: number of Lloyd iterations to perform (default: 50)
- `fixed_indices`: vector of vertex indices to keep fixed during iterations (default: empty)
"""
function surface_lloyd(coords, triangles; iterations=50, fixed_indices=Int[], repel::Bool=false, ϵ=1e-15, λₗ=0.0001)
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

        if repel
            # move points that are too closed to each other
            for i in 1:n
                i in fixed_set && continue
                move = zeros(3)
                for j in adj[i]
                    diff = pts[i] - pts[j]
                    dist2= dot(diff,diff) + ϵ
                    move .+= diff ./ dist2
                end
                pts[i] .+= λₗ * move
            end
            
            # project to implicit surface
            for i in 1:n
                i in fixed_set && continue    
                px,py,pz = project_to_surf(pts[i]...)
                pts[i] .= (px,py,pz)
            end
        end
    end
    # reshape vertex list appropriately
    newcoords = hcat(pts...)
    return newcoords, triangles
end
