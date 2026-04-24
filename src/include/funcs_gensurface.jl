#=
Source:
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
                        [i1;     new_idx; other],
                        [new_idx; i2;    other])
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

"""
Lloyd's algorithm
    surface_lloyd(coords, triangles; iterations, fixed_indices)
- `coords`: 3 x N matrix of vertex coordinates
- `triangles`: 3 x M matrix of triangle vertex indices
- `iterations`: number of Lloyd iterations to perform (default: 50)
- `fixed_indices`: vector of vertex indices to keep fixed during iterations (default: empty)
"""
function surface_lloyd(coords, triangles; iterations=50, fixed_indices=Int[])
    # convert coordinate matrix into a vector
    pts = [coords[:, i] for i in axes(coords, 2)]
    fixed_set = Set(fixed_indices)
    for _ in 1:iterations
        mass = zeros(length(pts))
        centroid = [zeros(3) for _ in 1:length(pts)]
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
        for i in eachindex(pts)
            # ignore fixed vertices
            i in fixed_set && continue
            #ensure no division by 0
            mass[i] == 0 && continue
            #new vertex is computed as an area-weighted centroid
            newp = centroid[i] / mass[i]
            #project to implicit surface
            px,py,pz = project_to_surf(newp...)
            #updates pts vector
            pts[i] .= (px,py,pz)
        end
    end
    #reshape vertex list appropriately
    newcoords = hcat(pts...)
    return newcoords, triangles
end