"""
Project to surface
    project_to_surf(x,y,z; max_iters=200, tol=1e-15, 
                dist_func::Function=given_surface,
                grad_func::Function=given_grad_surface)
- `x,y,z` : coordinates of points to project onto surface
- `max_iters` : maximum number of iterations to run projection
- `tol` : how close to the surface is close enough
- `dist_func` : distance function describing implicit surface
- `grad_func` : gradient of dist_func
"""
function project_to_surf(
    x,y,z; max_iters::Int=200, tol=1e-15,
    dist_func::Function=given_surface,
    grad_func::Function=given_grad_surface
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
