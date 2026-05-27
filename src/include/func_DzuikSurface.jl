#=
Source for surface function:
Gerhard Dziuk. Finite Elements for the Beltrami operator on arbitrary surfaces, 
pages 142–155. Springer Berlin Heidelberg, Berlin, Heidelberg, 1988.
=#
"""
Dziuk_surface(x,y,z) -> dₛ::Float64
"""
function Dziuk_surface(x,y,z)
    return z^2 + y^2 + (x - z^2)^2 - 1.0
end

"""
gradDziuk_surface(x,y,z) -> (∇d₁, ∇d₂, ∇d₃)::Tuple{Float64, Float64, Float64}
"""
function gradDziuk_surface(x,y,z)
    gx = 2.0 *(x - z^2)
    gy = (2.0 *y)
    gz = 2.0 *z - (4.0 *z*(x - (z^2)))
    gx,gy,gz
end

@info "Loaded functions Dziuk_surface(x,y,z), gradDziuk_surface(x,y,z)"