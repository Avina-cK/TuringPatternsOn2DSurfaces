##Set up system conditions and functions
struct SchnakenbergMaterial{T}
    Dᵤ::T
    Dᵥ::T
    a::T
    b::T
end

function reaction1(u, v, a, b, γ)
    return γ*(a - u + (u*u*v))
end

function reaction2(u, v, a, b, γ)
    return γ*(b - (u*u*v))
end

function steadystate(a,b)
    return [a+b; b/((a+b)^2)]
end

function fᵤ(a,b)
    return (b-a)/(b+a)
end

function fᵥ(a,b)
    return (a+b)^2
end

function gᵤ(a,b)
    -2*b/(a+b)
end

function gᵥ(a,b)
    -fᵥ(a,b)
end

function d_critical(a,b)
    f_u = fᵤ(a,b)
    f_v = fᵥ(a,b)
    g_u = gᵤ(a,b)
    g_v = gᵥ(a,b)
    return (((f_u * g_v) - (2.0 * f_v * g_u)) + sqrt(((f_u*g_v) - (2.0*f_v*g_u))^2 - (f_u*f_u*g_v*g_v)))/(f_u*f_u)
end

function γ_critical(a,b, k)
    d_c = d_critical(a,b)
    f_u = fᵤ(a,b)
    g_v = gᵥ(a,b)
    return (2.0*k*k*d_c)/((d_c*f_u)+g_v)
end
