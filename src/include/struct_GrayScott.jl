##Set up system conditions and functions
struct GrayScottMaterial{T}
    Dᵤ::T
    Dᵥ::T
    α::T
    β::T
end

function reaction1(u, v, α, β)
    return (α * (1.0 - u)) - (u * v^2)
end

function reaction2(u, v, α, β)
    return (-(α + β) * v) + (u * v^2)
end
