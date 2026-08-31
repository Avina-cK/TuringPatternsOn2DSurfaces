using LinearAlgebra
using Symbolics
using Latexify

# Define the variables (do not use commas!)
@variables ζ η

Dζ = Differential(ζ)
Dη = Differential(η)

X = [ (η^2)+(cos(ζ)*sqrt(1-(η^2)))
      sin(ζ)*sqrt(1-(η^2))
      η ]

latexify(X)

∂X∂η = expand_derivatives.(Dη.(X))
∂X∂ζ = expand_derivatives.(Dζ.(X))

XζxXη = cross(∂X∂ζ, ∂X∂η)
latexify(XζxXη)

## Check is X_ζ × X_η ≠ 0 ∀ (ζ, η)

XζxXη = simplify.(XζxXη)

XζxXη_biult = build_function(XζxXη, ζ, η; expression=Val{false})[1]

ζ_range = range(0.0, 2π; length=500)
η_range = range(-1.0, 1.0; length=500)

min_norm  = Inf
min_ζη = (NaN, NaN)
zero_pts  = Tuple{Float64,Float64}[]

for ζ0 in ζ_range, η0 in η_range
    v = XζxXη_biult(ζ0, η0)
    norm_v = norm(v)
    if norm_v < min_norm
        global min_norm, min_ζη = norm_v, (ζ0, η0)
    end
    if norm_v < 1e-8
        push!(zero_pts, (ζ0, η0))
    end
end

println("min |X_ζ × X_η| = $min_norm at (ζ, η) = $min_ζη")
println(isempty(zero_pts) ? "No zero crossings found." :
        "Zero (or near-zero) at $(length(zero_pts)) grid points, at $(zero_pts)")

## First fundamental form, G

#κ₁ = norm(∂X∂ζ)
κ₁ = sqrt(1 - (η^2))
κ₁² = 1 - (η^2)

κ₂ = sqrt(1.0 + (η^2)*(1/(1-(η^2)) + 4 - 2*(cos(ζ)/sqrt(1-(η^2)))))
κ₂² = 1.0 + (η^2)*(1/(1-(η^2)) + 4 - 2*(cos(ζ)/sqrt(1-(η^2))))
dot_XζXη = simplify(dot(∂X∂ζ, ∂X∂η))
dot_XζXη = -0.5 * η * sin(ζ) * sqrt(1-(η^2))

latexify(dot_XζXη)

G = [κ₁² dot_XζXη;
    dot_XζXη κ₂²]
detG = det(G)

latexify(detG)

# H̄ = |G|G⁻¹
H̄ = [κ₂² -dot_XζXη;
    -dot_XζXη κ₁²]

∂H11∂x1 = expand_derivatives(Dζ(κ₂²))
#∂H12∂x2 = 0.5*sin(ζ)*expand_derivatives(Dη(sqrt(1 - (η^2))*η))
∂H12∂x2 = 0.5*sin(ζ)*(1 - (2*(η^2)))/(sqrt(1 - (η^2)))
∂H21∂x1 = expand_derivatives(Dζ(-dot_XζXη))
∂H22∂x2 = expand_derivatives(Dη(κ₁²))

∂detG∂x1 = expand_derivatives(Dζ(detG))

#∂detG∂x2 = expand_derivatives(Dη(detG))
∂detG∂x2 = (η/(sqrt(1- (η^2))))*(2*(3*(η^2) -2)*cos(ζ) + sqrt(1- (η^2))*(1 - 2*(η^2))*(8 - 0.5*(sin(ζ)^2)))
latexify(∂detG∂x2)

ã = κ₂²
b̃ = -dot_XζXη
c̃ = κ₁²

#d̃ = (1/detG)*(∂H11∂x1 + ∂H12∂x2 -0.5*((κ₂²*∂detG∂x2) - (dot_XζXη*∂detG∂x1)))
# ∂H11∂x1 + ∂H12∂x2 = sin(ζ)*(η^2 +0.5)/(sqrt(1 - η^2))

d̃ = (1/detG)*((sin(ζ)*(η^2 +0.5)/(sqrt(1 - η^2))) -0.5*((κ₂²*∂detG∂x2) - (dot_XζXη*∂detG∂x1)))

sum_dHe = ∂H21∂x1 + ∂H22∂x2
sum_2e =  ((κ₁²)*∂detG∂x1) - (dot_XζXη*∂detG∂x2)

ẽ = (1/detG)*(sum_dHe - 0.5*sum_2e)

