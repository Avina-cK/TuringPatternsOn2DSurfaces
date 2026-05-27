using LinearAlgebra
using Ferrite

"""
Φₜ(cell_nodes, ξ)

Φₜ: ̂T -> T
maps a point in the reference triangle to one on the surface

"""
function Φₜ(cell_nodes::Vector{Vec{3, Float64}}, ξ::Vec{2})
    # Get shape functions
    N₁ = ξ[1]
    N₂ = ξ[2]
    N₃ = 1 - ξ[1] - ξ[2]

    # coordinates of nodes
    q1 = cell_nodes[1]
    q2 = cell_nodes[2]
    q3 = cell_nodes[3]

    # Map reference point to physical point
    x = (N₁ .* q1) + (N₂ .* q2) + (N₃ .* q3)
    return x
end

# Jacobian of Φₜ
function J_Φ(cell_nodes::Vector{Vec{3, Float64}})
    # Derivatives of shape functions
    dN₁dξ = Vec{2, Float64}((1.0, 0.0))
    dN₂dξ = Vec{2, Float64}((0.0, 1.0))
    dN₃dξ = Vec{2, Float64}((-1.0, -1.0))    

    # coordinates of nodes
    q1 = cell_nodes[1]
    q2 = cell_nodes[2]
    q3 = cell_nodes[3]

    # Jacobian matrix ∈ ℝ³ˣ²
    J = zeros(3,2)
    J[:,1] = dN₁dξ[1]*q1 + dN₂dξ[1]*q2 + dN₃dξ[1]*q3
    J[:,2] = dN₁dξ[2]*q1 + dN₂dξ[2]*q2 + dN₃dξ[2]*q3

    return J
end

# first fundamental form: G
function matrix_G(J::Matrix{Float64})
    G = J' * J
    return G
end

function matrix_invG(J::Matrix{Float64})
    G = matrix_G(J)
    detG = G[1,1]*G[2,2] - G[1,2]*G[2,1]
    invG = (1.0/detG) * [G[2,2] -G[2,1];
                            -G[1,2] G[1,1]]
    return invG
end

# area of cell
function dΩₑ(G::Matrix{Float64})
    return sqrt(det(G))
end

# Custom Mapping Type
# Reference: ξ ∈ RefTriangle ⊂ ℝ² -> Physical: x ∈ ℝ³
"""
    SurfaceMapping

Non-identity mapping for a scalar interpolation defined on a 2D reference
triangle (ξ ∈ ℝ²) whose physical nodes live in ℝ³.

The Jacobian J = ∂x/∂ξ ∈ ℝ³ˣ² is rectangular, so gradient mapping uses
the first fundamental form G = Jᵀ J ∈ ℝ²ˣ² instead of a square inverse.
"""
struct SurfaceMapping end

#--- Required(?): tell Ferrite which mapping LTI uses -----------#
#Ferrite.mapping_type(::LTI) = SurfaceMapping()

# Mapping Values Container
# Precomputed at each quadrature point during reinit!
"""
    SurfaceMappingValues

Holds all geometry quantities computed from the physical cell nodes
at a single reference point ξ:

  - J    : Jacobian ∂x/∂ξ ∈ ℝ³ˣ²
  - G    : First fundamental form  G = Jᵀ J  ∈ ℝ²ˣ²
  - invG : Inverse  G⁻¹ ∈ ℝ²ˣ²
  - detJ : Area element  √det(G)  (scalar)
"""
struct SurfaceMappingValues
    J:: Matrix{Float64}     # 3×2
    G:: Matrix{Float64}     # 2×2
    invG::Matrix{Float64}   # 2×2
    detJ::Float64           # √det(G)
end

"""
    compute_surface_mapping(cell_nodes, ξ) -> SurfaceMappingValues

Compute all mapping quantities for the surface element at reference
point ξ ∈ RefTriangle given physical node coordinates in ℝ³.
"""
function compute_surface_mapping(cell_nodes::Vector{Vec{3, Float64}}, ξ::Vec{2})
    J = J_Φ(cell_nodes)   # ℝ³ˣ² Jacobian (already defined above)
    G = matrix_G(J)        # Jᵀ J  ∈  ℝ²ˣ²
    invG = matrix_invG(J)     # G⁻¹
    detJ = dΩₑ(G)             # √det(G)
    return SurfaceMappingValues(J, G, invG, detJ)
end

# Value Mapping  (identity for scalar interpolations)
"""
    map_shape_value(::SurfaceMapping, φ_hat, mv) -> φ

Scalar shape function values are the same in reference and physical space.
"""
function map_shape_value(::SurfaceMapping, φ_hat::Number,mv::SurfaceMappingValues)
    return φ_hat  # φ(x) = φ_hat(ξ) - identity mapping
end

#-- Gradient Mapping -------------------------------------------#
# Derivation
# ----------
# The physical coordinates satisfy  x = Φₜ(ξ),  so
#
#     φ(x) = φ_hat(ξ)   =>   ∇_x φ · (∂x/∂ξ) = ∇_ξ φ_hatᵀ
#
# Writing J = ∂x/∂ξ ∈ ℝ³ˣ², the left-hand side gives:
#
#     Jᵀ ∇_x φ = ∇_ξ φ_hat
#
# The least-squares solution (pseudoinverse) is:
#
#     ∇_x φ = (Jᵀ)⁺ ∇_ξ φ_hat= J (Jᵀ J)⁻¹ ∇_ξ φ_hat= J G⁻¹ ∇_ξ φ_hat ∈ ℝ³
#
# This is the surface gradient — the component of ∇φ that lies in the
# tangent plane of the surface.

"""
    map_shape_gradient(::SurfaceMapping, ∇ξ_φ_hat, mv) -> ∇_x φ ∈ ℝ³

Map a reference-space gradient ∇_ξ φ_hat∈ ℝ² to the physical surface
gradient ∇_x φ ∈ ℝ³ via:

    ∇_x φ = J G⁻¹ ∇_ξ φ_hat

where J ∈ ℝ³ˣ² is the Jacobian and G = Jᵀ J is the first fundamental form.
"""
function map_shape_gradient(::SurfaceMapping, ∇ξ_φ_hat::Vec{2, Float64}, mv::SurfaceMappingValues)
    # Step 1:  G⁻¹ ∇_ξ φ_hat ∈ ℝ²
    v = mv.invG * [∇ξ_φ_hat[1], ∇ξ_φ_hat[2]]   # 2×2 times 2-vector
    # Step 2:  J · v  ∈ ℝ³
    physical_grad = mv.J * v               # 3×2 times 2-vector
    return Vec{3, Float64}((physical_grad[1], physical_grad[2], physical_grad[3]))
end


# Area Element
"""
    get_detJ(::SurfaceMapping, mv) -> √det(G)

Return the area scaling factor for integration on the surface element.

    ∫_Ωₑ f dΩ = ∫_T̂ f(Φₜ(ξ)) √det(G) dξ
"""
get_detJ(::SurfaceMapping, mv::SurfaceMappingValues) = mv.detJ
