struct SurfaceMappingValues{T<:AbstractFloat}
    J :: Matrix{T} # 3×2 Jacobian
    G :: Matrix{T} # 2×2 first fundamental form JᵀJ
    invG :: Matrix{T} # 2×2 inverse of G
    detJ :: T # √det(G) = ‖col₁ × col₂‖
end

"""
    compute_surface_mapping(nodes, ξ) -> SurfaceMappingValues

Compute all geometry quantities at a single reference point ξ ∈ ℝ².
This is the per-quadrature-point workhorse called inside reinit!.
"""
function compute_surface_mapping(
        nodes :: AbstractVector{Vec{3,T}},
        ξ     :: Vec{2,T}
    ) where {T}
    J    = J_Φ(collect(Vec{3,T}, nodes))   # 3×2 Matrix
    G    = matrix_G(J)
    invG = matrix_invG(J)
    col1 = Vec{3,T}((J[1,1], J[2,1], J[3,1]))
    col2 = Vec{3,T}((J[1,2], J[2,2], J[3,2]))
    detJ = norm(col1 × col2)
    return SurfaceMappingValues{T}(J, G, invG, detJ)
end

"""
    get_detJ(mapping, mv) -> T

Extract the embedding determinant from a SurfaceMappingValues.
The `mapping` argument is accepted for API consistency but unused - all data lives in `mv`.
"""
@inline get_detJ(::SurfaceTriangleMapping, mv::SurfaceMappingValues) = mv.detJ


"""
    map_shape_gradient(mapping, ∇φ_ref, mv) -> Vec{3,T}

Push a reference-space gradient ∇φ_hat∈ ℝ² forward to a physical-space surface gradient ∈ ℝ³ using the covariant formula:

    ∇ₛφ = J · G⁻¹ · ∇φ_hat
    
where J ∈ ℝ³ˣ² is the Jacobian and G⁻¹ = (JᵀJ)⁻¹ ∈ ℝ²ˣ².
"""
function map_shape_gradient(
        ::SurfaceTriangleMapping,
        ∇φ_ref :: Vec{2,T},
        mv     :: SurfaceMappingValues{T}
    ) :: Vec{3,T} where {T}
    # G⁻¹ ∇φ_hat (2-vector in reference space)
    g = mv.invG * collect(∇φ_ref)           # 2×2 · [2]->[2]
    # J · (G⁻¹ ∇φ_hat) (3-vector in physical space)
    v = mv.J * g                             # 3×2 · [2]->[3]
    return Vec{3,T}((v[1], v[2], v[3]))
end
