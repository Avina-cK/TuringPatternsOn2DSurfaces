using Ferrite, Tensors

## Helper functions
"""
    Φₜ(cell_nodes, ξ) -> Vec{3,Float64}
 
Map reference point `ξ ∈ ℝ²` to the physical point `x ∈ ℝ³` using the three linear Lagrange shape functions of the reference triangle.
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
    return (N₁ .* q1) + (N₂ .* q2) + (N₃ .* q3)
end

"""
    J_Φ(cell_nodes) -> Matrix{Float64}  (3x2)
 
Jacobian of the map Φₜ.  Because N₃ = 1 - ξ₁ - ξ₂ the gradient contributions from node 3 are (-1, -1), giving the constant matrix
    J = [q₁ - q₃ | q₂ - q₃]   (columns are ∂x/∂ξ₁ and ∂x/∂ξ₂)
"""
function J_Φ(cell_nodes::Vector{Vec{3, Float64}})
    # coordinates of nodes
    q1 = cell_nodes[1]
    q2 = cell_nodes[2]
    q3 = cell_nodes[3]
    col1 = q1 - q3   # Vec{3}
    col2 = q2 - q3   # Vec{3}

    J = zeros(Float64, 3, 2)
    J[:, 1] = collect(col1)
    J[:, 2] = collect(col2)
    return J
end

"""
    matrix_G(J) -> Matrix{Float64}  (2x2)
 
First fundamental form  G = Jᵀ J.
"""
function matrix_G(J::Matrix{Float64})
    return J' * J
end

"""
    matrix_invG(J) -> Matrix{Float64}  (2x2)

Analytic inverse of the first fundamental form.
"""
function matrix_invG(J::Matrix{Float64})
    G    = matrix_G(J)
    detG = G[1,1]*G[2,2] - G[1,2]*G[2,1]
    return (1.0/detG) * [G[2,2]  -G[2,1];
                        -G[1,2]   G[1,1]]
end

"""
    dΩₑ(G) -> Float64
 
Surface-element weight: √det(G) = ‖∂x/∂ξ₁ x ∂x/∂ξ₂‖.
"""
function dΩₑ(G::Matrix{Float64})
    return sqrt(det(G))
end

##
# https://github.com/Ferrite-FEM/Ferrite.jl/blob/4eec1a347ebbaebd72db9c3d7e67b97d120b08b9/src/FEValues/GeometryMapping.jl

"""
    SurfaceTriangleMapping{T}
 
Pre-computed geometry mapping for a linear triangle in 3D (surface element, rdim = 2, sdim = 3).

∵ the Jacobian J ∈ ℝ³ˣ² is rectangular, the volume-element weight uses the embedding determinant  dΩ = √det(Jᵀ J) = ‖col₁ x col₂‖₂.

Stored per quadrature point (index `q`):
  * `Jmat[q]`   - plain 3x2 Julia matrix 
  * `Jtens[q]`  - `MixedTensor2{3,2,T}` (Ferrite/Tensors type)
  * `detJ[q]`   - embedding determinant √det(Jᵀ J)
  * `invG[q]`   - inverse first fundamental form (Jᵀ J)⁻¹ (2x2 matrix)
  * `x[q]`      - spatial coordinate of the quadrature point in ℝ³

All entries are filled by `reinit!(mapping, cell_nodes)`.
"""
struct SurfaceTriangleMapping{T<:AbstractFloat}
    qr :: QuadratureRule # quadrature rule on the triangle
    #TODO: replace with LTI?
    ip :: Lagrange{RefTriangle, 1} # linear scalar interpolation
 
    # Mutable state filled during reinit!
    Jmat :: Vector{Matrix{T}}             # 3×2 plain matrices
    Jtens :: Vector{MixedTensor2{3,2,T}}   # Ferrite Tensors type
    detJ :: Vector{T}
    invG :: Vector{Matrix{T}}
    x :: Vector{Vec{3,T}}
end

"""
    SurfaceTriangleMapping{T}(qr) -> SurfaceTriangleMapping
 
Construct the mapping, pre-allocating storage for every quadrature point in `qr`.
"""
function SurfaceTriangleMapping{T}(
        qr::QuadratureRule
    ) where {T<:AbstractFloat}

    ip     = Lagrange{RefTriangle,1}()
    n_qp   = Ferrite.getnquadpoints(qr)
    zero33 = zeros(T, 3, 2)
    zero22 = zeros(T, 2, 2)
    Jmat   = [copy(zero33) for _ in 1:n_qp]
    Jtens  = [zero(MixedTensor2{3,2,T}) for _ in 1:n_qp]
    detJ_v = zeros(T, n_qp)
    invG_v = [copy(zero22) for _ in 1:n_qp]
    x_v    = [zero(Vec{3,T}) for _ in 1:n_qp]
    return SurfaceTriangleMapping{T}(qr, ip, Jmat, Jtens, detJ_v, invG_v, x_v)
end

SurfaceTriangleMapping(qr::QuadratureRule) = SurfaceTriangleMapping{Float64}(qr)

## reinit! (to update mapping for a new cell)
@inline function getspatialcoord(
        ip :: Lagrange{RefTriangle, 1},
        ξ :: Vec{2, T},
        cell_nodes :: AbstractVector{Vec{3, T}}
    ) :: Vec{3, T} where {T}

    n_basefuncs = getnbasefunctions(ip)
    sx = zero(Vec{3, T})
    for i in 1:n_basefuncs
        M  = Ferrite.reference_shape_value(ip, ξ, i)    # scalar Mᵢ(ξ)
        sx += M * cell_nodes[i] # ∑ Mᵢ x̂ᵢ
    end
    return sx
end
"""
    reinit!(m::SurfaceTriangleMapping, cell_nodes::Vec)
 
Recompute all geometry quantities for the cell defined by `cell_nodes`(a length-3 vector of `Vec{3,T}`), called once per cell inside the element loop.
"""
function reinit!(
        m :: SurfaceTriangleMapping{T},
        cell_nodes :: AbstractVector{Vec{3,T}}
    ) where {T}

    qr_points = Ferrite.getpoints(m.qr)
    
    for (q,ξ) in enumerate(qr_points)
        J = J_Φ(cell_nodes)
        m.Jmat[q] = J

        # For in MixedTensor{3,2} form for Ferrite
        m.Jtens[q] = MixedTensor2{3,2,T}(
            (
                J[1,1], J[2,1], J[3,1],
                J[1,2], J[2,2], J[3,2]
            )
        )

        # embedding determinant: dΩ
        col1 = Vec{3,T}((J[1,1], J[2,1], J[3,1]))
        col2 = Vec{3,T}((J[1,2], J[2,2], J[3,2]))
        m.detJ[q] = norm(col1 × col2)

        m. invG[q] = matrix_invG(J)
        m.x[q] = getspatialcoord(m.ip, ξ, cell_nodes)
    end
    return m
end

@inline getjacobian(m::SurfaceTriangleMapping, q::Int) = m.Jtens[q]
@inline getjacobianmat(m::SurfaceTriangleMapping, q::Int) = m.Jmat[q]

@inline getdetJ(m::SurfaceTriangleMapping, q::Int) = m.detJ[q]

@inline getinvG(m::SurfaceTriangleMapping, q::Int) = m.invG[q]

@inline getnquadpoints(m::SurfaceTriangleMapping) = length(m.detJ)

##
"""
    integrate_scalar(f, m::SurfaceTriangleMapping) -> T
 
Compute  ∑_q  w_q · f(x_q) · detJ_q
where `f(x)` is a user-supplied scalar function of the physical coordinate.
"""
function integrate_scalar(
        f::Function,
        m::SurfaceTriangleMapping{T}
    ) where {T}

    weights = Ferrite.getweights(m.qr)
    n_qp    = Ferrite.getnquadpoints(m.qr)
    result  = zero(typeof(f(m.x[1])))   # infer return type from f
    for q in 1:n_qp
        result += weights[q] * f(m.x[q]) * m.detJ[q]
    end
    return result
end

"""
    surface_normal(m::SurfaceTriangleMapping{T}, q::Int) -> Vec{3,T}
 
Outward unit normal at quadrature point `q`.
 
    n = (∂x/∂ξ₁ × ∂x/∂ξ₂) / ‖∂x/∂ξ₁ × ∂x/∂ξ₂‖
"""
function surface_normal(m::SurfaceTriangleMapping{T}, q::Int) where {T}
    J    = m.Jmat[q]
    col1 = Vec{3,T}((J[1,1], J[2,1], J[3,1]))
    col2 = Vec{3,T}((J[1,2], J[2,2], J[3,2]))
    cross_vec = col1 × col2
    return cross_vec / norm(cross_vec)
end
 
##
#=
function _selftest()
    println("SurfaceTriangleMapping self-test\n")
 
    # -- T1: Flat triangle in the xy-plane: area = 0.5 ---#
    nodes_flat = [
        Vec{3,Float64}((0.0, 0.0, 0.0)),
        Vec{3,Float64}((1.0, 0.0, 0.0)),
        Vec{3,Float64}((0.0, 1.0, 0.0)),
    ]
    qr   = QuadratureRule{RefTriangle}(3)   # degree-3 rule
    m = SurfaceTriangleMapping(qr)
    reinit!(m, nodes_flat)
 
    area_flat = integrate_scalar(x -> 1.0, m)
    println("Flat  triangle area (expected 0.5): ", area_flat)
    @assert isapprox(area_flat, 0.5; atol=1e-12) "Area test failed"
 
    n_flat = surface_normal(m, 1)
    println("Flat  triangle normal  (expected [0,0,±1]): ", n_flat)
    @assert isapprox(abs(n_flat[3]), 1.0; atol=1e-12) "Normal test failed"
 
    # -- T2. Lifted triangle in 3D: vertices (0,0,0),(1,0,0),(0,1,1) ------#
    nodes_3d = [
        Vec{3,Float64}((0.0, 0.0, 0.0)),
        Vec{3,Float64}((1.0, 0.0, 0.0)),
        Vec{3,Float64}((0.0, 1.0, 1.0)),
    ]
    reinit!(m, nodes_3d)
 
    # Edge vectors: e1 = (1,0,0), e2 = (0,1,1)
    # Cross product: e1 × e2 = (0·1−0·1, 0·0−1·1, 1·1−0·0) = (0,−1,1)
    # Area = ½‖(0,−1,1)‖ = ½√2
    expected_area_3d = 0.5 * sqrt(2.0)
    area_3d = integrate_scalar(x -> 1.0, m)
    println("Lifted triangle area (expected $(round(expected_area_3d; digits=8))): ",
            round(area_3d; digits=8))
    @assert isapprox(area_3d, expected_area_3d; atol=1e-12) "3-D area test failed"
 
    # -- T3. Centroid location -----------------------------------------#
    # Centroid = (q1+q2+q3)/3 = (1/3, 1/3, 1/3)
    area_3d  = integrate_scalar(x -> 1.0, m) # Float64
    moment   = integrate_scalar(x -> x,   m) # Vec{3,Float64}
    centroid = moment / area_3d 
    println("Centroid (expected [1/3, 1/3, 1/3]): ", centroid)
    @assert isapprox(centroid, Vec{3,Float64}((1/3, 1/3, 1/3)); atol=1e-12)
    println("\nAll tests passed ✓")
end
 
_selftest()

=#