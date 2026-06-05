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

## SurfaceTriangleMapping
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

##SurfaceMappingValues

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


##

const LTI = Lagrange{RefTriangle,1}

"""
    SurfaceCellValues{T, QR, IP} <: Ferrite.AbstractCellValues
 
Cell values for a scalar field defined on a surface (2-D triangle embedded in ℝ³).  Wraps a `SurfaceTriangleMapping` so that `reinit!` populates:
 
  * `N_surf[i,q]` : shape value  φᵢ(ξ_q)          (same as reference)
  * `∇N_surf[i,q]` : surface gradient  ∇ₛφᵢ ∈ ℝ³   via J·G⁻¹·∇̂φᵢ
  * `detJdV[q]` : area-element weight  √det(G)·w_q
"""

struct SurfaceCellValues{T <: AbstractFloat,
                         QR <: QuadratureRule{RefTriangle},
                         IP
                        } <: Ferrite.AbstractCellValues
    N_ref :: Matrix{T}
    ∇N_ref :: Matrix{Vec{2,T}}
    N_surf :: Matrix{T}
    ∇N_surf :: Matrix{Vec{3,T}}
    detJdV :: Vector{T}
    qr::QR
    ip::IP
    mapping::SurfaceTriangleMapping{T}
end

function SurfaceCellValues(
        ::Type{T},
        qr::QR,
        ip::IP
    ) where {T <: AbstractFloat,
            QR <: QuadratureRule{RefTriangle},
            IP <: Ferrite.ScalarIterpolation}
    
    no_q = length(Ferrite.getweights(qr))
    no_bfs = Ferrite.getnbasefunctions(ip)
    mapping = SurfaceTriangleMapping{T}(qr)
    
    N_ref = zeros(T, no_bfs, no_q)
    ∇N_ref = zeros(Vec{2,T}, no_bfs, no_q)

    for (q,ξ) in enumerate(Ferrite.getpoints(qr))
        for i in 1:no_bfs
            N_ref[i,q] = T(Ferrite.reference_shape_value(ip,ξ,i))
            ∇N_ref[i,q] = Vec{2,T}(Tuple(Ferrite.reference_shape_gradient(ip,ξ,i)))
        end
    end

    N_surf = zeros(T, no_bfs, no_q)
    ∇N_surf = zeros(Vec{3,T}, no_bfs, no_q)
    detJdV = zeros(T, no_q)

    return SurfaceCellValues{T,QR,IP}(
        N_ref, ∇N_ref,
        N_surf, ∇N_surf, detJdV,
        qr, ip, mapping
    )
end

SurfaceCellValues(qr::QuadratureRule{RefTriangle}, ip) = SurfaceCellValues(Float64, qr, ip)

Ferrite.getnquadpoints(cv::SurfaceCellValues) = length(cv.detJdV)
Ferrite.getnbasefunctions(cv::SurfaceCellValues) = size(cv.N_surf, 1)
Ferrite.getdetJdV(cv::SurfaceCellValues, q::Int) = cv.detJdV[q]
Ferrite.shape_value(cv::SurfaceCellValues, q::Int, i::Int) = cv.N_surf[i, q]
Ferrite.shape_gradient(cv::SurfaceCellValues, q::Int, i::Int) = cv.∇N_surf[i, q]
Ferrite.shape_value_type(::SurfaceCellValues{T}) where {T} = T
Ferrite.shape_gradient_type(::SurfaceCellValues{T}) where {T} = Vec{3,T}
Ferrite.function_interpolation(cv::SurfaceCellValues) = cv.ip
Ferrite.get_mapping(cv::SurfaceCellValues) = cv.mapping
Ferrite.geometric_interpolation(cv::SurfaceCellValues) = cv.mapping.ip

function Ferrite.reinit!(
        cv   :: SurfaceCellValues{T},
        cell :: Ferrite.CellIterator
    ) where {T}
    Ferrite.reinit!(cv, Ferrite.getcoordinates(cell))
    return cv
end

function Ferrite.reinit!(
        cv        :: SurfaceCellValues{T},
        cellnodes :: AbstractVector{<:Vec{3}}
    ) where {T}
 
    nodes  = collect(Vec{3,T}, cellnodes)
    no_bfs = Ferrite.getnbasefunctions(cv)
 
    for (q, (ξ, w)) in enumerate(
            zip(Ferrite.getpoints(cv.qr), Ferrite.getweights(cv.qr)))
 
        # per-quadrature-point geometry (from cf_02)
        mv = compute_surface_mapping(nodes, ξ) # SurfaceMappingValues{T}
 
        cv.detJdV[q] = mv.detJ * T(w) # √det(G) · w_q
 
        for i in 1:no_bfs
            cv.N_surf[i, q]  = cv.N_ref[i, q] 
            cv.∇N_surf[i, q] = map_shape_gradient(cv.mapping, cv.∇N_ref[i, q], mv)
        end
    end
    return cv
end

##
function setup_surface_cellvalues(grid::Ferrite.Grid{3, Triangle, Float64})
 
    qr = QuadratureRule{RefTriangle}(1)     # 1-point centroid rule (degree 1)
    ip = LTI()                              # Lagrange{RefTriangle,1}
    cv = SurfaceCellValues(Float64, qr, ip)
 
    println("SurfaceCellValues constructed:")
    println("  n_quadpoints  = ", Ferrite.getnquadpoints(cv))
    println("  n_basefuncs   = ", Ferrite.getnbasefunctions(cv))
    println("  grid cells    = ", Ferrite.getncells(grid))
    println("  grid nodes    = ", Ferrite.getnnodes(grid))
    println()
 
    # -- test loop: reinit! every cell and accumulate total surface area ----
    total_area = 0.0
    for cell in CellIterator(grid)
        Ferrite.reinit!(cv, cell)
        for q in 1:Ferrite.getnquadpoints(cv)
            total_area += Ferrite.getdetJdV(cv, q)   # ∑ √det(G)·w_q over all cells
        end
    end
    println("  Total surface area = ", total_area)
 
    # -- spot-check on cell 1 ----------------------------------------------
    cell1_nodes = Ferrite.getcoordinates(grid, 1)
    Ferrite.reinit!(cv, cell1_nodes)
 
    println("\nCell 1 spot-check (q=1):")
    println("  detJdV[1]    = ", cv.detJdV[1])
    println("  N_surf[:,1]  = ", cv.N_surf[:, 1])
    println("  ∇N_surf[1,1] = ", cv.∇N_surf[1, 1],
            "  (∈ ℝ³, surface gradient of φ₁)")
    println("  ∇N_surf[2,1] = ", cv.∇N_surf[2, 1])
    println("  ∇N_surf[3,1] = ", cv.∇N_surf[3, 1])
 
    # -- partition-of-unity check: ∑ N_surf[i,q] == 1 for all q ----------
    for q in 1:Ferrite.getnquadpoints(cv)
        s = sum(cv.N_surf[i,q] for i in 1:Ferrite.getnbasefunctions(cv))
        @assert isapprox(s, 1.0; atol=1e-12) "PoU failed at q=$q: ∑Nᵢ = $s"
    end
    println("\n  Partition-of-unity check passed ✓")
 
    return cv
end
