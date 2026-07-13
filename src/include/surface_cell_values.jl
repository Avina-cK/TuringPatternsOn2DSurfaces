using Ferrite, LinearAlgebra
#=
ferrite-fem.github.io/Ferrite.jl/stable/devdocs/reference_cells/
Reference triangle: 
----------------+--------------------
Vertex numbers: | Vertex coordinates:
    2           |
    | \         | v1: 𝛏 = (1.0, 0.0)
    |   \       | v2: 𝛏 = (0.0, 1.0)
ξ₂^ |     \     | v3: 𝛏 = (0.0, 0.0)
  | 3-------1   |
  +--> ξ₁       |
----------------+--------------------
Edge numbers:   | Edge identifiers:
    +           |
    | \         | e1: (v1, v2)
    2   1       | e2: (v2, v3)
    |     \     | e3: (v3, v1)
    +---3---+   |
----------------+--------------------
Face numbers:   | Face identifiers:
    +           |
    | \         |
    |   \       | f1: (v1, v2, v3)
    |  1  \     |
    +-------+   | 
----------------+-------------------- 
=#

## Mapping primitives (generalised Φₜ / J_Φ / matrix_G / matrix_invG / dΩₑ)
"""
    Φₜ(cell_nodes, ξ) -> Vec{3,Float64}

Map reference point `ξ ∈ ℝ²` to the physical point `x ∈ ℝ³` 
using the three linear Lagrange shape functions of the reference triangle.
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
    x_surf = (N₁ .* q1) + (N₂ .* q2) + (N₃ .* q3)
    return Vec(x_surf[1], x_surf[2], x_surf[3])
end


"""
    J_Φ(cell_nodes) -> Matrix{Float64}  (3x2)

Jacobian of the map Φₜ. 
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
    matrix_G(J::Matrix{Float64}) -> Matrix{Float64}

First fundamental form `G = JᵀJ` of a Jacobian `J`.
"""
matrix_G(J::Matrix{Float64}) = J' * J

"""
    matrix_invG(J::Matrix{Float64}) -> Matrix{Float64}

Analytic inverse of the 2x2 first fundamental form `G = JᵀJ`.
"""
function matrix_invG(J::Matrix{Float64})
    G = matrix_G(J)
    detG = G[1, 1] * G[2, 2] - G[1, 2] * G[2, 1]
    return (1.0 / detG) * [G[2, 2] -G[2, 1]; -G[1, 2] G[1, 1]]
end

"""
    dΩₑ(G::Matrix{Float64}) -> Float64

Surface-element weight `√det(G) = ‖∂x/∂ξ₁ × ∂x/∂ξ₂‖`.
"""
function dΩₑ(G::Matrix{Float64})
    detG = G[1,1]*G[2,2] - G[1,2]*G[2,1]
    return sqrt(detG)
end

"""
    surface_normal(J::Matrix{T}) -> Vec{3,T}

Unit surface normal `(∂x/∂ξ₁ × ∂x/∂ξ₂)/‖·‖` for a 3x2 Jacobian `J`.
"""
function surface_normal(J::Matrix{T}) where {T}
    c1 = Vec{3, T}((J[1, 1], J[2, 1], J[3, 1]))
    c2 = Vec{3, T}((J[1, 2], J[2, 2], J[3, 2]))
    n = c1 × c2
    #n = Vec{3, T}((abs.(n[1]), abs.(n[2]), abs.(n[3])))
    return n / norm(n)
end

## SurfaceCellValues

"""
    SurfaceMapping{T}

Per-quadrature-point mapping data for the current cell, populated by
`reinit!`. `x[q]` is what `func_assembleKF_manual.jl` reads as `cv.mapping.x[q]`.
"""
struct SurfaceMapping{T}
    x::Vector{Vec{3, T}}      # spatial (physical) coordinate
    J::Vector{Matrix{T}}      # 3x2 Jacobian ∂x/∂ξ
    G::Vector{Matrix{T}}      # 2x2 first fundamental form
    invG::Vector{Matrix{T}}   # 2x2 inverse of G
    dΩ::Vector{T}             # √det(G)
    normal::Vector{Vec{3, T}} # unit surface normal
end

# https://ferrite-fem.github.io/Ferrite.jl/stable/topics/FEValues/#SimpleCellValues
"""
    SurfaceCellValues(qr::QuadratureRule{RefTriangle}, ip::Interpolation{RefTriangle}; T=Float64)

`CellValues`-style object for a scalar field on triangular surface elements
embedded in ℝ³. Reference shape values/gradients (`M`, `dNdξ`) are
precomputed once at construction; physical/tangential shape gradients
(`dNdx`), integration weights (`detJdV`), and the per-qpoint `mapping` are
(re)computed by `reinit!` for each cell.
"""
mutable struct SurfaceCellValues{IP <: Interpolation, T}
    ip::IP
    qr::QuadratureRule
    M::Matrix{T}             # n_basefuncs x n_qpoints, reference shape values
    dNdξ::Matrix{Vec{2, T}}  # n_basefuncs x n_qpoints, reference shape gradients
    dNdx::Matrix{Vec{3, T}}  # n_basefuncs x n_qpoints, physical/tangential shape gradients
    detJdV::Vector{T}        # n_qpoints, dΩₑ * w_q
    mapping::SurfaceMapping{T}
end

function SurfaceCellValues(qr::QuadratureRule, ip::Interpolation; T::Type = Float64)
    n_basefuncs = Ferrite.getnbasefunctions(ip)
    n_qpoints = Ferrite.getnquadpoints(qr)

    M = zeros(T, n_basefuncs, n_qpoints)
    dNdξ = zeros(Vec{2, T}, n_basefuncs, n_qpoints)
    for (q, ξ) in enumerate(Ferrite.getpoints(qr))
        for i in 1:n_basefuncs
            g, v = Ferrite.reference_shape_gradient_and_value(ip, ξ, i)
            dNdξ[i, q] = g
            M[i, q] = v
        end
    end

    dNdx = zeros(Vec{3, T}, n_basefuncs, n_qpoints)
    detJdV = zeros(T, n_qpoints)
    mapping = SurfaceMapping{T}(
        zeros(Vec{3, T}, n_qpoints),
        [zeros(T, 3, 2) for _ in 1:n_qpoints],
        [zeros(T, 2, 2) for _ in 1:n_qpoints],
        [zeros(T, 2, 2) for _ in 1:n_qpoints],
        zeros(T, n_qpoints),
        zeros(Vec{3, T}, n_qpoints),
    )
    return SurfaceCellValues{typeof(ip), T}(ip, qr, M, dNdξ, dNdx, detJdV, mapping)
end

Ferrite.getnbasefunctions(cv::SurfaceCellValues) = size(cv.M, 1)
Ferrite.getnquadpoints(cv::SurfaceCellValues) = size(cv.M, 2)
Base.@propagate_inbounds Ferrite.shape_value(cv::SurfaceCellValues, q::Int, i::Int) = cv.M[i, q]
Base.@propagate_inbounds Ferrite.shape_gradient(cv::SurfaceCellValues, q::Int, i::Int) = cv.dNdx[i, q]

"""
    getdetJdV(cv::SurfaceCellValues, q::Int) -> Float64

Integration weight `dΩₑ(ξ_q) * w_q` at quadrature point `q` for the last
`reinit!`-ed cell, as consumed by `assemble_elementKF!`.
"""
getdetJdV(cv::SurfaceCellValues, q::Int) = cv.detJdV[q]
getnormal(cv::SurfaceCellValues, q::Int) = cv.mapping.normal[q]

"""
    Ferrite.reinit!(cv::SurfaceCellValues, cell)

Reinitialize `cv` for `cell` (a `CellCache` from `CellIterator`). 
    For every quadrature point:
    * builds the Jacobian `J_Φ = [q1-q3 | q2-q3]` for linear triangles,
    * forms `G`, `G⁻¹`, `dΩₑ` via functions `matrix_G`,`matrix_invG`, `dΩₑ`,
    * pushes the reference shape gradients forward onto the physical tangent
      plane, `∇ₓNₐ = J G⁻¹ ∇ξNₐ`,
    * stores the physical coordinate, integration weight and unit normal.
"""
function Ferrite.reinit!(cv::SurfaceCellValues{IP, T}, cell) where {IP, T}
    cell_nodes = Ferrite.getcoordinates(cell)
    n_basefuncs = size(cv.M, 1)
    length(cell_nodes) == n_basefuncs ||
        throw(ArgumentError("cell has $(length(cell_nodes)) nodes but ip has $n_basefuncs basis functions"))

    weights = Ferrite.getweights(cv.qr)
    @inbounds for q in 1:Ferrite.getnquadpoints(cv)
        J = zeros(T, 3, 2)
        x = zero(Vec{3, T})
        for a in 1:n_basefuncs
            xa = cell_nodes[a]
            Na = cv.M[a, q]
            ga = cv.dNdξ[a, q]
            x += Na * xa
            J[1, 1] += xa[1] * ga[1]; J[1, 2] += xa[1] * ga[2];
            J[2, 1] += xa[2] * ga[1]; J[2, 2] += xa[2] * ga[2];
            J[3, 1] += xa[3] * ga[1]; J[3, 2] += xa[3] * ga[2]
        end
        G = matrix_G(J)
        invG = matrix_invG(J)
        dΩ = dΩₑ(G)

        cv.mapping.x[q] = x
        cv.mapping.J[q] = J
        cv.mapping.G[q] = G
        cv.mapping.invG[q] = invG
        cv.mapping.dΩ[q] = dΩ
        cv.mapping.normal[q] = surface_normal(J)
        cv.detJdV[q] = dΩ * weights[q]

        for a in 1:n_basefuncs
            ga = cv.dNdξ[a, q]
            g = invG * collect(ga)
            v = J * g
            cv.dNdx[a, q] = Vec{3, T}((v[1], v[2], v[3]))
        end
    end
    return cv
end

function Ferrite.function_value(cv::SurfaceCellValues{IP, T}, q_point::Int, u::AbstractVector) where {IP, T}
    n_basefuncs = Ferrite.getnbasefunctions(cv)
    length(u) == n_basefuncs ||
        throw(ArgumentError("length(u) = $(length(u)) but cv has $n_basefuncs basis functions"))
    val = zero(T)
    @inbounds for i in 1:n_basefuncs
        val += Ferrite.shape_value(cv, q_point, i) * u[i]
    end
    return val
end