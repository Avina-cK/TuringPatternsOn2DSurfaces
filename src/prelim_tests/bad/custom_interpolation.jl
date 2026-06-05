#=
For RefTriangle
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
#using Ferrite
struct LTI <: ScalarInterpolation{RefTriangle, 1}
end

function Ferrite.reference_shape_value(ip::LTI, ξ::Vec{2}, shape_number::Int)
    ξ₁ = ξ[1]
    ξ₂ = ξ[2]
    shape_number == 1 && return ξ₁  #v1: 1 at ξ = (1,0)
    shape_number == 2 && return ξ₂  #v2: 2 at ξ = (0,1)
    shape_number == 3 && return 1 - ξ₁ - ξ₂  #v3: 3 at ξ = (0,0)
    throw(ArgumentError("no shape function $shape_number for interpolation $ip"))
end

function Ferrite.reference_coordinates(::LTI)
    return[
        Vec{2, Float64}((1.0, 0.0)), # v1
        Vec{2, Float64}((0.0, 1.0)), # v2
        Vec{2, Float64}((0.0, 0.0)), # v3
    ]    
end

Ferrite.vertexdof_indices(::LTI) = ((1,), (2,), (3,))

# Edge interior DOFs : linear triangle has no edge-interior DOFs
Ferrite.edgedof_interior_indices(::LTI) = ((), (), ())

# All edge DOFs (including vertices)
Ferrite.edgedof_indices(::LTI) = ((1, 2,), (2, 3,), (3, 1,))

# Face interior DOFs : no face-interior DOFs for linear triangle
Ferrite.facedof_interior_indices(::LTI) = ((),)

# All face DOFs (including vertices and edges)
Ferrite.facedof_indices(::LTI) = ((1, 2, 3),)

# Volume interior DOFs : none for 2D element
Ferrite.volumedof_interior_indices(::LTI) = ()

# Total number of base functions
Ferrite.getnbasefunctions(::LTI) = 3

# No need to adjust DOFs during distribution for linear elements
Ferrite.adjust_dofs_during_distribution(::LTI) = false

# Conformity : H1, aka continuous function values across cell borders
Ferrite.conformity(::LTI) = Ferrite.H1Conformity()

#= shape gradients
function Ferrite.reference_shape_gradient(ip::LTI, ξ::Vec{2}, shape_number::Int)
    shape_number == 1 && return Vec{2, Float64}((1.0, 0.0))
    shape_number == 2 && return Vec{2, Float64}((0.0, 1.0))
    shape_number == 3 && return Vec{2, Float64}((-1.0, -1.0))
    throw(ArgumentError("no shape function $shape_number for interpolation $ip"))
end
=#

## Tests
#=
using Ferrite, Test

# Create an instance
ip = LTI()

# Test basic properties
@test getnbasefunctions(ip) == 3
@test getrefshape(ip) == RefTriangle
@test Ferrite.getorder(ip) == 1
@test Ferrite.conformity(ip) == Ferrite.H1Conformity()

# Test shape function values at vertices
ξ1 = Vec{2, Float64}((1.0, 0.0))  # vertex 1
@test Ferrite.reference_shape_value(ip, ξ1, 1) == 1.0
@test Ferrite.reference_shape_value(ip, ξ1, 2) == 0.0
@test Ferrite.reference_shape_value(ip, ξ1, 3) == 0.0

ξ2 = Vec{2, Float64}((0.0, 1.0))  # vertex 2
@test Ferrite.reference_shape_value(ip, ξ2, 1) == 0.0
@test Ferrite.reference_shape_value(ip, ξ2, 2) == 1.0
@test Ferrite.reference_shape_value(ip, ξ2, 3) == 0.0

ξ3 = Vec{2, Float64}((0.0, 0.0))  # vertex 3
@test Ferrite.reference_shape_value(ip, ξ3, 1) == 0.0
@test Ferrite.reference_shape_value(ip, ξ3, 2) == 0.0
@test Ferrite.reference_shape_value(ip, ξ3, 3) == 1.0

# Test partition of unity
ξ_center = Vec{2, Float64}((1/3, 1/3))
sum_values = sum(i -> Ferrite.reference_shape_value(ip, ξ_center, i), 1:3)
@test sum_values ≈ 1.0
=#