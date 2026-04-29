#=
This is a test to simulate the transient heat equation on a plane. The equation is given by:
    ∂u/∂t - ∇ · (k ∇u) = f on Ω × (0, Tf]
    u = 0 on ∂Ω × (0, Tf]
where u is the temperature, k is the thermal conductivity, and f is a source term. 
The test will use a manufactured solution to verify the accuracy of the numerical method.
=#

using Ferrite
using SparseArrays
using Random
using WriteVTK
include("../include/funcs_error_analysis.jl")

k = 1.0e-3
Δt = 0.1
L₁ = 1.0
L₂ = 1.0  # domain = [-L₁,L₁] x [-L₂,L₂]

