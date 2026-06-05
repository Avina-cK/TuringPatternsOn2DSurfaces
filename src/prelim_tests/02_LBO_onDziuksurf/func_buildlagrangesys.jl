"""
    build_lagrange_system(Kₕ, F, C) → (A, b)

Assemble the (n+1) × (n+1) saddle-point system from:
  - Kₕ :: SparseMatrixCSC{Float64}  stiffness matrix (n × n)
  - F  :: Vector{Float64}            RHS force vector (n)
  - C  :: Vector{Float64}            constraint row (n), i.e. C[j] = ∫φ_j dΩ

Returns a dense augmented system A and RHS b.
"""
function build_lagrange_system(
    Kₕ::AbstractMatrix{Float64},   # n × n
    F::AbstractVector{Float64},   # n
    C::AbstractVector{Float64}    # n  (row vector stored as 1D array)
)
    n  = length(F)
    @assert size(Kₕ) == (n, n) "K_h must be n × n"
    @assert length(C) == n "C must have length n"

    # -- Augmented matrix A ∈ ℝ^(n+1 × n+1) -------------------------#
    #   A = [ K_h   C ]   <- last column is Cᵀ (stored as column)
    #       [  Cᵀ   0 ]   <- last row is C, bottom-right is 0
    C_col = sparse(C)
    A = [Kₕ        C_col        ;
         transpose(C_col) sparse([0.0])]

    # -- Augmented RHS b ∈ ℝ^(n+1) ------------------------------------#
    #   b = [ F; 0 ]
    b = [F; 0.0]
    return A, b
end
