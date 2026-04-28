"""
    L2_norm(u_given::Vector, cellvalues::CellValues, dh::DofHandler)
    Compute the L-2 norm of the given function
"""
function L2_norm(u_given::Vector, cellvalues::CellValues, dh::DofHandler)
    norm_sq = 0.0

    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        dofs = celldofs(cell)

        u_given_local = u_given[dofs]

        for q_point in 1:getnquadpoints(cellvalues)
            dΩᵣ = getdetJdV(cellvalues, q_point)
            uh_q = function_value(cellvalues, q_point, u_given_local)
            norm_sq += (uh_q)^2 * dΩᵣ
        end
    end

    return sqrt(norm_sq)
end

"""
    compute_L2_error(u_h::Vector, u_exact::Vector, cellvalues::CellValues, dh::DofHandler)
Compute the L-2 error between the simulated and true solution
    ‖ u - uₕ ‖_(L²)

"""
function compute_L2_error(u_h::Vector, u_exact::Vector, cellvalues::CellValues, dh::DofHandler)
    error_sq = 0.0

    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        dofs = celldofs(cell)

        u_h_local     = u_h[dofs]
        u_exact_local = u_exact[dofs]

        for q_point in 1:getnquadpoints(cellvalues)
            dΩᵣ = getdetJdV(cellvalues, q_point)
            uh_q = function_value(cellvalues, q_point, u_h_local)
            ue_q = function_value(cellvalues, q_point, u_exact_local)
            error_sq += (uh_q - ue_q)^2 * dΩᵣ
        end
    end

    return sqrt(error_sq)
end

"""
    compute_relative_error(u_h::Vector, u_exact::Vector, cellvalues::CellValues, dh::DofHandler)
        ‖(u-uₕ)‖_(L²) / ‖u‖_(L²)
"""
function relative_L2_error(u_h::Vector, u_exact::Vector, cellvalues::CellValues, dh::DofHandler)
    num = compute_L2_error(u_h, u_exact, cellvalues, dh)
    denom = L2_norm(u_exact, cellvalues, dh)
    return num/denom
end

"""
    EOC(h₁, h₂, L2err1, L2err2)
    Experimental Order of Convergence
           ⎧ ‖(u-uₕ₁)‖  ⎫      
        log⎥ —————————— ⎥ / log(h1/h2) 
           ⎩ ‖(u-uₕ₂)‖  ⎭
"""
function EOC(h₁, h₂, L2err1, L2err2)
    num = L2err1/L2err2
    denom = h₁/h₂
    return log2(num)/log2(denom)
end

@info "Loaded functions defining errors"