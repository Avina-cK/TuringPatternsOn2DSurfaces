using JLD2

"""Save final state so it can be reloaded as an initial condition for a restart run."""
function save_state(uₜ::Vector{Float64}, t::Float64, iₜ::Int, filepath::String)
    jldsave(filepath; u=uₜ, t=t, iter=iₜ)
end

"""Load a previously saved state."""
function load_state(filepath::String)
    data = load(filepath)
    return data["u"], data["t"], data["iter"]
end