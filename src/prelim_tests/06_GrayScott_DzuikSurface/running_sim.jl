include("func_main_sim.jl")

##Parameters
ref=7
T = 10000.0
Δt = 1.0

Dᵤ, Dᵥ= 0.00016, 0.00008
α = 0.07
β = 0.062

## Simulation
given_m = GrayScottMaterial{Float64}(Dᵤ, Dᵥ, α, β)

sim_gs_on_Dziukmesh(given_m, ref, Δt, T)

##

continue_T = 10000.0

resultsfolder_name = "GS_Dziuk_surf_$(ref)_$(fts(Dᵤ))_$(fts(Dᵥ))_$(fts(α))_$(fts(β))"
saved_state = "GS_Dziuk_surf_state_$(Int(T)).jld2"
restart_sol = "results\\$(resultsfolder_name)\\$(saved_state)"
sim_gs_on_Dziukmesh(
    given_m, ref, Δt, continue_T;
    restart_from = restart_sol)
