include("func_final_sim.jl")
include("../../include/func_convergence_study.jl")

##
refinements = 0:6
T  = 5.0
Δt = 0.05

using JLD2
cd(@__DIR__)
@load "../../Dziuk_surf_meshes/Dzuik_surf_meshsizes.jld2" Dzuik_surf_Hs

Hs = Dzuik_surf_Hs[refinements.+1]

# Δt ∝ h^2 
Δts = min.(Δt, Δt .* (Hs.^2))
Δts = floor.(Δts, sigdigits=1)

k_given = 1.0
write_output = false

sim_func = (r, T, Δt) -> sim_THE_onDziuksurf(r, T, Δt, k_given; tosavesol=write_output)

##
cd(@__DIR__)
cd("results/")
results, eoc, T_sim = convergence_study(sim_func, refinements, T, Δts; results_folder="results/")
@save "results_data.jld2" results, eoc, T_sim

##
pretty_table([results[1] results[2] Δts results[3] eocs];
        column_labels = ["Refinement", "Mesh Size", "Δt", "L2 Error", "EOC"],
        source_notes = "-Δₛu = f on Dzuik surface. uₜᵣᵤₑ(x, y, z,t) = x * y * sin(t + π/2), T=$(T)",
        backend = :markdown
    )

L2_errors = results[3]
Hs = results[2]

using Plots
plot(Hs, L2_errors,
    marker=:circle,
    xlabel="Mesh size, h", ylabel="L₂-error",
    title="L₂-error vs Mesh size",
    legend=false, )

plot(Hs, L2_errors,
    marker=:circle, scale=:log10,
    xlabel="Log of mesh size, log₂(h)", ylabel="Log of L₂-error",
    title="Convergence of L₂-error for -Δₛu = f on a sphere",
    legend=false, )
