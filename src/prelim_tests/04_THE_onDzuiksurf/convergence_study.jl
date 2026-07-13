include("func_final_sim.jl")
include("../../include/func_convergence_study.jl")

##
refinements = 0:5
T  = 2.0 * π
Δt = 0.1

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
results = convergence_study(sim_func, refinements, T, Δts; results_folder="results/")

##
#results = refs, hs, Δt, errs, eoc, T_sim
@save "results_data.jld2" results

pretty_table([results[1] results[2] results[3] results[4] results[5] results[6]];
        column_labels = ["Refinement", "Mesh Size", "Δt", "L2 Error", "EOC", "Tf"],
        source_notes = "-Δₛu = f on Dzuik surface. uₜᵣᵤₑ(x, y, z,t) = x * y * sin(t + π/2), T=$(T)",
        backend = :markdown
    )

L2_errors = results[4]
Hs = results[2]

using Plots
plot(Hs, L2_errors,
    marker=:circle,
    xlabel="Mesh size, h", ylabel="L₂-error",
    title="L₂-error vs Mesh size",
    legend=false, )
savefig("results_L2vsH.png")

plot(Hs, L2_errors,
    marker=:circle, scale=:log2,
    xlabel="Log of mesh size, log₂(h)", ylabel="Log of L₂-error",
    title="Error convergence for ∂ₜu -Δₛu = f on Dzuik surface",
    legend=false, )
savefig("results_logL2vslogH.png")