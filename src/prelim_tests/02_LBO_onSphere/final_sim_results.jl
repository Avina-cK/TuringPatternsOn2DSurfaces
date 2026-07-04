include("final_sim_func.jl")

## Parameters

u_exact(x, y, z) = z    # Manufactured solution

f_rhs(x, y, z) = 2 * z  # RHS function: f = -Δu = 2z

refs=[0, 1, 2, 3]   # Refinement levels

## Run simulation for all refinement levels
cd(@__DIR__)
mkpath("results")
cd("results")
dir_to_sphere_meshes = "../../../sphere_meshes"

L2errors = zeros(length(refs))
Hs = zeros(length(refs))

for (i, ref) in enumerate(refs)
    Hs[i], L2errors[i] = solve_laplace_eq_on_sphere(ref, u_exact, f_rhs, false; dir_to_sphere_meshes=dir_to_sphere_meshes)
end

EOCs = [log(L2errors[i-1]/L2errors[i])/log(Hs[i-1]/Hs[i]) for i in 2:length(refs)]
EOCs = [NaN; EOCs]

## Display results

using PrettyTables
pretty_table([refs Hs L2errors EOCs]; 
            column_labels= ["Refinement", "Mesh Size", "L2 Error", "EOC"],
            source_notes = "-Δₛu = f on a sphere. uₜᵣᵤₑ(x, y, z) = z, f(x, y, z) = 2z",
            backend = :markdown
            )

## Save results
using PrettyTables
open("results.txt", "w") do f
            pretty_table(f, [refs Hs L2errors EOCs]; 
            column_labels= ["Refinement", "Mesh Size", "L2 Error", "EOC"],
            source_notes = "-Δₛu = f on a sphere. uₜᵣᵤₑ(x, y, z) = z, f(x, y, z) = 2z"
            )
end

using Plots
plot(Hs, L2errors,
    marker=:circle,
    xlabel="Mesh size, h", ylabel="L₂-error",
    title="L₂-error vs Mesh size",
    legend=false, )

plot(Hs, L2errors,
    marker=:circle, scale=:log2,
    xlabel="Log of mesh size, log₂(h)", ylabel="Log of L₂-error",
    title="Convergence of L₂-error for -Δₛu = f on a sphere",
    legend=false, )
