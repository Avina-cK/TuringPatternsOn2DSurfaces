include("final_sim_func.jl")

Refs = [0,1,2,3,4,5, 6, 7, 8]
n_refs = length(Refs)
Hs = zeros(n_refs)
L2_errors = copy(Hs)
i=1
for i in 1:n_refs
    _, Hs[i], L2_errors[i] = solve_LB_eq_on_Dzuiksurf(Refs[i])
end

results = hcat(Refs, Hs, L2_errors)

EOCs = [log(L2_errors[i-1]/L2_errors[i])/log(Hs[i-1]/Hs[i]) for i in 2:n_refs]
EOCs = [NaN; EOCs]

## Display results

using PrettyTables
pretty_table([Refs Hs L2_errors EOCs]; 
            column_labels= ["Refinement", "Mesh Size", "L2 Error", "EOC"],
            source_notes = "-Δₛu = f on the Dzuik surface. uₜᵣᵤₑ(x, y, z) = xy",
            #backend = :markdown
            )

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
