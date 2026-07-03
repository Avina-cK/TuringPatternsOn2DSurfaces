#=
    This file solves the stationary problem
        -Δₛu = f 
    on a unit sphere using a manufactured solution,
        u_exact(x,y,z) = z
    and its corresponding f, RHS function,
        f_rhs(x,y,z) = 2z
    
    It outputs the results for 4 meshes as a table 
    containing the following:
    * refinement level
    * mesh size (longest edge length)
    * L2-error
    * Error of Convergence (EOC).
=#

include("final_sim_func.jl")

# Manufactured solution
u_exact(x, y, z) = z
# RHS function: f = -Δu = 2z
f_rhs(x, y, z) = 2 * z

cd(@__DIR__)
mkpath("results")
cd("results")
dir_to_sphere_meshes = "../../../sphere_meshes"
refs=[0, 1, 2, 3]
L2errors = zeros(length(refs))
Hs = zeros(length(refs))
for (i, ref) in enumerate(refs)
    Hs[i], L2errors[i] = solve_laplace_eq_on_sphere(ref, u_exact, f_rhs, false; dir_to_sphere_meshes=dir_to_sphere_meshes)
end
EOCs = [log(L2errors[i-1]/L2errors[i])/log(Hs[i-1]/Hs[i]) for i in 2:length(refs)]
EOCs = [NaN; EOCs]
using PrettyTables
pretty_table([refs Hs L2errors EOCs]; 
            column_labels= ["Refinement", "Mesh Size", "L2 Error", "EOC"],
            source_notes = "-Δₛu = f on a sphere. uₜᵣᵤₑ(x, y, z) = z, f(x, y, z) = 2z",
            #formatters = [fmt__latex_sn(1)], backend = :latex
            )
open("results.txt", "w") do f
            pretty_table(f, [refs Hs L2errors EOCs]; 
            column_labels= ["Refinement", "Mesh Size", "L2 Error", "EOC"],
            source_notes = "-Δₛu = f on a sphere. uₜᵣᵤₑ(x, y, z) = z, f(x, y, z) = 2z"
            )
end

using Plots
plot(Hs, L2errors,
    marker=:circle, scale=:log2,
    xlabel="Log of mesh size, log₂(h)", ylabel="Log of L₂-error",
    title="Convergence of L₂-error for -Δₛu = f on a sphere",
    legend=false, )