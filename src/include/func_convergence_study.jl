using Printf, PrettyTables

function convergence_study(sim_func::Function, refinements::AbstractVector{Int}, T::Float64, Δt::AbstractVector{Float64};
                            write_output::Bool=false, write_resuts::Bool=true, results_folder::String="")

    n = length(refinements)
    refs = Vector{Int}(undef, n)
    hs   = Vector{Float64}(undef, n)
    errs = Vector{Float64}(undef, n)
    T_sim = [min(T, dt * 20000) for dt in Δt]
    T_sim = floor.(T_sim, sigdigits=2)

    for (i, r) in enumerate(refinements)
        ref, h, err = sim_func(r, T_sim[i], Δt[i])
        refs[i] = ref
        hs[i]   = h
        errs[i] = err
        @info "refinement=$ref  h=$h  L2error=$err  T_final=$(T_sim)"
    end

    # Compute experimental order of convergence between consecutive levels
    eoc = fill(NaN, n)
    for i in 2:n
        eoc[i] = log(errs[i-1] / errs[i]) / log(hs[i-1] / hs[i])
    end

    eoc_str = [isnan(e) ? "-" : @sprintf("%.3f", e) for e in eoc]

    if write_resuts
        if any(T_sim .!= T)
            footnote_T="T=$(T_sim)"
        else footnote_T="T=$(T)"
        end
        if !isempty(results_folder)
            cd(results_folder)
        end
        open("results.md", "w") do f
            pretty_table(f, [refs hs Δt errs eoc_str];
                column_labels = ["Refinement", "Mesh Size", "Δt", "L2 Error", "EOC"],
                backend=:markdown,
                source_notes=footnote_T
            )
        end

    end
    
    return [refs, hs, Δt, errs, eoc, T_sim]
end