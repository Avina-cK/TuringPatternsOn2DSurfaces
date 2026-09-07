include("func_main_sim.jl")

γ(α, β) = sqrt(Complex(1.0 - (4.0*((α+β)^2)/α)))

w0 = [1.0; 0.0]
w1(α, β) = [0.5*(1.0+γ(α, β)); 0.5*(α/(α+β))*(1.0-γ(α, β))]
w2(α, β) = [0.5*(1.0-γ(α, β)); 0.5*(α/(α+β))*(1.0+γ(α, β))]

α_vals = 0.005:0.005:0.080
β_vals = 0.045:0.00125:0.070
β_vals_half = 0.045625:0.00125:0.070625
##Parameters
ref=6
T = 80000.0
Δt = 1.0
no_frames = 400
Dᵤ, Dᵥ= 0.00016, 0.00008
α = 0.030
βvals = [0.054375]

## Simulation
α_len = length(α_vals)
β_len = length(βvals)

for i in 1:β_len
    given_m = GrayScottMaterial{Float64}(Dᵤ, Dᵥ, α, βvals[i])
    sim_gs_on_Dziukmesh(given_m, ref, Δt, T, true, no_frames;
                     #init_const=[uv_const_given]
                     )
end

##
#
continue_T = 80000.0
β_val = 0.065
resultsfolder_name = "GS_Dziuk_surf_$(ref)_$(fts(Dᵤ))_$(fts(Dᵥ))_$(fts(α))_$(fts(β_val))"
saved_state = "GS_Dziuk_surf_state_$(Int(4*T)).jld2"
restart_sol = "results\\$(resultsfolder_name)\\$(saved_state)"
given_m = GrayScottMaterial{Float64}(Dᵤ, Dᵥ, α, β_val)
sim_gs_on_Dziukmesh(
    given_m, ref, Δt, continue_T, true, no_frames;
    restart_from = restart_sol)

#
##

#=


uv_steadystate = w1(α, β)
if imag(uv_steadystate)==[0.0;0.0]
    u_steadystate = real(w1(α, β)[1])
    v_steadystate = real(w1(α, β)[2])
    uv_const_given = [u_steadystate; v_steadystate]

else
    @info "steady state is complex"
    uv_const_given = [0.0; 0.0]
end

=#