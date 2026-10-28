include("func_main_sim.jl")

ref= 6
T = 30000.0
Δt = 0.2
no_frames = 600
Dᵤ, Dᵥ= 0.00016, 0.00008

α_val = 0.020
β_vals = 0.0425:0.00125:0.06
β_len = length(β_vals)

b_ins = [0.8, 1.0]

for b_in in b_ins
    for β in β_vals
        given_m = GrayScottMaterial{Float64}(Dᵤ, Dᵥ, α_val, β)
        sim_gs_on_BloodCellmesh(given_m, ref, b_in, Δt, T, true, no_frames;
                        init_cond="polar",
                        polar_cap_given = 0.48
                        #noise_amp_given=0.005,
                        #init_const=[uv_const_given]
                    )
    end
end
