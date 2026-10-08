#=
using Random

surf_max = 1.26

u_stationary(x) = x[1]*x[2]

pert_centres = [
    [1.0, 0.0, 0.0],
    [-1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, -1.0, 0.0],
    #[0.3626239883908428, 0, -0.8970119084817058],
    #[0.3626239883908428, 0, 0.8970119084817058]
]

σ = 0.5
amp_u = 0.5
amp_v = 0.5
noise_amp = 0.0

#Random.seed!(314)

function u0_func(x1, x2, x3)
    x = [x1,x2,x3]
    bump_u = sum(exp(-norm(x - c)^2 / (2.0 * (σ^2))) for c in pert_centres)
    u0_val = 1.0 - amp_u * bump_u + noise_amp * (rand() - 0.5)
    return u0_val
end

function v0_func(x1, x2, x3)
    x = [x1,x2,x3]
    bump_v = sum(exp(-norm(x - c)^2 / (2σ^2)) for c in pert_centres)
    v0_val = amp_v * bump_v + noise_amp * (rand() - 0.5)
    return v0_val
end
=#
##

function setup_const_init_conditions!(u₀::Vector, cellvalues, dh::DofHandler, uv_const; noise_amp::Float64=0.01)
    u_const = uv_const[1]
    v_const = uv_const[2]
    n_basefuncs = getnbasefunctions(cellvalues)
    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        dofs = celldofs(cell)
        uₑ = @view u₀[dofs]
        rv₀ₑ = reshape(uₑ, (2, n_basefuncs))
        for i in 1:n_basefuncs
            rv₀ₑ[1, i] = u_const
            rv₀ₑ[2, i] = v_const
        end
    end
    Random.seed!(314)
    u₀ .+= noise_amp * rand(ndofs(dh))
    return
end;