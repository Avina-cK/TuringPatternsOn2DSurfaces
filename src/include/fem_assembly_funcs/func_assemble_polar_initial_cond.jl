function setup_polar_initial_conditions!(u₀::Vector, cellvalues, dh::DofHandler; noise_amp::Float64=0.01, polar_cap=0.9)
    u₀ .= ones(ndofs(dh))
    u₀[2:2:end] .= 0.0
    n_basefuncs = getnbasefunctions(cellvalues)
    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        coords = getcoordinates(cell)
        dofs = celldofs(cell)
        uₑ = @view u₀[dofs]
        rv₀ₑ = reshape(uₑ, (2, n_basefuncs))
        Random.seed!(314)
        noise_u = noise_amp * rand(n_basefuncs)
        for i in 1:n_basefuncs
            if coords[i][3] > polar_cap
                rv₀ₑ[1, i] = 0.5 - noise_u[i]
                rv₀ₑ[2, i] = 0.25 + noise_u[i]
            else
                rv₀ₑ[1, i] -= noise_u[i]
                rv₀ₑ[2, i] += noise_u[i]
            end
        end
    end
    #u₀ .+= noise_amp * rand(ndofs(dh))
    return
end;

function setup_equator_initial_conditions!(u₀::Vector, cellvalues, dh::DofHandler; noise_amp::Float64=0.01, equiv_band=0.1)
    u₀ .= ones(ndofs(dh))
    u₀[2:2:end] .= 0.0
    n_basefuncs = getnbasefunctions(cellvalues)
    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        coords = getcoordinates(cell)
        dofs = celldofs(cell)
        uₑ = @view u₀[dofs]
        rv₀ₑ = reshape(uₑ, (2, n_basefuncs))
        Random.seed!(314)
        noise_u = noise_amp * rand(n_basefuncs)
        for i in 1:n_basefuncs
            if coords[i][3] < equiv_band/2 && coords[i][3] > -equiv_band/2
                rv₀ₑ[1, i] = 0.5 - noise_u[i]
                rv₀ₑ[2, i] = 0.25 + noise_u[i]
            else
                rv₀ₑ[1, i] -= noise_u[i]
                rv₀ₑ[2, i] += noise_u[i]
            end
        end
    end
    #u₀ .+= noise_amp * rand(ndofs(dh))
    return
end;

