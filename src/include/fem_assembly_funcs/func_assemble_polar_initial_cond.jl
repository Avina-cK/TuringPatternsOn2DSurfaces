function setup_polar_initial_conditions!(u₀::Vector, cellvalues, dh::DofHandler; noise_amp::Float64=0.01)
    u₀ .= ones(ndofs(dh))
    u₀[2:2:end] .= 0.0
    n_basefuncs = getnbasefunctions(cellvalues)
    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        coords = getcoordinates(cell)
        dofs = celldofs(cell)
        uₑ = @view u₀[dofs]
        rv₀ₑ = reshape(uₑ, (2, n_basefuncs))
        for i in 1:n_basefuncs
            if coords[i][3] > 0.9
                rv₀ₑ[1, i] = 0.5
                rv₀ₑ[2, i] = 0.25
            end
        end
    end
    Random.seed!(314)
    u₀ .+= noise_amp * rand(ndofs(dh))
    return
end;
