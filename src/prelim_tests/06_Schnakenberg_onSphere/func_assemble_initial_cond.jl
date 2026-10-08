function setup_random_initcond!(u₀::Vector, cellvalues, dh::DofHandler,
                                 u0::Float64, v0::Float64;
                                 noise_amp::Float64=0.05)
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
                rv₀ₑ[1, i] = u0
                rv₀ₑ[2, i] = v0
        end
    end
    u₀ .+= noise_amp * rand(ndofs(dh))
    return
end;

function setup_initial_conditions!(u₀::Vector, cellvalues, dh::DofHandler)
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
    u₀ .+= 0.05 * rand(ndofs(dh))
    return
end;

# ratio (l-m)! / (l+m)!  for 0 <= m <= l
function norm_ratio(l::Int, m::Int)
    @assert 0 <= m <= l "norm_ratio requires 0 <= m <= l"
    m == 0 && return 1.0
    r = 1.0
    # (l-m)!/(l+m)! = 1 / [(l-m+1)(l-m+2)...(l+m)]
    for k in (l - m + 1):(l + m)
        r /= k
    end
    return r
end

function setup_rand_sh_initcond!(
            u₀::Vector, cellvalues, dh::DofHandler,
            u0::Float64, v0::Float64;
            noise_amp::Float64=0.01,
            N::Int=8, n::Int=4
         )
    u₀ .= ones(ndofs(dh))
    u₀[2:2:end] .= 0.0
    n_basefuncs = getnbasefunctions(cellvalues)

    lm_pairs = [(l, m) for l in [n] for m in 1-N:N-1 if abs(m) <= n]
    print(lm_pairs)
    n_modes = length(lm_pairs)

    # seperate coeffs for seperate noise
    coeffs_u = randn(n_modes)
    coeffs_v = randn(n_modes)

    #https://mathworld.wolfram.com/AssociatedLegendrePolynomial.html
    function P_lm(l::Int, m::Int, x::Float64)
        pmm = 1.0
        if m > 0
            somx2 = sqrt((1.0 - x) * (1.0 + x))
            fact = 1.0
            for i in 1:m
                pmm *= -fact * somx2
                fact += 2.0
            end
        end
        l == m && return pmm
        pmmp1 = x * (2m + 1) * pmm
        l == m + 1 && return pmmp1
        pll = 0.0
        for ll in (m + 2):l
            pll = ((2ll - 1) * x * pmmp1 - (ll + m - 1) * pmm) / (ll - m)
            pmm, pmmp1 = pmmp1, pll
        end
        return pll
    end

    # real spherical harmonic Yˡₘ(θ, φ)
    function real_Yˡₘ(l::Int, m::Int, θ::Float64, φ::Float64)
        am  = abs(m)
        Plm = P_lm(l, am, cos(θ))
        norm = sqrt((2l + 1) / (4π) * norm_ratio(l, am))
        if m == 0
            return norm * Plm
        elseif m > 0
            return sqrt(2) * norm * Plm * cos(m * φ)
        else
            return sqrt(2) * norm * Plm * sin(am * φ)
        end
    end

    # truncated spherical-harmonic sum evaluated at a point on the sphere
    function sh_noise(x::Float64, y::Float64, z::Float64, coeffs::Vector{Float64})
        r = sqrt(x^2 + y^2 + z^2)
        θ = acos(clamp(z / r, -1.0, 1.0))
        φ = atan(y, x)
        s = 0.0
        for (k, (l, m)) in enumerate(lm_pairs)
            s += coeffs[k] * real_Yˡₘ(l, m, θ, φ)
        end
        return s
    end

    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        coords = getcoordinates(cell)
        dofs = celldofs(cell)
        uₑ = @view u₀[dofs]
        rv₀ₑ = reshape(uₑ, (2, n_basefuncs))
        for i in 1:n_basefuncs
            x, y, z = coords[i][1], coords[i][2], coords[i][3]
            rv₀ₑ[1, i] = u0 + noise_amp * sh_noise(x, y, z, coeffs_u)
            rv₀ₑ[2, i] = v0 + noise_amp * sh_noise(x, y, z, coeffs_v)
        end
    end
    return nothing
    
end


