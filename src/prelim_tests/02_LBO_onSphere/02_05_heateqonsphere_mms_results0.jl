using Ferrite, FerriteGmsh, SparseArrays, LinearAlgebra, Printf

## Manufactured solution:
ϕ(x) = x[1]*x[2]

const k_decay = 0.5
u_ms(x,t)  =  exp(-k_decay * t) * ϕ(x)
ut_ms(x,t) = -k_decay * exp(-k_decay * t) * ϕ(x)

## Elemental matrix assembly
function assemble_element_mass!(Me, cv::CellValues)
    fill!(Me, 0.0)
    n_bfs = getnbasefunctions(cv)
    for q_point in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, q_point)
        for i in 1:n_bfs
            φᵢ = shape_value(cv, i, q_point)
            for j in 1:n_bfs
                φⱼ = shape_value(cv, j, q_point)
                Me[i, j] += φᵢ * φⱼ * dΩ
            end
        end
    end
end

function assemble_element_stiffness!(Ke, cv::CellValues, k)
    n_bfs = getnbasefunctions(cv)
    fill!(Ke, 0)
    for q_point in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, q_point)
        for i in 1:n_bfs
            ∇φᵢ = shape_gradient(cv, i, q_point)
            for j in 1:n_bfs
                ∇φⱼ = shape_gradient(cv, j, q_point)
                Ke[i, j] += k * (∇φᵢ ⋅ ∇φⱼ) * dΩ
            end
        end
    end
end

## Global assembly

function assemble_system!(M, K, cv, dh, k_sys::Float64)
    n_bfs = getnbasefunctions(cv)
    Me = zeros(n_bfs, n_bfs)
    Ke = zeros(n_bfs, n_bfs)
    aM = start_assemble(M)
    aK = start_assemble(K)

    for cell in CellIterator(dh)
        reinit!(cv, cell)
        assemble_element_mass!(Me, cv)
        assemble_element_stiffness!(Ke, cv, k_sys)
        assemble!(aM, celldofs(cell), Me)
        assemble!(aK, celldofs(cell), Ke)
    end
    return nothing
end

## RHS

function build_mms_rhs!(F, M, K, uvec, utvec)
    mul!(F, M, utvec)   # F = M*u_t
    F .+= K * uvec       # F = F + K*u
end

## Interpolation to DOFs
function interpolate_to_dofs!(vec, dh, ufun, t)
    fill!(vec, 0)
    for cell in CellIterator(dh)
        cds = celldofs(cell)
        coords = getcoordinates(cell)
        for (i, dof) in enumerate(cds)
            vec[dof] = ufun(coords[i], t)
        end
    end
end

## Backward Euler

function solve_mms(dh, cv; k=1.0, Δt=0.05, T=0.5)
    n_dofs = ndofs(dh)
    M = allocate_matrix(dh)
    K = allocate_matrix(dh)

    assemble_system!(M, K, cv, dh, k)

    A = Symmetric(M + Δt .* K)
    cholA = cholesky(A)

    u  = zeros(n_dofs)
    uₙ = zeros(n_dofs)

    u_exact  = zeros(n_dofs)
    ut_exact = zeros(n_dofs)
    F = zeros(n_dofs)

    interpolate_to_dofs!(uₙ, dh, u_ms, 0.0)
    u .= uₙ

    t = 0.0
    nsteps = round(Int, T/Δt)

    for n in 1:nsteps
        t += Δt
        interpolate_to_dofs!(u_exact, dh, u_ms, t)
        interpolate_to_dofs!(ut_exact, dh, ut_ms, t)
        build_mms_rhs!(F, M, K, u_exact, ut_exact)
        rhs = M * uₙ + Δt .* F
        u .= cholA \ rhs
        uₙ .= u
    end
    return u, M
end

## Error

function compute_L2_error(u, M, dh, t)
    u_ex = zeros(ndofs(dh))
    interpolate_to_dofs!(u_ex, dh, u_ms, t)
    e = u .- u_ex
    return sqrt(real(e' * M * e))
end

## Final trial
include("02_01_createSphere.jl")

function run_convergence_test(final_ref::Int=2)
    ip = Lagrange{RefTriangle, 1}()
    qr = QuadratureRule{RefTriangle}(2)

    for ref in 0:final_ref
        Ω = create_embedded_sphere(ref)
        dh = DofHandler(Ω)
        add!(dh, :u, ip)
        close!(dh)

        cv = CellValues(qr, ip, ip^3)

        for Δt in (0.1, 0.025)
            u, M = solve_mms(dh, cv; k=1.0, Δt=Δt, T=0.5)
            err = compute_L2_error(u, M, dh, 0.5)
            @printf("%-6d %-8.3f %-12d %-15.6e\n",
                    ref, Δt, ndofs(dh), err)
        end
    end
end

run_convergence_test(3)

## EOC sweep over mesh refinement at fixed Δt
function run_test(ref::Int=1; k_g=1.0, Δtg::Float64=0.05, T_g::Float64=0.5)
    Ω = create_embedded_sphere(ref)
    ip = Lagrange{RefTriangle, 1}()
    qr = QuadratureRule{RefTriangle}(2)

    dh = DofHandler(Ω)
    add!(dh, :u, ip)
    close!(dh)

    cellvalues = CellValues(qr, ip, ip^3)

    u, M = solve_mms(dh, cellvalues; k=k_g, Δt=Δtg, T=T_g)
    err = compute_L2_error(u, M, dh, T_g)
    @printf("%-6d %-8.3f %-12d %-15.6e\n",
            ref, Δtg, ndofs(dh), err)
    return ref, Δtg, ndofs(dh), err
end

Refs = 0:1:3
Δts  = 0.05
Ts   = 0.5

n_dofs_coll = zeros(length(Refs))
errs        = zeros(length(Refs))

for refs in Refs
    results = run_test(refs; k_g=1.0, Δtg=Δts, T_g=Ts)
    n_dofs_coll[refs+1] = results[3]
    errs[refs+1]        = results[4]
end

function EOC(e1, e2, h1, h2)
    return log(e2/e1) / log(h2/h1)
end

# h ~ n_dofs^(-1/2)
h_approx(nd) = nd^(-0.5)
Hs = h_approx.(n_dofs_coll)

EOCerros = zeros(length(Refs))
for i in 2:length(Refs)
    EOCerros[i] = EOC(errs[i-1], errs[i], Hs[i-1], Hs[i])
end

println("\nRefinement  Δt      n_dofs   L2 error        EOC")
for i in 1:length(Refs)
    @printf("%-10d  %-7.3f %-8d %-15.6e %s\n",
            Refs[i], Δts, n_dofs_coll[i], errs[i],
            i == 1 ? "  --" : @sprintf("%.3f", EOCerros[i]))
end
