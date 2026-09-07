# Linear analysis of Gray-Scott system
using Plots
using LinearAlgebra
using MeshGrid
using LaTeXStrings
## Steady state solutions

γ(α, β) = sqrt(Complex(1.0 - (4.0*((α+β)^2)/α)))

w0 = [1.0; 0.0]
w1(α, β) = [0.5*(1.0+γ(α, β)); 0.5*(α/(α+β))*(1.0-γ(α, β))]
w2(α, β) = [0.5*(1.0-γ(α, β)); 0.5*(α/(α+β))*(1.0+γ(α, β))]

# α ∈ (0.0, ¼)
# β ∈ (0.0, ½(√α - 2α))

## Jacobian of the system

J_GS(u, v, α, β) = [-α-(v*v) -2.0*u*v;
                    v*v (2.0*u*v)-(α+β)]

## Plotting eigenvalues of Jacobian

function J0(a,b)
    u_s = w0[1]
    v_s = w0[2]
    Jacob = J_GS(u_s, v_s, a, b)
    return Jacob
end

function J1(a,b)
    u_s = w1(a,b)[1]
    v_s = w1(a,b)[2]
    Jacob = J_GS(u_s, v_s, a, b)
    return Jacob
end

function J2(a,b)
    u_s = w2(a,b)[1]
    v_s = w2(a,b)[2]
    Jacob = J_GS(u_s, v_s, a, b)
    return Jacob
end

## Helper functions
function classify_real_complex_patterns(Z)
    Z_prop = Array{Int}(undef, size(Z))

    for i in eachindex(Z)
        r = Z[i]
        r1, r2 = r[1], r[2]
        re1, re2 = real(r1), real(r2)
        im1, im2 = imag(r1), imag(r2)

        if im1 == 0 && im2 == 0
            # purely real
            if re1 > re2 && re2 > 0
                Z_prop[i] = 6
            elseif re2 > re1 && re1 > 0
                Z_prop[i] = 6
            elseif re2 < re1 && re2 < 0 && re1 < 0
                Z_prop[i] = -6
            elseif re1 < re2 && re1 < 0 && re2 < 0
                Z_prop[i] = -6
            elseif re2 < 0 && re1 > 0
                Z_prop[i] = 2
            elseif re1 < 0 && re2 > 0
                Z_prop[i] = 2
            elseif re1 == re2 && re1 > 0
                Z_prop[i] = 4
            elseif re1 == re2 && re1 < 0
                Z_prop[i] = -4
            else
                Z_prop[i] = -2 #r1=r2=0
            end
        else
            # complex numbers
            if re1 > re2 && re2 > 0
                Z_prop[i] = 5
            elseif re2 > re1 && re1 > 0
                Z_prop[i] = 5
            elseif re2 < re1 && re2 < 0 && re1 < 0
                Z_prop[i] = -5
            elseif re1 < re2 && re1 < 0 && re2 < 0
                Z_prop[i] = -5
            elseif re2 < 0 && re1 > 0
                Z_prop[i] = 1
            elseif re1 < 0 && re2 > 0
                Z_prop[i] = 1
            elseif re1 == re2 && re1 > 0
                Z_prop[i] = 3
            elseif re1 == re2 && re1 < 0
                Z_prop[i] = -3
            else
                Z_prop[i] = -1 # catch other complex cases
            end
        end
    end
    return Z_prop
end

function round_complex_matrix!(Z, precision::Float64=1.0e-4)
    for i in eachindex(Z)
        Z[i] = ComplexF64[
            complex(round(real(z), digits=8), round(imag(z), digits=4)) for z in Z[i]
        ]
    end
    return Z
end
function find_valid_integers(M::AbstractMatrix{<:Integer})
    allowed = Set(-6:-1) ∪ Set(1:6)   # {-6,...,-1,1,...,6}, excludes 0
    found = Set{Int}()
    for x in M
        if x in allowed
            push!(found, x)
        end
    end
    return sort(collect(found))
end


les7 = cgrad([colorant"rgb(213, 45, 0)"
                    , colorant"rgb(239, 118, 39)"
                    , colorant"rgb(255, 154, 86)"
                    , colorant"rgb(255,255,255)"
                    , colorant"rgb(209, 98, 164)"
                    , colorant"rgb(181, 86, 144)"
                    , colorant"rgb(163, 2, 98)"]
                , 12, categorical=true)

##
a_min = 0.00005;
a_max = 0.26;
b_min = 0.00005;
b_max = 0.067;
δgrid = 0.0001;

a_vals = a_min:δgrid:a_max;
b_vals = b_min:δgrid:b_max;
A, B = meshgrid(a_vals, b_vals);
Jacob = J2.(A,B)
Z_evals = eigvals.(Jacob)
Z_evals_r = round_complex_matrix!(Z_evals)
Z_prop = classify_real_complex_patterns(Z_evals_r)


## Plot
plot_title=L"\mathrm{Eigenvalues\ of}\ \mathbf{J}(\mathit{\mathbf{w}}^*_2) : \lambda_{1/2} = r_{1/2} \pm c_{1/2}i"
p_main = heatmap(b_vals,a_vals, transpose(Z_prop)
    , c=les7
    , title=plot_title
    , xlabel=L"b", ylabel=L"a"
    , colorbar=true
    #, clim=(-absmaxZ,absmaxZ)
    , clim=(-6.8, 6.8)
    , fontfamily = "Computer Modern"
    , aspectratio=b_max/a_max
    , xlim=(b_min, b_max)
    , ylim=(a_min, a_max)
    , size = (900, 800)
    , guidefont = (14, :black), tickfont = (12, :black)
    , titlefont = (16, :black)
    )
contour!(p_main, b_vals,a_vals, transpose(Z_prop);
    levels=[6], linewidth=1.0, color=:red)
contour!(p_main, b_vals,a_vals, transpose(Z_prop);
    levels=[-1], linewidth=1.0, color=:blue)
    
#savefig("src/illustrations/GS_DziukSurface_lin_analysis_eigvals_w2.png")

##
ticklabels = [L"r_1>r_2>0,"*"\n"*L"c_{1/2}=0",
             L"r_1=r_2>0,"*"\n"*L"c_{1/2}=0",
             L"r_{1/2}>0,"*"\n"*L"c_{1/2}\ne 0",
             L"r_{1}<0<r_{2},"*"\n"*L"c_{1/2}= 0",
             L"r_{1/2}=0,"*"\n"*L"c_{1/2}\ne 0",
             L"r_{1/2}<0,"*"\n"*L"c_{1/2}\ne 0",
             L"r_{1}=r_{2}<0,"*"\n"*L"c_{1/2} = 0",
             L"r_1<r_2<0,"*"\n"*L"c_{1/2}= 0",]
nbins = 9
clim_lo, clim_hi = -4.5, 4.5
tickpos = range(clim_lo, clim_hi, length=nbins+1)
tickcenters = (tickpos[1:end-1] .+ tickpos[2:end]) ./ 2
cbar_y = range(clim_lo, clim_hi, length=200)
cbar_z = reshape(cbar_y, :, 1)

p_cbar = heatmap(1:1, cbar_y, cbar_z
    , c = les7
    , clim = (clim_lo, clim_hi)
    , legend = false
    , colorbar = false
    , xticks = false
    , yticks = (tickcenters, ticklabels)
    , ymirror = true            # put labels on the right, like a normal colorbar
    , tickfont = (12, :black)
    , framestyle = :box
    , grid = false
    )

# --- combine into one figure ---
plt = plot(p_main, p_cbar
    , layout = grid(1, 2, widths=[0.95, 0.05])
    , size = (1000, 800)
    )

plt