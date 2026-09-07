using CSV, DataFrames
cd(@__DIR__)
results = CSV.read("log.csv", DataFrame)

αs = results[:," α"]
βs = results[:," β"]
patterns = results[:," pattern"]
possible_patterns = [
    "const.",
    "swirls",
    "chaos",
    "periodic",
    "spots",
    "spots + blob-strip",
    "spots + blob-strips",
    "spots + smooth-strips",
    "smooth-strips",
]

pattern_map = Dict(p => i for (i, p) in enumerate(possible_patterns))
patterns_no = [pattern_map[p] for p in patterns]

# Split into per-pattern vectors
alphas = [αs[patterns_no .== i] for i in 1:length(possible_patterns)]
betas  = [βs[patterns_no .== i] for i in 1:length(possible_patterns)]

# Pick 9 distinguishable shapes and colors
shapes = [:rect, :star5, :star8, :pentagon, :circle, :utriangle, :rtriangle, :ltriangle, :dtriangle]
colors = [:white, :yellow, :red, :green, :black, :orange, :blue, :purple, :pink]

include("../../../include/plot_defaults.jl")
plt = Plots.scatter(
     betas[1],alphas[1],
    label=possible_patterns[1],
    markershape=shapes[1],
    color=colors[1],
    ylabel="α", xlabel="β",
    #ylim = (0.0435,0.0615),
    #xlim = (0.005, 0.065)
)
patterns_len = length(possible_patterns)
for i in 2:patterns_len
    Plots.scatter!(
        plt,
        betas[i],alphas[i], 
        label=possible_patterns[i],
        markershape=shapes[i],
        color=colors[i],
        markersize=8
    )
end

plot!(legend=:outerright,
    title="Gray-Scott on Dziuk surface (Du,Dv)=(1.6, 0.8)e-4")

#savefig("GS_Dziuk_surf_initresults.png")


## Pearson (and Munafo) Classes
# http://www.mrob.com/pub/comp/xmorphia/pearson-classes.html#alp

# incomplete list
PM_Classes = [
    "R"   "uniform high u";
    "B"   "uniform low u";
    "α"   "chaotic swirls";
    "β"   "chaos";
    "θ"   "stripes (labyrinth)";
    "δ"   "spots tessalation";
    "η"   "spots + stripes";
    "ι"   "low u + spots";
    "μ"   "stripes (worms) + spots";
    "ξ"   "sustained swirls";
    "π"   "low u + patterns";
    "ρ"   "soap bubbles"  
]

PM_Classification = DataFrame([:type => PM_Classes[:,1], :description => PM_Classes[:,2]])
