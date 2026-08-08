using Plots
 
"""
    export_colourmap_xml(grad::ColorGradient, filename::String;
                         name::String = "gradient",
                         space::String = "CIELAB",
                         nan_rgb::Tuple{Float64,Float64,Float64} = (0.25, 0.0, 0.0))
 
Write a Plots.jl/PlotUtils `cgrad` ColorGradient out to an XML file using the
<ColorMaps><ColorMap>...<Point .../></ColorMap></ColorMaps> schema
(e.g. as used by ParaView-style colormap files).
 
One <Point> is written for each color stop in the gradient (`grad.values`
gives the x-position of each stop, `grad.colors` gives the corresponding
color). Opacity is always written as 1.
"""
function export_colourmap_xml(grad, filename::String;
                              name::String = "gradient",
                              space::String = "CIELAB",
                              nan_rgb::Tuple{Float64,Float64,Float64} = (0.25, 0.0, 0.0))
    xs = grad.values    # stop positions, 0..1
    cs = grad.colors.colors # RGB colours at each stop
 
    open(filename, "w") do io
        println(io, "<ColorMaps>")
        println(io, "<ColorMap space=\"$space\" indexedLookup=\"false\" name=\"$name\">")
        for (x, c) in zip(xs, cs)
            rc = red(c); gc = green(c); bc = blue(c)
            println(io, "<Point x=\"$(round(x, digits=6))\" o=\"1\" " *
                        "r=\"$(round(rc, digits=6))\" " *
                        "g=\"$(round(gc, digits=6))\" " *
                        "b=\"$(round(bc, digits=6))\"/>")
        end
        nr, ng, nb = nan_rgb
        println(io, "<NaN r=\"$nr\" g=\"$ng\" b=\"$nb\"/>")
        println(io, "</ColorMap>")
        println(io, "</ColorMaps>")
    end
    return filename
end
 
"""
    export_colormaps_xml(grads::Dict{String,<:Any}, outdir::String = "."; kwargs...)
 
Batch version: loop over a name => gradient dict and write one XML file per
gradient, named "<key>.xml", into `outdir`. Extra `kwargs` (space, nan_rgb)
are forwarded to `export_colourmap_xml` for every gradient.
"""
function export_colormaps_xml(grads::Dict{String,<:Any}, outdir::String = "."; kwargs...)
    isdir(outdir) || mkpath(outdir)
    for (name, grad) in grads
        filename = joinpath(outdir, "$(name).xml")
        export_colourmap_xml(grad, filename; name = name, kwargs...)
        println("wrote $filename")
    end
end

## 
include("QColourGradients.jl")

all_grads = Dict(
    "gm7" => gm7,
    "les7" => les7,
    "nb4_bpyw" => nb4_bpyw,
    "rainbow6" => rainbow6,
    "tragen3" => tragen3,
    "bis3" => bis3,
    "ace4" => ace4,
    "pan3" => pan3,
    "aro5" => aro5,
    "aroace5" => aroace5,
    "genque3" => genque3
    )
export_colormaps_xml(all_grads, "colormaps")