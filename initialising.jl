using Pkg
Pkg.activate(".")
Pkg.add([
    "Ferrite",
    "FerriteGmsh",
    "Gmsh",
    "LinearAlgebra",
    "NaNMath",
    "Plots",
    "RuntimeGeneratedFunctions",
    "Symbolics"
    ])

read("Project.toml", String)