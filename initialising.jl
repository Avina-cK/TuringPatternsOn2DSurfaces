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
    "SparseArrays",
    "Symbolics",
    "Tensors"
    ])

read("Project.toml", String)