using LinearAlgebra
using Roots
using Gmsh.gmsh

include("func_BloodCellSurface.jl")
include("../func_floattostring.jl")
## 

"""
    gen_init_bloodcell_mesh(d_in, b_in, h_in)

# Examples
```julia-repl
julia> d_in = 4.0; h_in = 1.0; b_in = 0.8
julia> gen_init_bloodcell_mesh(d_in, b_in, h_in)
Info    : Writing '../blood_cell_meshes/Blood_cell_mesh_4-00_0-80_1-00_0.msh'...
Info    : Done writing '../blood_cell_meshes/Blood_cell_mesh_4-00_0-80_1-00_0.msh'
[ Info: Blood cell mesh generated for d=4.0, b=0.8, h=1.0.
```
"""
function gen_init_bloodcell_mesh(d_in, b_in, h_in)
    @assert h_in>=b_in "max. height, h, must not be less than dimple height, b."

    p0 = 0.5*b_in
    p2 = d_in/(2*√2)

    if h_in>b_in
        pt = 0.5*h_in
        x_v_func(x) = blood_cell(x,0.0,h_in/2; d=d_in,b=b_in,h=h_in)
        p1 = Roots.find_zero(x_v_func, d_in/4)
    elseif h_in==b_in
        p1 = 0.25*d_in
        z_v_func(z) = blood_cell(0.25*d_in,0.0,z; d=d_in,b=b_in,h=h_in)
        pt = Roots.find_zero(z_v_func, (0.45*h_in))
    end
    vertices0 = [
        0.0,    0.0,    p0,   #1
        0.0,    0.0,   -p0,   #2
        p1,     0.0,    pt,   #3
        p1,     0.0,   -pt,   #4
        -p1,    0.0,   -pt,   #5
        -p1,    0.0,    pt,   #6
        0.0,    p1,     pt,   #7
        0.0,    -p1,    pt,   #8
        0.0,    p1,    -pt,   #9
        0.0,    -p1,   -pt,   #10
        p2,     p2,     0.0,  #11
        p2,     -p2,    0.0,  #12
        -p2,    p2,     0.0,  #13
        -p2,    -p2,    0.0,  #14
    ]
    vertices0_numbering = collect(1:(length(vertices0)÷3))

    triangles0 = [
        1, 3, 7,
        1, 7, 6,
        1, 6, 8,
        1, 8, 3,
        3, 8, 12,
        3, 12, 11,
        3, 11, 7,
        7, 13, 6,
        7, 11, 13,
        6, 13, 14,
        6, 14, 8,
        8, 14, 12,
        4, 11, 12,
        4, 12, 10,
        10, 12, 14,
        10, 14, 5,
        5, 14, 13,
        5, 13, 9,
        9, 13, 11,
        9, 11, 4,
        2, 4, 10, 
        2, 10, 5,
        2, 5, 9,
        2, 9, 4
    ]
    triangles0_numbering = collect(1:(length(triangles0)÷3))

    gmsh.initialize(append!(["gmsh"], ARGS))

    # add a model
    gmsh.model.add("Blood_cell_mesh_$(fts(b_in, 2))")

    # add a discrete entity of dimension 2 and tag 1
    gmsh.model.addDiscreteEntity(2, 1)

    # add nodes to the model
    gmsh.model.mesh.addNodes(2, 1, vertices0_numbering, vertices0)

    # add triangles to the model
    gmsh.model.mesh.addElements(2, 1,
                                [2],    # elementType 2 ↦ 3-node triangle
                                [triangles0_numbering], # elementTags
                                [triangles0])   # nodeTags describing triangles

    # Write mesh into file
    cd(@__DIR__)
    mesh_name="Blood_cell_mesh_$(fts(d_in, 2))_$(fts(b_in, 2))_$(fts(h_in, 2))_0"
    gmsh.write("../../blood_cell_meshes/$(mesh_name).msh")

    # finalise and close gmsh
    gmsh.finalize()
    @info "Blood cell mesh generated for d=$(d_in), b=$(b_in), h=$(h_in)."
end


gen_init_bloodcell_mesh(4.0, 1.0, 1.0)