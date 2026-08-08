## Helper function
include("func_mesh_size.jl")
##
using Ferrite, FerriteGmsh
function create_embedded_sphere(refinements; tosave=false)
    gmsh.initialize()

    gmsh.model.occ.addSphere(0.0, 0.0, 0.0, 1.0)
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(2)

    for _ in 1:refinements
        gmsh.model.mesh.refine()
    end

    if tosave
        gmsh.write("$(dir_to_spheres)/embedded_sphere_$(refinements).msh")
    end
    
    nodes = tonodes()
    elements, _ = toelements(2)
    Ω = Grid(elements, nodes)

    max_edge_length = mesh_size()

    gmsh.finalize()
    return Ω, max_edge_length
end

function load_sphere(refinements, dir_to_spheres="../../sphere_meshes")
    gmsh.initialize()
    gmsh.open("$(dir_to_spheres)/embedded_sphere_$(refinements).msh")
    nodes = tonodes()
    elements, _ = toelements(2)
    Ω = Grid(elements, nodes)
    max_edge_length = mesh_size()
    gmsh.finalize()
    return Ω, max_edge_length
end

## Main reference: ferrite-fem.github.io/Ferrite.jl/stable/tutorials/reactive_surface/