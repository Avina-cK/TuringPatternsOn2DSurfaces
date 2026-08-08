using Gmsh.gmsh

## Initial mesh vertices and triangles
# initial surface has 6 points:
vertices0 = [1.25, 0.0, -sqrt(3)/2, #1 : bottom pole
             1.0, 0.0, 0.0,         #2 : equator point 1
             1.25, 0.0, sqrt(3)/2,  #3 : top pole
             0.0, 1.0, 0.0,         #4 : equator point 2
             -1.0, 0.0, 0.0,        #5 : equator point 3
             0.0, -1.0, 0.0         #6 : equator point 4
            ]
vertices0_numbering = collect(1:(length(vertices0)÷3))

# Initial surface has 8 triangles:
triangles0 = [1, 4, 2,
              1, 5, 4,
              1, 6, 5,
              1, 2, 6,
              3, 2, 4,
              3, 4, 5,
              3, 5, 6,
              3, 6, 2]
triangles0_numbering = collect(1:(length(triangles0)÷3))

## Create mesh using Gmsh
# initialise gmsh
gmsh.initialize(append!(["gmsh"], ARGS))

# add a model
gmsh.model.add("Dzuik_mesh_0")

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
gmsh.write("../../Dziuk_surf_meshes/Dzuik_mesh_0.msh")

# finalise and close gmsh
gmsh.finalize()