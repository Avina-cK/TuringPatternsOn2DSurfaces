# Brief overview of file contents

## funcs_error_analysis.jl
Defines three functions to calculate $L^2$-errors on a grid, $\Omega$:
1. `L2_norm(u_given::Vector, cellvalues::CellValues, dh::DofHandler)`: Computes the discrete version of the square of the $L^2$ norm of a function $u$ (`u_given`) as 
```math
∥u∥_{L^2}​^2= \int_{\Omega}​u^2 d\Omega
```
2. `compute_L2_error(u_h::Vector, u_exact::Vector, cellvalues::CellValues, dh::DofHandler)`:
Computes the discrete version of the square of the $L^2$ error of a function, $u_h$ (`u_h`), with respect to another function, $u_{exact}$ (`u_exact`) as 
```math
∥u_{exact} - u_h∥_{L^2}^2 = \int_{\Omega}​{u_{exact} - u_h}^2 d\Omega 
```
4. `compute_relative_error(u_h::Vector, u_exact::Vector, cellvalues::CellValues, dh::DofHandler)`: Computes the discrete relative error of a function, $u_h$ (`u_h`), with respect to another function, $u_{exact}$ (`u_exact`) as 
```math
\frac{∥u_{exact} - u_h∥_{L^2}^2}{∥u_{exact}∥_{L^2}^2}
```
5. `EOC(h₁, h₂, L2err1, L2err2)`: Computes the experimental order of convergence between two meshes, of size $h_1$ and $h_2$ as
```math
     EOC\left(h_1, h_2 \right) =\frac{ \log \left( {\frac{∥u_{exact} - u_{h_1}∥_{L^2}^2}{∥u_{exact} - u_{h_2}∥_{L^2}^2}}\right)} {\log \left({\frac{h_1}{h_2}} \right) }
```

## funcs_gensurface.jl
Defines six functions used to generate the Dziuk implicit surface mesh, given by $\left{ (x_1, x_2, x_3):d_S(x_1, x_2, x_3)=0 \right}$.
1. `Dziuk_surface(x,y,z)`
```math
d_S (x,y,z) = z^2 + y^2 + (x - z^2)^2 - 1.0
```
2. `gradDziuk_surface(x,y,z)`: Computes the gradient of $d_S$
```math
\nabla d_S (x,y,z) = 2.0 × \left[ x-z^2, y, z - 2z \left( x - z^2 \right) \right]
```
3. `project_to_surf(x,y,z)`: Projects a point $(x,y,z)$ to a point on the Dziuk surface by running the following mapping iteratively
```math
(x,y,z) \mapsto (x,y,z) - \frac{d_S (x,y,z)}{\nabla d_S (x,y,z)\cdot\nabla d_S (x,y,z)}
```
4. `break_direct_edge!(vertices, triangles, i1, i2)`: If vertex indices `i1` and `i2` share one or more triangle edges, this function splits each such triangle by inserting a midpoint vertex (projected back onto the implicit surface).
5. `edge_length_stats(vertices::Matrix{Float64}, triangles::Matrix{<:Integer}) -> (min_len, max_len)`: Calculates the minimum and maximum edge length of a mesh.
6. `surface_lloyd(coords, triangles; iterations, fixed_indices)`: Implements LLoyd's algorithm on a mesh for a given number of iterations and projects the new vertices back to the Dzuik surface (in each iteration). If `fixed indices` contains the indices of vertices that should not be moved, then the function ensures that. The Lloyd's algorithm basically calculates the area, $A_i$ , and centroid, $c_i$ , of all triangles, and then uses the following formula to calculate the new coordinates of each vertex, $v_i$:
```math
v_i^{new} = \frac{1}{\sum_j{A_j}} \cdot \sum_j {A_j \cdot c_j}
```
where $j$ is runs over all triangles touching $v_i$.
