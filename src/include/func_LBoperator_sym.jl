using Symbolics

# Define the variables
@variables x₁ x₂ x₃ tᵥ

"""
∇ₛ(u::Function, d::Function)
    - u(x,y,z) is the function for which we want to compute the Laplace-Beltrami operator
    - d(x,y,z) is the level set function that defines the surface (d=0 on the surface)
    The function returns a callable that computes the Laplace-Beltrami operator at any point (x,y,z).
The Laplace-Beltrami operator, on a surface, d=0, is given by:
    ∇ₛ(u) = ∇u - (∇u ⋅ n) n, where n = ∇d / |∇d| is the unit normal vector to the surface.
Example:
    d(x,y,z) = z^2 + y^2 + (x - z^2)^2 - 1.0
    u(x,y,z) = x*y
    ∇ₛ_u_given = ∇ₛ(u,d)
    ∇ₛ_u_given(0.0, -0.50, -1.0)
Then you can evaluate ∇ₛ_u_given at any point (x,y,z) to get the value of the Laplace-Beltrami operator at that point.
"""
function ∇ₛ(u::Function, d::Function)
    u_sym = u(x₁, x₂, x₃);
    d_sym = d(x₁, x₂, x₃);
    gradu = Symbolics.gradient(u_sym, [x₁, x₂, x₃]);
    gradd = Symbolics.gradient(d_sym, [x₁, x₂, x₃]);
    norm_n = sqrt(sum(gradd.^2));
    n = gradd ./ norm_n;
    LB_sym = gradu .- (sum(gradu .* n) .* n);
    LB_func, _ = Symbolics.build_function(LB_sym, x₁, x₂, x₃; expression = Val{false})
    return (x,y,z) -> LB_func(x,y,z)
end
