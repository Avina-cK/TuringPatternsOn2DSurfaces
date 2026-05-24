using Symbolics

# Define the variables (do not use commas!)
@variables x₁ x₂ x₃

"""
normal_atpt(d::Function)
    - d(x,y,z) is the level set function that defines the surface (d=0 on the surface)
    The function returns a callable that computes the normal at any point (x,y,z).
        n = ∇d / |∇d|
Example:
    d(x,y,z) = z^2 + y^2 + (x - z^2)^2 - 1.0
    n_d = normal_atpt(d)
    n_d(0.0, -0.50, -1.0)
Then you can evaluate n_d at any point (x,y,z) to get the value of the normal at that point.
"""
function normal_atpt(d::Function)
    d_sym = d(x₁, x₂, x₃)
    gradd = Symbolics.gradient(d_sym, [x₁, x₂, x₃])
    norm_n = sqrt(sum(gradd.^2))
    normalfuncsym = gradd ./ norm_n
    normalfunc_out,_ = build_function(normalfuncsym, x₁, x₂, x₃;
                                 expression = Val{false}
                                 , simplify=Val{true}
                                 , target=Symbolics.JuliaTarget())
    return (x,y,z) -> normalfunc_out(x,y,z)
end

"""
∇ₛ(u::Function, d::Function)
    - u(x,y,z) is the function for which we want to compute the tangential gradient
    - d(x,y,z) is the level set function that defines the surface (d=0 on the surface)
    The function returns a callable that computes the tangential gradient at any point (x,y,z).
The tangential gradient, on a surface,d=0 , is given by:
    ∇ₛ(u) = ∇u - (∇u ⋅ n) n, where n = ∇d / |∇d| is the unit normal vector to the surface.
Example:
    d(x,y,z) = z^2 + y^2 + (x - z^2)^2 - 1.0
    u(x,y,z) = x*y
    ∇ₛ_u_given = ∇ₛ(u,d)
    ∇ₛ_u_given(0.0, -0.50, -1.0)
Then you can evaluate ∇ₛ_u_given at any point (x,y,z) to get the value of the tangential gradient at that point.
"""
function ∇ₛ(u::Function, d::Function)
    u_sym = u(x₁, x₂, x₃);
    d_sym = d(x₁, x₂, x₃);
    gradu = Symbolics.gradient(u_sym, [x₁, x₂, x₃]);
    gradd = Symbolics.gradient(d_sym, [x₁, x₂, x₃]);
    norm_n = sqrt(sum(gradd.^2));
    n = gradd ./ norm_n;
    TG_sym = gradu .- (dot(gradu, n) * n);
    TG_func, _ = build_function(TG_sym, x₁, x₂, x₃;
                                 expression = Val{false}
                                 , simplify=Val{true}
                                 , target=Symbolics.JuliaTarget())
    return (x,y,z) -> TG_func(x,y,z)
end

"""
Δₛ(u::Function, d::Function)
    - u(x,y,z) is the function for which we want to compute the tangential gradient
    - d(x,y,z) is the level set function that defines the surface (d=0 on the surface)
    The function returns a callable that computes the Laplace-Beltrami operator on u at any point (x,y,z).
The Laplace-Beltrami operator on u (on a surface, d=0) is given by:
    Δₛu = ∇ₛ ⋅ ∇ₛ(u)
Example:
    d(x,y,z) = z^2 + y^2 + (x - z^2)^2 - 1.0
    u(x,y,z) = x*y
    Δₛ_u_given = Δₛ(u, d)
    Δₛ_u_given(0.0, -0.50, -1.0)
Then you can evaluate ∇ₛ_u_given at any point (x,y,z) to get the value of the tangential gradient at that point.
"""
function Δₛ(u::Function, d::Function)
    u_sym = u(x₁, x₂, x₃)
    d_sym = d(x₁, x₂, x₃)
    gradu = Symbolics.gradient(u_sym, [x₁, x₂, x₃])
    gradd = Symbolics.gradient(d_sym, [x₁, x₂, x₃])
    norm_n = sqrt(sum(gradd.^2))
    n = gradd ./ norm_n
    TG_sym = gradu .- (sum(gradu .* n) .* n)

    # Compute jacobian of tangential gradient
    JTG_sym = Symbolics.jacobian(TG_sym, [x₁, x₂, x₃])
    Δₛ_sym = JTG_sym[1,1] + JTG_sym[2,2] + JTG_sym[3,3]  # trace = divergence

    Δₛ_func = build_function(Δₛ_sym, x₁, x₂, x₃;
                                 expression = Val{false}
                                 , simplify=Val{true}
                                 , target=Symbolics.JuliaTarget())
   return (x,y,z) -> Δₛ_func(x,y,z)
end

@info "Loaded symbolically defined functions: normal_atpt(d), ∇ₛ(u,d) and Δₛ(u,d)."