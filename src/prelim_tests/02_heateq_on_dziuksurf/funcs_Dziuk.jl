using LinearAlgebra

# ---- normal on surface: n(x) --------------------------- #
#                ((x₁ - x₃²), (x₂), (x₃(1 - 2(x₁ - x₃²)))) 
#  n(x₁,x₂,x₃) = —————————————————————————————————————————
#                       (1 + 4x₃²(1 - x₁ - x₂²))²

# helper denominator function
function η(vp)
    x1 = vp[1];
    x2 = vp[2];
    x3 = vp[3];
    val = 1.0 + (4.0 * x3 * x3 * (1.0 - x1 - (x2*x2)));
    return val
end

"""
n1(x): 1st element of the normal vector
"""
function n1(vp)
    x1 = vp[1];
    x2 = vp[2];
    x3 = vp[3];
    η_val = η(vp)
    scale = 1.0 / (sqrt(η_val));
    return scale * (x1 - (x3*x3));
end

"""
n2(x): 2nd element of the normal vector
"""
function n2(vp)
    x1 = vp[1];
    x2 = vp[2];
    x3 = vp[3];
    η_val = η(vp)
    scale = 1.0 / (sqrt(η_val));
    return scale * x2;
end

"""
n3(x): 3rd element of the normal vector
"""
function n3(vp)
    x1 = vp[1];
    x2 = vp[2];
    x3 = vp[3];
    η_val = η(vp)
    scale = 1.0 / (sqrt(η_val));
    return scale * x3 * (1.0 - 2.0*(x1 - (x3*x3)));
end

"""
n(x): normal vector
    n(x) = (n1(x), n2(x), n3(x))
"""
function n(vp)
    n₁ = n1(vp);
    n₂ = n2(vp);
    n₃ = n3(vp);
    final_n = vcat(n₁, n₂, n₃);
    return final_n
end

"""
∂₁n₁(x) 
        = ∂n₁/∂x₁ (at x)
        = (-2x₃²(x₁ + 2x₂² - 2)- 2x₃⁴+1)/(η^(3/2))
"""
function ∂₁n₁(vp)
    x1 = vp[1];
    x2 = vp[2];
    x3 = vp[3];
    η_val = η(vp);
    denom = η_val^(-3/2);
    num11 = -2.0 * x3 * x3
    num12 = x1 + (2.0 * x2 * x2) + (x3 * x3) - 2.0
    num = (num11 * num12) + 1.0;
    return denom * num;
end

""" 
∂₂n₁(x) 
        = ∂n₁/∂x₂ (at x)
        = -(4x_{2}x_{3}^{2}(x_{3}^{2}-x_{1}))/(η^(3/2))
"""
function ∂₂n₁(vp)
    x1 = vp[1];
    x2 = vp[2];
    x3 = vp[3];
    η_val = η(vp);
    denom = η_val^(-3/2);
    num1 = 4.0 * x2 * x3 * x3;
    num2 = (x3*x3) - x1;
    num = -1.0 * num1 * num2;
    return num *denom;
end

"""
∂₃n₁(x) 
        = ∂n₁/∂x₃ (at x)
        = 2x₃( 2x₁² -1 + 2x₁(x₂² + x₃² - 1) + 2x₃²(x₂² -1))/(η^(3/2))
"""
function ∂₃n₁(vp)
    x1 = vp[1];
    x2 = vp[2];
    x3 = vp[3];
    η_val = η(vp);
    denom = η_val^(-3/2);
    num1 = 2.0*x3;
    num2 = (2.0 * x1*x1) - 1.0;
    num3 = (2.0*x1)*((x2*x2) + (x3*x3) -1.0);
    num4 = (2.0*x3*x3)*((x2*x2)-1.0);
    num = num1*(num2 + num3 + num4);
    return num * denom;
end

"""
grad_n₁(x) 
"""
function grad_n1(vp)
    gn1 = ∂₁n₁(vp);
    gn2 = ∂₂n₁(vp);
    gn3 = ∂₃n₁(vp);
    return vcat(gn1, gn2, gn3);
end

"""
∂₁n₂(x) 
        = ∂n₂/∂x₁ (at x)
        = (2x₂x₃²)/(η^(3/2))
"""
function ∂₁n₂(vp)
    x1 = vp[1];
    x2 = vp[2];
    x3 = vp[3];
    η_val = η(vp);
    denom = η_val^(-3/2);
    num = 2.0 * x2 * x3 *x3;
    return num * denom;
end

"""
∂₂n₂(x) 
        = ∂n₂/∂x₂ (at x)
        = (1 - 4x₃²(x₁ - 1))/(η^(3/2))
"""
function ∂₂n₂(vp)
    x1 = vp[1];
    x2 = vp[2];
    x3 = vp[3];
    η_val = η(vp);
    denom = η_val^(-3/2);
    num = 1.0 - ((4.0 * x3 *x3)*(x1 - 1.0));
    return num * denom;
end

"""
∂₃n₂(x) 
        = ∂n₂/∂x₃ (at x)
        = (4x₂x₃(x₁ + x₂² - 1))/(η^(3/2))
"""
function ∂₃n₂(vp)
    x1 = vp[1];
    x2 = vp[2];
    x3 = vp[3];
    η_val = η(vp);
    denom = η_val^(-3/2);
    num1 = 4.0 * x3 *x2;
    num2 = x1 + (x2*x2) - 1.0;
    num = num1*num2;
    return num * denom;
end

"""
grad_n₂(x) 
"""
function grad_n2(vp)
    gn1 = ∂₁n₂(vp);
    gn2 = ∂₂n₂(vp);
    gn3 = ∂₃n₂(vp);
    return vcat(gn1, gn2, gn3);
end


"""
∂₁n₃(x) 
        = ∂n₃/∂x₁ (at x)
        = (x₃³(4x₁ + 8x₂² - 6) + 4x₃⁵ - 2x₃ ) /(η^(3/2))
"""
function ∂₁n₃(vp)
    x1 = vp[1];
    x2 = vp[2];
    x3 = vp[3];
    η_val = η(vp);
    denom = η_val^(-3/2);
    num11 = (x3)^3; 
    num12 = (4.0*x1) + (8.0*x2*x2) - 6.0;
    num1 = num11*num12;
    num2 = 4.0 * ((x3)^5);
    num3 = -2.0 * x3;
    num = num1 + num2 + num3;
    return num * denom;
end

"""
∂₂n₃(x) 
        = ∂n₃/∂x₂ (at x)
        = ( ) /(η^(3/2))
"""
function ∂₂n₃(vp)
    x1 = vp[1];
    x2 = vp[2];
    x3 = vp[3];
    η_val = η(vp);
    denom = η_val^(-3/2);
    num1 = 4.0 * x2 * (x3^3);
    num2 = 1.0 + (2.0*x3*x3) - (2.0*x1);
    num = num1*num2;
    return num * denom;
end

"""
∂₃n₃(x) 
        = ∂n₃/∂x₃ (at x)
        = ( ) /(η^(3/2))
"""
function ∂₃n₃(vp)
    x1 = vp[1];
    x2 = vp[2];
    x3 = vp[3];
    η_val = η(vp);
    denom = η_val^(-3/2);
    num1 = -2.0*x1*(1.0 + (8.0 * (x3^4)));
    num2 = -16.0 * (x3^4) * ((x2*x2) - 1);
    num3 = 6.0*x3*x3;
    num = num1 + num2 + num3 +1.0;
    return num * denom;
end

"""
grad_n₃(x) 
"""
function grad_n3(vp)
    gn1 = ∂₁n₃(vp);
    gn2 = ∂₂n₃(vp);
    gn3 = ∂₃n₃(vp);
    return vcat(gn1, gn2, gn3);
end

"""
div_n(x)
"""
function div_n(vp)
    d1 = ∂₁n₁(vp);
    d2 = ∂₂n₂(vp);
    d3 = ∂₃n₃(vp);
    return d1 + d2 + d3
end

# ---- end of normal & functions ---- #

# Problem to solve: -Δₛu = f on S
# True solution: u(x)=x₁x₂
function true_u(vp)
    return vp[1]*vp[2]
end
function ∇u_true(vp)
    ∇u₁ = vp[2];
    ∇u₂ = vp[1];
    ∇u₃ = 0.0;
    final_∇u_true = [∇u₁,∇u₂,∇u₃];
    return final_∇u_true
end

# -------------------------------------
#= RHS: f = -∇ₛ⋅v = -∇⋅v + Σ₁³(∇vⱼ⋅n)nⱼ
    where v = ∇u - (∇u ⋅ n)n
=#

function v_f(vp)
    gradu = ∇u_true(vp);
    normal = n(vp);
    p21 = dot(gradu, normal);
    p2 = p21 * normal;
    final_v_f = gradu - p2;
    return final_v_f
end

# ∇vⱼ = 
function grad_v_f_j(vp, j::Int)
    x1 = vp[1];
    x2 = vp[2];
    x3 = vp[3]; 
    n₁ = n1(vp);
    n₂ = n2(vp);
    n₃ = n3(vp);
    grad_n₁ = grad_n1(vp);
    grad_n₂ = grad_n2(vp);
    grad_n₃ = grad_n3(vp);
    
    if j==1
        p1 = vcat(-n₁*n₂, 1.0 - (n₁*n₁), 0.0);
        s1 = (2.0 * x2 * n₁) + (x1 * n₂);
        s2 = x1*n₁
        final_grad_v_f_j = p1 - (s1 * grad_n₁) - (s2 * grad_n₂)
    elseif  j==2
        p1 = vcat(1.0 - (n₂*n₂),-1.0*n₁*n₂,  0.0);
        s1 = n₂*x2;
        s2 = (n₁*x2) + (2.0*n₂*x1);
        final_grad_v_f_j = p1 - (s1 * grad_n₁) - (s2 * grad_n₂)
    elseif j==3
        p1 = vcat(-1.0*n₂*n₃, -1.0*n₁*n₃, 0.0);
        s1 = x2*n₃;
        s2 = x1*n₃;
        s3 = (x2*n₁)+(x1*n₂);
        final_grad_v_f_j = p1 - (s1 * grad_n₁) - (s2 * grad_n₂) - (s3*grad_n₃)
    else print("j must belong to {1,2,3}")
    end
    return final_grad_v_f_j
end

# ∇⋅v = -(∇(∇u ⋅ n)⋅n) - (∇u ⋅ n)(∇⋅n)
function div_v_f(vp)
    gradu = ∇u_true(vp);
    normal = n(vp);
    graddot = vp[1]*grad_n2(vp) + vp[2]*grad_n1(vp) + vcat(n2(vp), n1(vp), 0.0);
    p1 = dot(graddot, normal);
    p2 = dot(gradu, normal) * div_n(vp);
    return -1.0 *(p1 + p2)
end

#f = -∇⋅v + Σ₁³(∇vⱼ⋅n)nⱼ
function rhs_f(vp)
    div_v = div_v_f(vp);
    normal = n(vp);
    sum1 = dot(grad_v_f_j(vp,1), normal) * n1(vp);
    sum2 = dot(grad_v_f_j(vp,2), normal) * n2(vp);
    sum3 = dot(grad_v_f_j(vp,3), normal) * n3(vp);
    sum_f = sum1+sum2+sum3;
    return sum_f - div_v
end
