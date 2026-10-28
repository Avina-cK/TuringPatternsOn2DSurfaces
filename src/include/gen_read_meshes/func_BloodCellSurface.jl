
P(d,b,h) = -(d*d)/2 + (h*h/2)*(((d^2)/(b^2)) -1.0)*(1.0 - sqrt(1 - ((b^2)/(h^2))))
Q(d,b,h) = P(d,b,h)*((d^2)/(b^2)) + ((b*b/4)*(((d/b)^4) - 1.0))
R(d,b,h) = -(d*d/4)*(P(d,b,h) + (d*d/4))

"""
    blood_cell(x,y,z; d=4.0,b=0.8,h=1.0) -> dₛ::Float64
        dₛ(x,y,z) = (x² + y² + z²)² + P(x² + y²) + Qz² + R
    where
        0<d and 0<b<h
        P(d,b,h) = -(d*d)/2 + (h*h/2)*(((d^2)/(b^2)) -1.0)*(1.0 - sqrt(1 - ((b^2)/(h^2))))
        Q(d,b,h) = P(d,b,h)*((d^2)/(b^2)) + ((b*b/4)*(((d/b)^4) - 1.0))
        R(d,b,h) = -(d*d/4)*(P(d,b,h) + (d*d/4))
"""
function blood_cell(x,y,z; d=4.0,b=0.8,h=1.0)
    @assert h>=b "max. height, h, must not be less than b"
    P_bc = P(d,b,h)
    Q_bc = Q(d,b,h)
    R_bc = R(d,b,h)
    return (x^2 + y^2 + z^2)^2 + P_bc*(x^2 + y^2) + Q_bc*z*z + R_bc
end

"""
    grad_blood_cell(x,y,z; d=4.0,b=0.8,h=1.0) -> (∇d₁, ∇d₂, ∇d₃)::Tuple{Float, Float, Float}
        ∇d₁ = 4(x² + y² + z²)x + 2Px
        ∇d₂ = 4(x² + y² + z²)y + 2Py
        ∇d₃ = 4(x² + y² + z²)z + 2Qz
    where
        0<d and 0<b<h
        P(d,b,h) = -(d*d)/2 + (h*h/2)*(((d^2)/(b^2)) -1.0)*(1.0 - sqrt(1 - ((b^2)/(h^2))))
        Q(d,b,h) = P(d,b,h)*((d^2)/(b^2)) + ((b*b/4)*(((d/b)^4) - 1.0))
"""
function grad_blood_cell(x,y,z; d=4.0,b=0.8,h=1.0)
    @assert h>=b "max. height, h, must not be less than b"
    P_bc = P(d,b,h)
    Q_bc = Q(d,b,h)
    ρ = x^2 + y^2 + z^2
    gx = 4.0*ρ*x + 2.0*P_bc*x
    gy = 4.0*ρ*y + 2.0*P_bc*y
    gz = 4.0*ρ*z + 2.0*Q_bc*z
    return gx,gy,gz
end

@info "Loaded functions blood_cell(x,y,z; d=4.0,b=0.8,h=1.0) and grad_blood_cell(x,y,z; d=4.0,b=0.8,h=1.0)"