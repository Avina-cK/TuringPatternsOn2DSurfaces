
# Reaction-Diffusion equation (RDE) problem
Given a domain $\Omega$, $k\in\mathbb{R}$ and two reaction functions $f,g$, find $(u,v):\left( \Omega, [0,\infty) \right)\rightarrow \mathbb{R} \times \mathbb{R}$ such that
```math
\frac{\partial}{\partial t} u(x,t) - D\_u \Delta_{\Omega} \left( u(x,t) \right) = f(u,v,x,t;p) \quad \forall\ x \in \Omega ,\ t\in[0,\infty)
```
```math
\frac{\partial}{\partial t} v(x,t) - D\_v \Delta_{\Omega} \left( v(x,t) \right) = g(u,v,x,t;p) \quad \forall\ x \in \Omega ,\ t\in[0,\infty)
```
where 
- $\Delta_{\Omega}$ is the Laplace-Beltrami operator
- $D_u $ and $D_v $ are positive diffusion coefficients
- $p$ is a set of constant parameters

The reaction functions are often only implicitly dependent on $(x,t)$ via $(u,v)$, i.e., $f(u,v,x,t;p) \equiv f(u,v;p)$ and $g(u,v,x,t;p) \equiv g(u,v;p)$.

## Gray-Scott system
```math
f(u,v; \alpha , \beta) = \alpha (1-u) - uv^2
```
```math
g(u,v; \alpha , \beta) =  -(\alpha +\beta)v + uv^2
```

# Domain(s)

## Dziuk surface
