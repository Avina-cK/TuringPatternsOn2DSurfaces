# Aim
The purpose of these preliminary tests is to solve specific partial differential equation (PDE) systems on specific domains, and ensure that the numerical method(s) work correctly.

## Laplace-Beltrami problem
The Laplace-Beltrami-Operator (LBO) problem is defined as follows: Given a function, $f:\Omega \rightarrow \mathbb{R}$, on a domain, $\Omega$, find $u:\Omega \rightarrow \mathbb{R}$ such that
```math
-\Delta_{\Omega} u = f \quad \mathrm{on} \ \Omega 
```
where $\Delta_{\Omega}$ is the Laplace-Beltrami operator.

## Time-dependent heat equation
An extension of the Laplace-Beltrami problem by a time dependent term gives one the time-dependent heat equation, i.e.,

Given a domain $\Omega$, $k\in\mathbb{R}$ and $f:\left( \Omega, [0,\infty) \right)\rightarrow \mathbb{R}$, find $u:\left( \Omega, [0,\infty) \right)\rightarrow \mathbb{R}$ such that
```math
\frac{\partial}{\partial t} u(x,t) - k\Delta_{\Omega} \left( u(x,t) \right) = f(x,t) \quad \forall\ x \in \Omega ,\ t\in[0,\infty)
```
where $\Delta_{\Omega}$ is the Laplace-Beltrami operator.

## Reaction-Diffusion equation (RDE) problem
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

### Gray-Scott system
```math
f(u,v; \alpha , \beta) = \alpha (1-u) - uv^2
```
```math
g(u,v; \alpha , \beta) =  -(\alpha +\beta)v + uv^2
```
### Schnakenberg system
```math
f(u,v; a,b, \gamma) = \gamma(a-u+u^2v)
```
```math
g(u,v;a,b, \gamma) = \gamma(b-u^2v)
```

# Cases considered:

1. Time-dependent heat equation on a plane (01_transient_heat_eqn_on_plane)
2. LBO problem on a unit sphere
3. LBO problem on the Dziuk surface
4. Time-dependent heat equation on the Dziuk surface
5. Gray-Scott RDE on a unit sphere
6. Schnakenberg RDE on a unit sphere

# Implementation

## Mesh creation

Using $d_\Omega$, an initial coarse triangular mesh, $\Omega_0$, with a mesh size of $h_0$, is created. This is then iteratively refined an $m$ number of times to create finer meshes, $\Omega_1, \Omega_2, ...,\Omega_{m}$, with reducing mesh sizes $h_1>h_2>...h_m$. 

```math
\mathbf{\Omega}_m := \{\Omega_i:i=0,1,...,m\},\quad H_m:=\{h_i:i=0,1,...,m \}
```
Equivalently, the notation $\Omega_h$ is sometimes used to directly link the mesh with its mesh size.

## Inputs

The main input for all problems is the domain $\rightarrowtail$ $\Omega\in\R^{n_d}$, in the form of a mesh.

### Cases 1, 2, 3 and 4
A manufactured solution, $u_{true}(x,t)$ is created analytically and from this, the function $f(x,t)$ is computed either analytically (in simple cases) or computationally (See [func_sym.jl](../include/funcs_sym.jl)).
An initial condidtion $u_{true}(x,0)$ is also set on the mesh.

### Cases 5 and 6
The exact solution to these problems were not calculated and different initial conditions (randomly perturbed spots or randomly perturbed linearly stable homogeneous steady states) were taken in diffrerent experiments.
Additionally, the functions of the reaction terms (See [struct_GrayScott.jl](../include/struct_GrayScott.jl) and [struct_Schnakenberg.jl](../include/struct_Schankenberg.jl)) are defined and their parameters are set in each experiment.

## Numerical method

The PDE systems are solved using Finite-Element-Methods (FEM).

<!--TODO: details-->

## Final simulation

### Cases 1, 2, 3 and 4
For each mesh, $\Omega_i \in \mathbf{\Omega}_m $, the PDE problem is solved and the $L_2$-error of the numerical solution (with respect to $u_{true}$) is calculated. Finally, to see whether the numerical method works correctly, the following values are printed and examined:
* refinement level
* mesh size (longest edge length)
* L2-error
* Error of Convergence (EOC).

The $L_2$-error is expected to descrease as $h$ reduces, and since the refinement stratergy is to divide each triangle into four smaller triangles, the EOC is expected to converge to 2.0.

### Case 5 : Gray-Scott RDE on a sphere
This section tests whether the custom `SurfaceCellValues` works similar to the tutorial found in Ferrite documentation ([Reactive Surfaces](https://ferrite-fem.github.io/Ferrite.jl/stable/tutorials/reactive_surface/)).

### Case 6
This case aimed to reproduce a few results from the paper _"Spatio-temporal pattern formation on spherical surfaces: numerical simulation and application to solid tumour growth"_ by M.A.J. Chaplain, M. Ganesh and I.G. Graham (2001). The Schankenberg system was simulated for three different initial conditions and parameter sets.

# Results

See the README files in the respective subfolders.