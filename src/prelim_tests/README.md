# Aim
The purpose of these preliminary tests is to solve specific partial differential equation (PDE) systems on specific domains, and ensure that the numerical method(s) work correctly.

## Laplace-Beltrami problem
The Laplace-Beltrami problem is defined as follows: Given a function, $f:\Omega \rightarrow \mathbb{R}$, on a domain, $\Omega$, find $u:\Omega \rightarrow \mathbb{R}$ such that
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

# Implementation

## Inputs and computation of $f$

* Domain $\rightarrowtail$ $\Omega\in\R^{n_d}$
    - (for embedded 2-dimensional surfaces) Signed distance function $\rightarrowtail$ $d_{\Omega}:\R^{3}\rightarrow\R$ such that $x\in\Omega\iff d_{\Omega}(x)=0$.
* Manufactured solution $\rightarrowtail$ $u_{true}(x,t)$

where $x=[x_1,...,x_{n_d}]$.

The function $f$ is computed either analytically (in simple cases) or computationally (See [func_sym.jl](../include/funcs_sym.jl)).

## Mesh creation

Using $d_\Omega$, an initial coarse triangular mesh, $\Omega_0$, with a mesh size of $h_0$, is created. This is then iteratively refined an $m$ number of times to create finer meshes, $\Omega_1, \Omega_2, ...,\Omega_{m}$, with reducing mesh sizes $h_1>h_2>...h_m$. 

```math
\mathbf{\Omega}_m := \{\Omega_i:i=0,1,...,m\},\quad H_m:=\{h_i:i=0,1,...,m \}
```

## Numerical method

The PDE systems are solved using Finite-Element-Methods (FEM). See respective subfolders for details.

## Final simulation
For each mesh, $\Omega_i \in \mathbf{\Omega}_m $, the PDE problem is solved and the $L_2$-error of the numerical solution (with respect to $u_{true}$) is calculated. Finally, to see whether the numerical method works correctly, the following values are printed and examined:
* refinement level
* mesh size (longest edge length)
* L2-error
* Error of Convergence (EOC).

The $L_2$-error is expected to descrease as $h$ reduces, and since the refinement stratergy is to divide each triangle into four smaller triangles, the EOC is expected to converge to 2.0.

# Results

See the README files in the respective subfolders.