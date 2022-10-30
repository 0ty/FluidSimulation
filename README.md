[![Demo gif](resources/demo.gif)](https://youtu.be/m-MrckE2E9k)\
[**Click here for demo video**](https://youtu.be/m-MrckE2E9k)

# Fluid Simulation (Compressible Navier-Stokes Solver)

Solves the **compressible Navier-Stokes equations** in 2D in real time (on the CPU for now):

1. $\partial_t (\rho)$ $= - \partial_x (\rho u)$ $+ \partial_y (\rho v)$
2. $\partial_t (\rho u)$ $= - \partial_x ( \rho u^2 )$ $- \partial_y ( \rho u v )$ $- \partial_x (p)$ $+ \partial_x (2\mu \partial_x (u))$ $+ \partial_y( \mu ( \partial_y(u)-\partial_x(v) ) )$ $+ \rho f_u$
3. $\partial_t (\rho v)$ $= - \partial_y ( \rho v^2 )$ $- \partial_x ( \rho u v )$ $- \partial_y (p)$ $+ \partial_y (2\mu \partial_y (v))$ $+ \partial_x( \mu ( \partial_x(v)-\partial_y(u) ) )$ $+ \rho f_v$
4. $\partial_t (E)$ $= - \partial_x ( E u )$ $- \partial_y ( E v )$ $- \partial_x (p u)$ $- \partial_y (p v)$ $+ \partial_x ( 2 \mu \partial_x (u) u)$ $+ \partial_y ( 2 \mu \partial_y (v) v )$ $+ \partial_y( \mu ( \partial_y(u)-\partial_x(v) ) v)$ $+ \partial_x( \mu ( \partial_x(v)-\partial_y(u) ) u )$ $+ \rho f_u u$ $+ \rho f_v v$

The discretization and algorithm is based on the **MacCormack scheme**, translated into a finite volume framework.

## Dependencies
All external libraries are available through [vcpkg](https://github.com/Microsoft/vcpkg):
- C++20 (+OpenMP)
- [Eigen](https://eigen.tuxfamily.org/)
- [Dear ImGui](https://github.com/ocornut/imgui)
- [ImGui-SFML](https://github.com/eliasdaler/imgui-sfml)
- [SFML](https://www.sfml-dev.org/)


## Remarks
- The method is fully explicit, so it might become unstable at times.
- What is different to other real time fluid solvers like [Jos Stam's Stable Fluids](https://d2f99xq7vri1nk.cloudfront.net/legacy_app_files/pdf/ns.pdf)? His algorithm solves the incompressible case, while here the more general compressible NS equations are solved. Also, the method here is fully conservative, meaning mass, momentum and energy are conserved exactly at every step.
- For clarity I stayed very close to the mathematical formulation of the algorithm as in the literature, so at many points the code is not nearly as optimized as it can be.
- For more info I highly recommend the book *Computational Fluid Mechanics and Heat Transfer* by *Anderson, Tannehill, Pletcher*.