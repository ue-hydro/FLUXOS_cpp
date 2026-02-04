Theoretical Foundation
==================================

This section describes the mathematical and numerical foundations of FLUXOS-OVERLAND.

Governing Equations
-------------------

Shallow Water Equations (Hydrodynamics)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

FLUXOS-OVERLAND solves the 2D shallow water equations (Saint-Venant equations) in conservative form:

**Mass Conservation:**

.. math::

   \frac{\partial h}{\partial t} + \frac{\partial (hu)}{\partial x} + \frac{\partial (hv)}{\partial y} = S_h

**Momentum Conservation (x-direction):**

.. math::

   \frac{\partial (hu)}{\partial t} + \frac{\partial}{\partial x}\left(hu^2 + \frac{1}{2}gh^2\right) + \frac{\partial (huv)}{\partial y} = gh(S_{0x} - S_{fx})

**Momentum Conservation (y-direction):**

.. math::

   \frac{\partial (hv)}{\partial t} + \frac{\partial (huv)}{\partial x} + \frac{\partial}{\partial y}\left(hv^2 + \frac{1}{2}gh^2\right) = gh(S_{0y} - S_{fy})

Where:

* :math:`h` = water depth [m]
* :math:`u, v` = depth-averaged velocities in x and y directions [m/s]
* :math:`g` = gravitational acceleration [m/s²]
* :math:`S_h` = source/sink term for water [m/s]
* :math:`S_{0x}, S_{0y}` = bed slopes in x and y directions [-]
* :math:`S_{fx}, S_{fy}` = friction slopes in x and y directions [-]

Advection-Dispersion-Reaction Equation (Solute Transport)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For solute transport, the 2D ADE is solved:

.. math::

   \frac{\partial (hC)}{\partial t} + \frac{\partial (huC)}{\partial x} + \frac{\partial (hvC)}{\partial y} = \frac{\partial}{\partial x}\left(hD_x\frac{\partial C}{\partial x}\right) + \frac{\partial}{\partial y}\left(hD_y\frac{\partial C}{\partial y}\right) + S_C

Where:

* :math:`C` = solute concentration [kg/m³]
* :math:`D_x, D_y` = dispersion coefficients [m²/s]
* :math:`S_C` = source/sink term for solute [kg/m³/s]

Numerical Methods
-----------------

Finite Volume Method
^^^^^^^^^^^^^^^^^^^^

FLUXOS-OVERLAND uses a cell-centered finite volume method on a structured Cartesian grid. The domain is discretized into rectangular cells where:

* State variables (h, qx, qy, C) are stored at cell centers
* Fluxes are computed at cell faces

Roe's Approximate Riemann Solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The numerical fluxes at cell interfaces are computed using Roe's approximate Riemann solver, which provides:

* Upwind discretization for advective terms
* Proper treatment of shock waves and discontinuities
* Robust handling of transcritical flows

The Roe flux is computed as:

.. math::

   F_{i+1/2} = \frac{1}{2}(F_L + F_R) - \frac{1}{2}\sum_{k=1}^{3}|\tilde{\lambda}_k|\tilde{\alpha}_k\tilde{r}_k

Where :math:`\tilde{\lambda}_k`, :math:`\tilde{\alpha}_k`, and :math:`\tilde{r}_k` are the Roe-averaged eigenvalues, wave strengths, and eigenvectors.

MUSCL Reconstruction
^^^^^^^^^^^^^^^^^^^^

Second-order spatial accuracy is achieved using MUSCL (Monotone Upstream-centered Schemes for Conservation Laws) reconstruction with slope limiters:

* Linear reconstruction of state variables at cell faces
* MinMod or Van Leer slope limiters to prevent spurious oscillations
* TVD (Total Variation Diminishing) property preservation

Time Integration
^^^^^^^^^^^^^^^^

**Adaptive Time Stepping:**

The time step is computed based on the CFL (Courant-Friedrichs-Lewy) condition:

.. math::

   \Delta t = CFL \cdot \min_{i,j}\left(\frac{\Delta x}{\sqrt{gh_{i,j}} + |u_{i,j}|}, \frac{\Delta y}{\sqrt{gh_{i,j}} + |v_{i,j}|}\right)

Where CFL < 1.0 (typically 0.5-0.9) ensures numerical stability.

**Explicit Euler Method:**

Time advancement uses an explicit Euler scheme:

.. math::

   U^{n+1} = U^n + \Delta t \cdot R(U^n)

Where :math:`U` is the vector of conserved variables and :math:`R` is the residual.

Wetting and Drying
^^^^^^^^^^^^^^^^^^

FLUXOS-OVERLAND includes robust wetting and drying treatment:

* Cells with depth below threshold (:math:`h_{dry}`) are treated as dry
* Momentum is set to zero in dry cells
* Flux limiters prevent negative depths
* Continuous tracking of wet/dry fronts

Friction Treatment
^^^^^^^^^^^^^^^^^^

Bed friction is modeled using Manning's equation:

.. math::

   S_{fx} = \frac{n^2 u \sqrt{u^2 + v^2}}{h^{4/3}}, \quad S_{fy} = \frac{n^2 v \sqrt{u^2 + v^2}}{h^{4/3}}

Where :math:`n` is Manning's roughness coefficient [s/m^(1/3)].

Boundary Conditions
-------------------

**Dirichlet (Fixed Value):**
  Specified water level or discharge at domain boundaries

**Neumann (Zero Gradient):**
  Free outflow boundaries with zero normal gradient

**Internal Weirs:**
  Weir flow equations for internal hydraulic structures

Parallel Computing Implementation
---------------------------------

Domain Decomposition
^^^^^^^^^^^^^^^^^^^^

For MPI parallelization, the domain is decomposed using a 2D Cartesian topology:

* The global domain is divided into rectangular subdomains
* Each MPI process handles one subdomain
* Ghost cells (halo regions) are used for inter-process communication

Ghost Cell Exchange
^^^^^^^^^^^^^^^^^^^

At each time step, ghost cells are exchanged between neighboring processes:

1. Pack boundary data into send buffers
2. Non-blocking MPI communication (Isend/Irecv)
3. Wait for communication completion
4. Unpack received data into ghost regions

This ensures each process has the necessary data from neighbors for flux computation.

OpenMP Parallelization
^^^^^^^^^^^^^^^^^^^^^^

Within each MPI process, OpenMP is used for shared-memory parallelism:

* Loop parallelization with ``#pragma omp parallel for``
* Collapse(2) for efficient nested loop parallelization
* Static scheduling for predictable load balancing
* Critical sections for global reductions (CFL computation)

References
----------

Costa, D., Shook, K., Spence, C., Elliott, J., Baulch, H., Wilson, H., & Pomeroy, J. W. (2020). Predicting variable contributing areas, hydrological connectivity, and solute transport pathways for a Canadian Prairie basin. Water Resources Research, 56, e2020WR027984. https://doi.org/10.1029/2020WR027984

Costa, D., Burlando, P., and Liong, S.-Y. (2016). Coupling spatially distributed river and groundwater transport models to investigate contaminant dynamics at river corridor scales. Environmental Modelling & Software, 86:91-110, https://doi.org/10.1016/j.envsoft.2016.09.009

Roe, P.L. (1981). Approximate Riemann solvers, parameter vectors, and difference schemes. Journal of Computational Physics, 43:357-372.

Van Leer, B. (1979). Towards the ultimate conservative difference scheme. V. A second-order sequel to Godunov's method. Journal of Computational Physics, 32:101-136.
