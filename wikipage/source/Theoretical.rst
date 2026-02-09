Theoretical Foundation
==================================

This section describes the mathematical and numerical foundations of FLUXOS.

Governing Equations
-------------------

Shallow Water Equations (Hydrodynamics)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

FLUXOS solves the 2D shallow water equations (Saint-Venant equations) in conservative form:

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

FLUXOS uses a cell-centered finite volume method. Two mesh types are supported:

**Regular Cartesian Mesh:**

* Structured grid with rectangular cells of uniform size ``dxy``
* State variables (h, qx, qy, C) stored at cell centers
* Fluxes computed at cell faces in x and y directions (directional sweep)

**Unstructured Triangular Mesh:**

* Arbitrary triangular cells with edge-based flux computation
* State variables stored at cell centroids
* Fluxes computed at each edge, rotated into the edge-normal frame
* Cell update: :math:`\Delta U_i = -\frac{\Delta t}{A_i} \sum_{e \in \partial \Omega_i} F_e \cdot L_e`

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

**Rotated Roe Solver (Triangular Mesh):**

For the unstructured triangular mesh, the Roe solver operates in an edge-normal reference frame:

1. MUSCL-reconstruct left/right states to the edge midpoint
2. Apply hydrostatic reconstruction: :math:`z_{b,face} = \max(z_{b,L}, z_{b,R})`
3. Rotate velocities into edge-normal frame: :math:`q_n = q_x n_x + q_y n_y`, :math:`q_t = -q_x n_y + q_y n_x`
4. Compute HLL/Roe flux in the 1D normal frame
5. Add turbulent stress contribution
6. Rotate flux back to global frame: :math:`F_x = F_n n_x - F_t n_y`, :math:`F_y = F_n n_y + F_t n_x`
7. Multiply by edge length

MUSCL Reconstruction
^^^^^^^^^^^^^^^^^^^^

Second-order spatial accuracy is achieved using MUSCL (Monotone Upstream-centered Schemes for Conservation Laws) reconstruction with slope limiters:

**Regular Mesh:**

* Linear reconstruction of state variables at cell faces
* MinMod or Van Leer slope limiters to prevent spurious oscillations
* TVD (Total Variation Diminishing) property preservation

**Triangular Mesh:**

* Least-squares gradient computation using neighbor cell values
* Precomputed 2x2 LSQ inverse matrix per cell for efficiency
* Barth-Jespersen slope limiter for monotonicity preservation on unstructured meshes:

.. math::

   \phi_i = \min_{e \in \partial \Omega_i} \begin{cases} \min\left(1, \frac{U_{\max} - U_i}{\tilde{U}_e - U_i}\right) & \text{if } \tilde{U}_e > U_i \\ \min\left(1, \frac{U_{\min} - U_i}{\tilde{U}_e - U_i}\right) & \text{if } \tilde{U}_e < U_i \\ 1 & \text{otherwise} \end{cases}

where :math:`\tilde{U}_e` is the reconstructed value at edge midpoint, and :math:`U_{\max}`, :math:`U_{\min}` are the maximum and minimum values among the cell and its neighbors.

Hydrostatic Reconstruction (C-Property)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The triangular mesh solver uses hydrostatic reconstruction to preserve the well-balanced property (C-property / lake-at-rest):

.. math::

   z_{b,face} = \max(z_{b,L}, z_{b,R}), \quad h_L^* = \max(0, z_L - z_{b,face}), \quad h_R^* = \max(0, z_R - z_{b,face})

This ensures that a lake at rest with varying bathymetry produces zero flux, preventing spurious currents over uneven beds.

Time Integration
^^^^^^^^^^^^^^^^

**Adaptive Time Stepping:**

The time step is computed based on the CFL (Courant-Friedrichs-Lewy) condition:

*Regular mesh:*

.. math::

   \Delta t = CFL \cdot \min_{i,j}\left(\frac{\Delta x}{\sqrt{gh_{i,j}} + |u_{i,j}|}, \frac{\Delta y}{\sqrt{gh_{i,j}} + |v_{i,j}|}\right)

*Triangular mesh (inradius-based):*

.. math::

   \Delta t = CFL \cdot \min_i \frac{r_i}{\sqrt{u_i^2 + v_i^2} + \sqrt{g h_i}}

Where :math:`r_i = 2 A_i / P_i` is the inradius (inscribed circle radius) of cell :math:`i`, precomputed from the mesh geometry. This ensures stability on cells of varying size and shape.

Where CFL < 1.0 (typically 0.5-0.9) ensures numerical stability.

**Explicit Euler Method:**

Time advancement uses an explicit Euler scheme:

.. math::

   U^{n+1} = U^n + \Delta t \cdot R(U^n)

Where :math:`U` is the vector of conserved variables and :math:`R` is the residual.

Wetting and Drying
^^^^^^^^^^^^^^^^^^

FLUXOS includes robust wetting and drying treatment:

* Cells with depth below threshold (:math:`h_{dry}`) are treated as dry
* Momentum is set to zero in dry cells
* Flux limiters prevent negative depths
* Continuous tracking of wet/dry fronts

**Velocity Desingularization (Triangular Mesh):**

To handle the singularity at :math:`h \to 0`, the triangular mesh solver uses:

.. math::

   u = \frac{2 h q}{h^2 + \max(h^2, \varepsilon^2)}

Where :math:`\varepsilon = h_{dry}`. This smoothly transitions velocities to zero as depth vanishes, preventing division-by-zero and large spurious velocities.

**Ritter Dry-Front Solution:**

At wet-dry interfaces on the triangular mesh, the analytical Ritter dam-break solution is used to compute the flux across the edge, providing a physics-based treatment of the advancing wet front.

Friction Treatment
^^^^^^^^^^^^^^^^^^

Bed friction is modeled using Manning's equation:

.. math::

   S_{fx} = \frac{n^2 u \sqrt{u^2 + v^2}}{h^{4/3}}, \quad S_{fy} = \frac{n^2 v \sqrt{u^2 + v^2}}{h^{4/3}}

Where :math:`n` is Manning's roughness coefficient [s/m^(1/3)].

Edge-Based ADE Solver (Triangular Mesh)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The advection-dispersion equation on the triangular mesh uses an edge-based formulation:

* **Advection:** Upwind scheme based on the hydrodynamic mass flux direction at each edge
* **Dispersion:** Fick's law with effective dispersion coefficient: :math:`F_{disp} = -D_{eff} \cdot \frac{C_R - C_L}{d_{LR}} \cdot L_e`
* **Mass conservation:** Outgoing flux limited to available mass in the donor cell
* **Bounding:** Concentrations bounded by neighbor min/max to prevent spurious oscillations

Boundary Conditions
-------------------

**Regular Mesh:**

* **Dirichlet (Fixed Value):** Specified water level or discharge at domain boundaries
* **Neumann (Zero Gradient):** Free outflow boundaries with zero normal gradient
* **Internal Weirs:** Weir flow equations for internal hydraulic structures

**Triangular Mesh:**

Boundary conditions are assigned via mesh physical group tags in the Gmsh/Triangle mesh file:

* **Wall (default, tag 1):** Reflective boundary - normal velocity reflected, tangential preserved
* **Outflow (tag 2):** Transmissive boundary - zero-gradient (free outflow)
* Custom boundary tags can be defined in the JSON configuration

Parallel Computing Implementation
---------------------------------

Domain Decomposition
^^^^^^^^^^^^^^^^^^^^

**Regular Mesh (MPI):**

* 2D Cartesian domain decomposition using MPI_Cart_create
* The global domain is divided into rectangular subdomains
* Each MPI process handles one subdomain
* Ghost cells (halo regions) used for inter-process communication

**Triangular Mesh (MPI):**

* Graph-based partitioning using METIS (if available) or naive block partitioning
* Halo cells: cells across partition-boundary edges exchanged as packed buffers
* Per-neighbor communication via ``MPI_Isend``/``MPI_Irecv``

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
* Collapse(2) for efficient nested loop parallelization (regular mesh)
* Edge-parallel and cell-parallel loops (triangular mesh)
* Static scheduling for predictable load balancing
* Reduction operations for global CFL computation

CUDA GPU Acceleration
^^^^^^^^^^^^^^^^^^^^^

**Regular Mesh Kernels:**

* 2D thread grid matching the domain (1 thread per cell)
* Separate kernels for Courant condition, hydrodynamics, and ADE
* Block-level reduction for global minimum (CFL time step)

**Triangular Mesh Kernels (7 specialized kernels):**

1. ``kernel_tri_wetdry`` - wet/dry classification (1 thread/cell)
2. ``kernel_tri_gradient`` - least-squares gradient per cell (1 thread/cell)
3. ``kernel_tri_limiter`` - Barth-Jespersen limiter per cell (1 thread/cell)
4. ``kernel_tri_edge_flux`` - rotated Roe flux at each edge (1 thread/edge) - main compute kernel
5. ``kernel_tri_accumulate`` - flux summation into cells using ``atomicAdd`` (1 thread/edge)
6. ``kernel_tri_update`` - state update + wet/dry clamping (1 thread/cell)
7. ``kernel_tri_courant`` - CFL reduction with block-level parallel minimum (1 thread/cell)

No race conditions: edge flux kernel writes only to its own edge, and ``atomicAdd`` ensures safe concurrent accumulation.

References
----------

Costa, D., Shook, K., Spence, C., Elliott, J., Baulch, H., Wilson, H., & Pomeroy, J. W. (2020). Predicting variable contributing areas, hydrological connectivity, and solute transport pathways for a Canadian Prairie basin. Water Resources Research, 56, e2020WR027984. https://doi.org/10.1029/2020WR027984

Costa, D., Burlando, P., and Liong, S.-Y. (2016). Coupling spatially distributed river and groundwater transport models to investigate contaminant dynamics at river corridor scales. Environmental Modelling & Software, 86:91-110, https://doi.org/10.1016/j.envsoft.2016.09.009

Roe, P.L. (1981). Approximate Riemann solvers, parameter vectors, and difference schemes. Journal of Computational Physics, 43:357-372.

Van Leer, B. (1979). Towards the ultimate conservative difference scheme. V. A second-order sequel to Godunov's method. Journal of Computational Physics, 32:101-136.
