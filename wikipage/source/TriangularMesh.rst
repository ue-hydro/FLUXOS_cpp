Triangular Mesh Guide
==================================

This section provides a comprehensive guide to using unstructured triangular meshes in FLUXOS-OVERLAND.

Overview
--------

FLUXOS-OVERLAND supports unstructured triangular meshes as an alternative to the default regular Cartesian grid. The triangular mesh solver uses an edge-based finite volume formulation with a rotated Roe/HLL flux solver, providing second-order accuracy on arbitrary mesh geometries.

**When to use triangular meshes:**

* Irregular domain boundaries (coastlines, river banks, infrastructure)
* Local mesh refinement needed (e.g., finer resolution near channels or structures)
* Complex channel geometries that are poorly represented by Cartesian grids
* Domains where uniform resolution is wasteful (large flat areas with small features)

**When to use regular meshes:**

* Rectangular domains with uniform features
* Maximum computational efficiency on simple geometries
* Compatibility with existing Cartesian-based workflows

Build Requirements
------------------

Triangular mesh support requires the ``USE_TRIMESH`` CMake option:

.. code-block:: bash

   cmake -DMODE_release=ON -DUSE_TRIMESH=ON ..
   make -j$(nproc)

For GPU acceleration of the triangular mesh solver:

.. code-block:: bash

   cmake -DMODE_release=ON -DUSE_TRIMESH=ON -DUSE_CUDA=ON ..
   make -j$(nproc)

Mesh Formats
------------

FLUXOS-OVERLAND supports two standard triangular mesh formats:

Gmsh .msh v2.2
^^^^^^^^^^^^^^^

`Gmsh <https://gmsh.info/>`_ is a popular open-source mesh generator. FLUXOS reads the ``.msh`` v2.2 ASCII format.

**Key sections parsed:**

* ``$Nodes``: Vertex coordinates (x, y, z)
* ``$Elements``: Triangles (type 2) and boundary lines (type 1)
* Physical group tags map to boundary condition types

**Creating a mesh in Gmsh:**

1. Define the domain geometry in Gmsh (points, lines, surfaces)
2. Assign physical groups to boundary lines (e.g., group 1 = wall, group 2 = outflow)
3. Generate the mesh: ``Mesh > 2D``
4. Save as ``.msh`` format version 2.2

.. code-block:: text

   // Example Gmsh .geo file
   Point(1) = {0, 0, 0, 1.0};
   Point(2) = {100, 0, 0, 1.0};
   Point(3) = {100, 50, 0, 1.0};
   Point(4) = {0, 50, 0, 1.0};

   Line(1) = {1, 2};    // bottom wall
   Line(2) = {2, 3};    // outflow
   Line(3) = {3, 4};    // top wall
   Line(4) = {4, 1};    // inflow wall

   Line Loop(1) = {1, 2, 3, 4};
   Plane Surface(1) = {1};

   Physical Line("wall") = {1, 3, 4};   // tag 1
   Physical Line("outflow") = {2};       // tag 2
   Physical Surface("domain") = {1};

Triangle .node/.ele
^^^^^^^^^^^^^^^^^^^

The `Triangle <https://www.cs.cmu.edu/~quake/triangle.html>`_ mesh generator produces three files:

* ``.node``: Vertex coordinates
* ``.ele``: Triangle connectivity
* ``.edge``: Edge information (optional)

Specify the base filename (without extension) in the configuration.

JSON Configuration
------------------

To use a triangular mesh, add these settings to your master JSON file:

.. code-block:: json

   {
       "MESH_TYPE": "triangular",
       "MESH_FILE": "path/to/domain.msh",
       "MESH_FORMAT": "gmsh",
       "BOUNDARY_CONDITIONS": {
           "1": {"type": "wall"},
           "2": {"type": "outflow"}
       }
   }

**Configuration options:**

.. list-table::
   :widths: 25 20 55
   :header-rows: 1

   * - Key
     - Values
     - Description
   * - ``MESH_TYPE``
     - ``"regular"`` or ``"triangular"``
     - Selects the mesh type (default: ``"regular"``)
   * - ``MESH_FILE``
     - file path
     - Path to the mesh file (relative to the JSON config directory)
   * - ``MESH_FORMAT``
     - ``"gmsh"`` or ``"triangle"``
     - Mesh file format
   * - ``BOUNDARY_CONDITIONS``
     - JSON object
     - Maps physical group tags to boundary types

**Boundary condition types:**

* ``"wall"`` (default): Reflective boundary - normal velocity is reflected, tangential velocity preserved
* ``"outflow"``: Transmissive boundary - zero-gradient free outflow

.. note::

   The DEM file (``DEM_FILE``) is still required even with triangular meshes. It provides topographic data that is interpolated to cell centroids during initialization.

Numerical Methods
-----------------

The triangular mesh solver employs the following algorithms:

**6-Phase Hydrodynamics Driver:**

Each timestep performs:

1. **Wet/dry classification**: Identify dry cells and update status
2. **Gradient computation**: Least-squares gradient using neighbor cell values with precomputed 2x2 inverse matrix
3. **Slope limiting**: Barth-Jespersen limiter applied per cell to ensure monotonicity
4. **Edge flux computation**: MUSCL-reconstructed left/right states at each edge, rotated Roe/HLL solver in edge-normal frame
5. **Flux accumulation**: Sum edge fluxes into cell residuals
6. **State update**: Advance state variables, apply friction, enforce wet/dry constraints

**CFL Condition:**

The time step uses the cell inradius (inscribed circle radius) for stability:

.. math::

   \Delta t = CFL \cdot \min_i \frac{r_i}{|\mathbf{u}_i| + \sqrt{g h_i}}, \quad r_i = \frac{2 A_i}{P_i}

This naturally accounts for varying cell sizes across the mesh.

**Edge-Based ADE Solver:**

Solute transport uses:

* Upwind advection based on hydrodynamic mass flux direction at each edge
* Fick's law dispersion with effective coefficient
* Mass conservation limiting (outgoing flux bounded by available mass)
* Concentration bounding by neighbor min/max

Output Format
-------------

Triangular mesh results are output as VTK XML Unstructured Grid files (``.vtu``):

* One ``.vtu`` file per output timestep
* A ``.pvd`` time series collection file for ParaView animation
* Cell-centered data arrays: h, z, zb, ux, uy, velocity_magnitude, qx, qy, us, conc_SW, twetimetracer, cell_area

See the :doc:`file` page for detailed output format documentation.

GPU Acceleration
----------------

When built with both ``USE_TRIMESH`` and ``USE_CUDA``, the triangular mesh solver uses 7 specialized CUDA kernels:

.. list-table::
   :widths: 30 20 50
   :header-rows: 1

   * - Kernel
     - Parallelism
     - Description
   * - ``kernel_tri_wetdry``
     - 1 thread/cell
     - Wet/dry classification and status update
   * - ``kernel_tri_gradient``
     - 1 thread/cell
     - Least-squares gradient computation
   * - ``kernel_tri_limiter``
     - 1 thread/cell
     - Barth-Jespersen slope limiter
   * - ``kernel_tri_edge_flux``
     - 1 thread/edge
     - Rotated Roe/HLL flux (main compute kernel)
   * - ``kernel_tri_accumulate``
     - 1 thread/edge
     - Flux accumulation with ``atomicAdd``
   * - ``kernel_tri_update``
     - 1 thread/cell
     - State update + wet/dry clamping
   * - ``kernel_tri_courant``
     - 1 thread/cell
     - CFL time step reduction

MPI Parallelization
-------------------

For distributed computing with triangular meshes:

* **Partitioning**: METIS graph partitioning (if available) or naive block partitioning
* **Halo exchange**: Cells adjacent to partition boundaries are exchanged between MPI ranks
* **Communication**: ``MPI_Isend``/``MPI_Irecv`` per neighbor rank

Build with MPI:

.. code-block:: bash

   cmake -DMODE_release=ON -DUSE_TRIMESH=ON -DUSE_MPI=ON ..
   make -j$(nproc)
   mpirun -np 4 ./bin/fluxos_mpi ./input/modset_triangular.json

Source Code Organization
------------------------

All triangular mesh code is contained in ``src/fluxos/trimesh/``:

.. list-table::
   :widths: 35 65
   :header-rows: 1

   * - File
     - Purpose
   * - ``TriMesh.h/.cpp``
     - Mesh topology, geometry, edge building, LSQ matrices
   * - ``TriSolution.h/.cpp``
     - Flat-array solution storage and allocation
   * - ``tri_mesh_io.h/.cpp``
     - Gmsh and Triangle format readers
   * - ``tri_initiate.h/.cpp``
     - Initialization, DEM interpolation, boundary conditions
   * - ``tri_gradient.h/.cpp``
     - LSQ gradient + Barth-Jespersen limiter
   * - ``tri_solver_wet.h/.cpp``
     - Rotated Roe/HLL flux with MUSCL reconstruction
   * - ``tri_solver_dry.h/.cpp``
     - Ritter dry-front solution
   * - ``tri_hydrodynamics.h/.cpp``
     - 6-phase driver + CFL condition
   * - ``tri_ade_solver.h/.cpp``
     - Edge-based advection-dispersion
   * - ``tri_forcing.h/.cpp``
     - Meteo and inflow forcing
   * - ``tri_write_results.h/.cpp``
     - VTK .vtu output + .pvd time series
   * - ``tri_mpi_domain.h/.cpp``
     - MPI decomposition + halo exchange
   * - ``tri_cuda_memory.h/.cu``
     - GPU memory management
   * - ``tri_cuda_kernels.h/.cu``
     - CUDA kernels (7 kernels)
