Overview
==================================

What is FLUXOS-OVERLAND?
------------------------

FLUXOS-OVERLAND is a high-performance, physics-based hydrodynamic and solute transport modelling tool designed for basin-scale, event-based simulations. Originally developed for Canadian Prairie hydrology, including snowmelt flooding analysis, the model has evolved into a comprehensive platform suitable for a wide range of environmental applications.

The model solves the 2D shallow water equations (Saint-Venant equations) coupled with advection-dispersion-reaction equations for solute transport, enabling integrated analysis of water flow and contaminant dynamics across complex landscapes.

Model Heritage
--------------

FLUXOS-OVERLAND is the successor to FLUXOS, which was originally developed as an integrated surface water-groundwater model. The current version focuses on overland flow processes while maintaining the rigorous numerical foundation of its predecessor.

Key Applications
----------------

* **Variable Contributing Areas**: Analysis of dynamically changing watershed contributing areas during hydrological events
* **Snowmelt Flooding**: Event-based simulation of spring snowmelt and associated flooding patterns
* **Hydrological Connectivity**: Assessment of surface water connectivity and flow pathways
* **Solute Transport**: Tracking of dissolved substances, nutrients, and contaminants across landscapes
* **Prairie Hydrology**: Specialized capabilities for pothole-dominated Prairie landscapes
* **Wetland Dynamics**: Fill-and-spill hydrology of depressional storage systems

Key Capabilities
----------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Feature
     - Description
   * - **Hydrodynamics**
     - 2D shallow water equations with wetting/drying, Roe's approximate Riemann solver
   * - **Unstructured Mesh**
     - Triangular mesh support with edge-based rotated Roe/HLL solver, MUSCL + Barth-Jespersen limiter, and hydrostatic reconstruction for well-balanced property
   * - **Solute Transport**
     - 2D Advection-Dispersion-Reaction (ADE) equation with source/sink terms (both Cartesian and edge-based formulations)
   * - **Numerical Scheme**
     - Finite volume method with MUSCL reconstruction for second-order accuracy; least-squares gradient on unstructured meshes
   * - **Time Stepping**
     - Adaptive time stepping based on CFL (Courant-Friedrichs-Lewy) condition; inradius-based CFL for triangular meshes
   * - **GPU Acceleration**
     - CUDA GPU offloading for both regular and triangular mesh solvers with 1-thread-per-cell/edge kernel design
   * - **Parallel Computing**
     - Hybrid MPI+OpenMP for HPC clusters, OpenMP for workstations; METIS graph partitioning for unstructured mesh decomposition
   * - **Boundary Conditions**
     - Dirichlet, Neumann, and internal weir boundary conditions; tag-based boundary conditions for triangular meshes (wall, outflow)
   * - **Output Formats**
     - ASCII text output for regular mesh; VTK XML Unstructured Grid (.vtu) with ParaView time series (.pvd) for triangular mesh

Mesh Types
----------

FLUXOS-OVERLAND supports two mesh types, selectable at runtime via the JSON configuration:

**Regular Cartesian Mesh (default)**

* Traditional structured grid with uniform cell size (``dxy``)
* Efficient row-column indexing with ``(irow, icol)``
* Input: ASCII-ESRI DEM files
* Output: ASCII text files

**Unstructured Triangular Mesh** (requires ``USE_TRIMESH`` build option)

* Arbitrary triangular elements for complex geometries
* Edge-based finite volume formulation (3 faces per cell, rotated Roe solver)
* Input: Gmsh ``.msh`` v2.2 or Triangle ``.node/.ele`` format
* Output: VTK XML UnstructuredGrid (``.vtu``) with ParaView time series (``.pvd``)
* Ideal for irregular domains, refined regions, and natural channel geometries

Computational Features
----------------------

FLUXOS-OVERLAND is optimized for high-performance computing:

**CUDA GPU Acceleration**

* Full GPU offloading of hydrodynamics and ADE solvers
* Regular mesh: Cartesian 2D thread grid (1 thread per cell)
* Triangular mesh: 1D thread indexing with 7 specialized kernels (wet/dry, gradient, limiter, edge flux, accumulate, update, Courant)
* Race-free flux accumulation via ``atomicAdd``
* Block-level parallel reduction for global CFL computation
* Typical speedup: 10-50x over single-core CPU for large domains

**Shared-Memory Parallelism (OpenMP)**

* Multi-threaded execution on multi-core workstations
* Automatic load balancing with static scheduling
* Loop collapse optimization for nested iterations
* Edge-parallel and cell-parallel loops for triangular mesh

**Distributed-Memory Parallelism (MPI)**

* Regular mesh: 2D Cartesian domain decomposition with ghost cell exchange
* Triangular mesh: graph-based partitioning (METIS optional, naive block fallback)
* Halo exchange for unstructured mesh neighbor cells
* Scalable to hundreds of processors on HPC clusters

**Hybrid MPI+OpenMP+CUDA**

* Combines distributed, shared memory, and GPU parallelism
* Optimal for modern HPC architectures with GPU nodes
* Reduces memory overhead compared to pure MPI

Scalability Guidelines
----------------------

The following table provides recommended configurations based on domain size:

.. list-table::
   :widths: 20 20 15 15 30
   :header-rows: 1

   * - Domain Size
     - Recommended Mode
     - MPI Processes
     - OpenMP Threads
     - Notes
   * - < 500x500 / < 100K cells
     - OpenMP or CUDA
     - 1
     - 4-8
     - GPU provides best speedup even at small scale
   * - 500x500 - 2000x2000
     - Hybrid or CUDA
     - 4-16
     - 2-4
     - Good balance of communication/computation
   * - 2000x2000 - 5000x5000
     - Hybrid+CUDA
     - 16-64
     - 2-4
     - Multi-GPU recommended
   * - > 5000x5000 / > 1M cells
     - Hybrid+CUDA
     - 64-256
     - 2-4
     - Large-scale HPC with GPU nodes

Target Environments
-------------------

* **Workstations**: Multi-core systems with 8-64 cores; NVIDIA GPU (Compute Capability 6.0+)
* **HPC Clusters**: SLURM, PBS, or LSF job schedulers; multi-GPU nodes
* **Cloud Computing**: AWS (P3/P4 instances), Azure (NC series), or GCP (A2 instances)

Typical Use Cases
-----------------

1. **Watershed Event Analysis**: Simulate individual storm or snowmelt events to understand flow patterns and contributing areas
2. **Nutrient Transport Studies**: Track agricultural nutrient movement across landscapes
3. **Flood Risk Assessment**: Evaluate flood extent and timing for infrastructure planning
4. **Climate Impact Studies**: Assess hydrological response to changing precipitation patterns
5. **Wetland Restoration Planning**: Model fill-and-spill dynamics for wetland design

References
----------

Costa, D., Shook, K., Spence, C., Elliott, J., Baulch, H., Wilson, H., & Pomeroy, J. W. (2020). Predicting variable contributing areas, hydrological connectivity, and solute transport pathways for a Canadian Prairie basin. Water Resources Research, 56, e2020WR027984. https://doi.org/10.1029/2020WR027984

Costa, D., Burlando, P., and Liong, S.-Y. (2016). Coupling spatially distributed river and groundwater transport models to investigate contaminant dynamics at river corridor scales. Environmental Modelling & Software, 86:91-110, https://doi.org/10.1016/j.envsoft.2016.09.009
