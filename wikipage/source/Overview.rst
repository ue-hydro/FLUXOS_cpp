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
   * - **Solute Transport**
     - 2D Advection-Dispersion-Reaction (ADE) equation with source/sink terms
   * - **Numerical Scheme**
     - Finite volume method with MUSCL reconstruction for second-order accuracy
   * - **Time Stepping**
     - Adaptive time stepping based on CFL (Courant-Friedrichs-Lewy) condition
   * - **Parallel Computing**
     - Hybrid MPI+OpenMP for HPC clusters, OpenMP for workstations
   * - **Boundary Conditions**
     - Dirichlet, Neumann, and internal weir boundary conditions

Computational Features
----------------------

FLUXOS-OVERLAND is optimized for high-performance computing:

**Shared-Memory Parallelism (OpenMP)**

* Multi-threaded execution on multi-core workstations
* Automatic load balancing with static scheduling
* Loop collapse optimization for nested iterations

**Distributed-Memory Parallelism (MPI)**

* 2D Cartesian domain decomposition
* Efficient ghost cell exchange for inter-process communication
* Scalable to hundreds of processors on HPC clusters

**Hybrid MPI+OpenMP**

* Combines distributed and shared memory parallelism
* Optimal for modern HPC architectures
* Reduces memory overhead compared to pure MPI

Scalability Guidelines
----------------------

The following table provides recommended configurations based on domain size:

.. list-table::
   :widths: 25 25 25 25
   :header-rows: 1

   * - Domain Size
     - Recommended Mode
     - MPI Processes
     - OpenMP Threads
   * - < 500x500
     - OpenMP only
     - 1
     - 4-8
   * - 500x500 - 2000x2000
     - Hybrid
     - 4-16
     - 2-4
   * - 2000x2000 - 5000x5000
     - Hybrid
     - 16-64
     - 2-4
   * - > 5000x5000
     - Hybrid
     - 64-256
     - 2-4

Target Environments
-------------------

* **Workstations**: Multi-core systems with 8-64 cores
* **HPC Clusters**: SLURM, PBS, or LSF job schedulers
* **Cloud Computing**: AWS, Azure, or GCP HPC instances

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
