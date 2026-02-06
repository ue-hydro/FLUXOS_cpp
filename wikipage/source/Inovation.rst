Inovation
==================================

    FLUXOS addresses a major modelling gap in Canada related to the characterisation of hydraulic properties at basin scales, which can be important to study complex phenomena such as variable contributing areas in areas found in the Canadian Prairies. Hydrological models tend to rely on simplified routing shemes to reduce computational cost, which do not allow a proper physically-based characterization of spatiotemporal patterns.

    FLUXOS solves the dynamic wave for surface water hydrodynamics, also known as the 2-D shallow-water partial differential equations (PDEs) derived from the depth integration of the Navier-Stokes PDEs (i.e., Saint Venant Equations) - that accounts for inertial (local and convective), pressure, gravity, and friction forces (momentum balance). Accounting for inertial forces is an essential, innovative aspect of this model. This is because, unlike with the commonly used kinematic or diffusion waves it captures back-flooding that is an important hydraulic phenomenon in low-relief terrains such as found in the Canadian Prairies. This unique ability has been explored to e.g., to examine hydrological and hydro-chemical connectivity, variable contributing areas, and predominant transport pathways.

Unstructured Triangular Mesh Support
--------------------------------------

FLUXOS now supports unstructured triangular meshes alongside the traditional regular Cartesian grid. This is a significant advancement that enables:

* **Complex geometries**: Triangular meshes can conform to irregular domain boundaries, river channels, and infrastructure without staircase approximations
* **Local mesh refinement**: Finer resolution where needed (e.g., near channels, weirs, or areas of interest) and coarser resolution elsewhere, reducing computational cost while maintaining accuracy
* **Edge-based finite volume formulation**: The solver uses a rotated Roe/HLL flux computation at each edge, providing accurate flux resolution across arbitrary cell orientations
* **Numerical dispersion mitigation**: Five techniques work together to minimize numerical diffusion on unstructured meshes:

  1. **MUSCL reconstruction** with least-squares gradients (second-order accuracy)
  2. **Barth-Jespersen limiter** (monotonicity-preserving, less diffusive than minmod)
  3. **Hydrostatic reconstruction** (preserves lake-at-rest / C-property)
  4. **Velocity desingularization** for robust wetting/drying
  5. **Ritter dry-front solution** at wet/dry interfaces

* **Standard mesh formats**: Supports Gmsh ``.msh`` v2.2 and Triangle ``.node/.ele`` formats
* **ParaView-compatible output**: VTK XML Unstructured Grid (``.vtu``) with time series (``.pvd``)

The triangular mesh is selected at runtime via the JSON configuration (``"MESH_TYPE": "triangular"``), with no changes needed to the existing regular mesh workflow.

CUDA GPU Acceleration
---------------------

FLUXOS now includes CUDA GPU acceleration, enabling significant speedups for large-scale simulations:

* **Regular mesh GPU solver**: Full GPU offloading of hydrodynamics and ADE solvers with 2D thread grid matching the domain
* **Triangular mesh GPU solver**: 7 specialized CUDA kernels optimized for unstructured mesh computation:

  - Wet/dry classification, gradient computation, Barth-Jespersen limiting
  - Edge-based rotated Roe flux (main compute kernel: 1 thread per edge)
  - Race-free flux accumulation via ``atomicAdd``
  - State update and CFL reduction

* **Typical performance**: 10-50x speedup over single-core CPU for domains with millions of cells
* **Minimal code changes**: GPU acceleration is enabled via a single CMake flag (``-DUSE_CUDA=ON``)

The GPU acceleration is complementary to the existing MPI+OpenMP parallelization, enabling multi-GPU execution on HPC clusters.

GeoTIFF DEM Support and Python Preprocessing
----------------------------------------------

FLUXOS now includes a Python preprocessing tool (``fluxos_setup.py``) that streamlines the workflow from raw GeoTIFF DEM data to simulation-ready inputs:

* **GeoTIFF DEM import**: Read GeoTIFF files directly using ``rasterio``, with automatic CRS validation and metadata inspection
* **DEM downscaling**: Resample DEMs to lower resolutions using bilinear interpolation for use with the regular Cartesian mesh
* **Slope-based adaptive mesh generation**: Automatically create triangular meshes from DEM topography using Gmsh, where steeper slopes produce finer mesh elements and flat areas produce coarser elements
* **DEM-in-mesh embedding**: DEM elevations are embedded as vertex z-coordinates in the Gmsh ``.msh`` file, eliminating the need for separate DEM interpolation during initialization
* **JSON config generation**: Automatically produce ``modset.json`` configuration files matching the C++ solver's expected format

This Python-only approach avoids adding GDAL as a C++ dependency. The C++ solver continues reading the same file formats (``.asc`` and ``.msh``), with backward-compatible support for mesh files with or without embedded vertex elevations.

Google Earth KML Animation
----------------------------

FLUXOS includes a KML exporter (``fluxos_viewer.py``) that converts simulation results to geo-referenced KMZ files for animated visualization in Google Earth:

* **Geo-referenced output**: UTM-projected results are automatically transformed to WGS84 lat/lon coordinates
* **Time-aware animation**: Each timestep is exported with KML ``TimeSpan`` tags, enabling playback via Google Earth's built-in time slider
* **Both mesh types**: Reads regular mesh ASCII output (``.txt``) and triangular mesh VTK output (``.vtu``/``.pvd``) directly
* **Multiple variables**: Water depth, velocity, or solute concentration can be exported with a blue colour ramp
* **Lightweight dependencies**: Only requires ``numpy`` and ``pyproj`` — no VTK or GUI libraries needed

See the :doc:`SupportingScripts` page for usage instructions and available options.
