Supporting Scripts
==================================

FLUXOS includes two collections of supporting scripts for preprocessing, post-processing, and visualization of simulation results.

.. contents:: In this section
   :local:
   :depth: 2

Python Preprocessing Tools (``fluxos_preprocessing/``)
-------------------------------------------------------

The ``fluxos_preprocessing/`` directory contains modern Python CLI tools for DEM preprocessing and results visualization. These tools are installed via:

.. code-block:: bash

   cd fluxos_preprocessing
   pip install -r requirements.txt

Google Earth KML Exporter (``fluxos_viewer.py``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Export simulation results as KMZ files for animated visualization in Google Earth. The exporter converts UTM-projected FLUXOS output to geo-referenced KML polygons with time spans, enabling playback via Google Earth's built-in time slider. Supports both regular and triangular mesh outputs.

.. code-block:: bash

   # Regular mesh results (auto-detect UTM zone)
   python fluxos_viewer.py --results-dir ./Results --dem ./terrain.asc

   # Triangular mesh results
   python fluxos_viewer.py --results-dir ./fluxos_out --dem ./terrain.asc --mesh-type triangular

   # Custom UTM zone, variable, and colour range
   python fluxos_viewer.py --results-dir ./Results --dem ./terrain.asc \
       --utm-zone 10 --variable velocity --clim 0 1.5

   # Reduce file size for regular mesh (export every 3rd cell)
   python fluxos_viewer.py --results-dir ./Results --dem ./terrain.asc --cell-skip 3

**Options:**

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Option
     - Default
     - Description
   * - ``--results-dir``
     - (required)
     - Directory containing output files (``.txt`` or ``.vtu``)
   * - ``--dem``
     - (required)
     - Path to ESRI ASCII DEM file (``.asc``)
   * - ``--mesh-type``
     - regular
     - ``regular`` (Cartesian ``.txt``) or ``triangular`` (VTK ``.vtu``)
   * - ``--variable``
     - h
     - Variable to color: ``h`` (water depth), ``velocity``, or ``conc_SW``
   * - ``--utm-zone``
     - auto
     - UTM zone number (1–60); auto-detected from DEM coordinates if omitted
   * - ``--clim``
     - auto
     - Colour range for the variable (e.g., ``--clim 0 2.0``)
   * - ``--opacity``
     - 180
     - Water polygon opacity (0–255)
   * - ``--cell-skip``
     - 1
     - (Regular mesh only) Export every Nth cell to reduce file size
   * - ``-o / --output``
     - fluxos_<type>.kmz
     - Output KMZ file path

**Viewing in Google Earth:**

1. Open the generated ``.kmz`` file in Google Earth (Pro or Web)
2. Use the **time slider** (top-left) to animate through simulation timesteps
3. Water depth is shown as blue polygons clamped to the terrain surface
4. Each timestep is a separate time-aware folder — drag the slider to scrub through them

GeoTIFF DEM Tool (``fluxos_setup.py dem``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Inspect, downscale, and convert GeoTIFF DEM files to ESRI ASCII format for use with the regular mesh solver.

.. code-block:: bash

   # Print GeoTIFF metadata
   python fluxos_setup.py dem --input terrain.tif --info

   # Convert GeoTIFF to ASC at native resolution
   python fluxos_setup.py dem --input terrain.tif --output-asc terrain.asc

   # Downscale DEM to 10m resolution and export to ASC
   python fluxos_setup.py dem --input terrain.tif --output-asc terrain_10m.asc --resolution 10

Adaptive Mesh Generator (``fluxos_setup.py mesh``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Generate slope-based adaptive triangular meshes from GeoTIFF DEM data using Gmsh. Steeper terrain produces finer mesh elements, flat areas produce coarser elements.

.. code-block:: bash

   python fluxos_setup.py mesh --input terrain.tif --output domain.msh \
       --min-size 5.0 --max-size 50.0 --slope-factor 1.0

See :doc:`Digital` for detailed documentation of the DEM and mesh generation tools.

Configuration Generator (``fluxos_setup.py config``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Generate ``modset.json`` configuration files for the FLUXOS solver.

.. code-block:: bash

   # Regular mesh config
   python fluxos_setup.py config --dem-file terrain.asc --mesh-type regular \
       --meteo-file forcing.fluxos --roughness 0.005 --output modset.json

   # Triangular mesh config
   python fluxos_setup.py config --dem-file domain_dem.asc --mesh-type triangular \
       --mesh-file domain.msh --meteo-file forcing.fluxos --output modset.json

Legacy Post-Processing Scripts (``fluxos_supporting_scripts/``)
----------------------------------------------------------------

The ``fluxos_supporting_scripts/`` directory contains the original collection of Python and MATLAB scripts for examining simulation results.

Python Scripts
^^^^^^^^^^^^^^^

.. list-table::
   :widths: 35 65
   :header-rows: 1

   * - Script
     - Description
   * - ``main.py``
     - Master post-processing script that orchestrates the tools below. Supports batch processing of multiple simulations with parallel execution (``joblib``).
   * - ``vtk_generator.py``
     - Converts regular mesh ``.txt`` output files to VTK XML format for visualization in ParaView. Uses ``vtktools.py`` for XML writing.
   * - ``vtktools.py``
     - Low-level VTK XML writer using Python DOM (``minidom``). Generates Unstructured Grid files from point-cloud data.
   * - ``googleearth_klmgen.py``
     - Generates KML files for overlaying inundation maps on Google Earth. Creates PNG raster images of simulation results and geo-references them.
   * - ``cross_section_extract.py``
     - Extracts cross-section data (e.g., flow, water depth) along user-defined transects through the simulation domain.
   * - ``data_management.py``
     - Data extraction utilities: column extraction from CSV (``xyz_extract_z_column``), conversion to matrix format (``xyz_to_matrix``), and observation data parsing (``obsextract``).
   * - ``graphing_functions.py``
     - Plotly-based graphing functions for plotting cross-section values, comparing simulated vs. observed data, and time series analysis.
   * - ``dem_sim_3dmap.py``
     - 3D surface map visualization using Plotly (``plotly.graph_objs.Surface``). Renders DEM and simulation results as interactive 3D surfaces.
   * - ``qmelt_tair-index.py``
     - Temperature-index snowmelt model for generating snowmelt forcing data.

MATLAB Scripts
^^^^^^^^^^^^^^^

.. list-table::
   :widths: 35 65
   :header-rows: 1

   * - Script
     - Description
   * - ``GeoMaping_results.m``
     - Plots simulation results or DEM data on a georeferenced map.
   * - ``Analyze_Results.m``
     - Comprehensive results analysis script for specific study configurations.
   * - ``Analyze_Results_Generic.m``
     - Generic version of the results analysis script for general use.
   * - ``post_processing.m``
     - Post-processing framework with support for multiple analysis types.
   * - ``Sensitivity_tests_tipping_points_paper.m``
     - Sensitivity analysis tools for parameter exploration.
   * - ``get_resultdir_list.m``
     - Utility to list and organize result directories from batch simulations.
   * - ``deg2utm.m`` / ``utm2deg.m``
     - Coordinate conversion between geographic (lat/lon) and UTM coordinates.
   * - ``othercolor.m``
     - Extended colormap definitions for MATLAB plotting.
   * - ``parfor_progressbar.m``
     - Progress bar utility for parallel MATLAB loops.

External Visualization Tools
-----------------------------

In addition to the built-in scripts, FLUXOS output can be visualized with:

* `ParaView <https://www.paraview.org/>`_ -- Open-source 3D visualization tool. Load ``.pvd`` time series files directly for triangular mesh output, or convert regular mesh output via ``vtk_generator.py``.
* `VisIt <https://wci.llnl.gov/simulation/computer-codes/visit>`_ -- Open-source interactive visualization and analysis tool from LLNL.
* `QGIS <https://qgis.org/>`_ -- Open-source GIS for spatial analysis and mapping of simulation results.
* `Google Earth <https://earth.google.com/>`_ -- Overlay inundation maps using KMZ files generated by ``fluxos_viewer.py`` (recommended) or the legacy ``googleearth_klmgen.py``.
