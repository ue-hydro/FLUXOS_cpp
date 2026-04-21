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

Model Configuration Template (``1_Model_Config/model_config_template.py``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The model configuration template replaces the old ``fluxos_setup.py`` CLI. A single editable Python file drives DEM preparation, Gmsh mesh generation, and ``modset.json`` creation in one run — and produces a self-contained HTML report that summarises the build and lists copy-paste Docker commands for running the model.

.. code-block:: bash

   cd fluxos_preprocessing/1_Model_Config
   # open model_config_template.py in your editor, fill in the _config dict, then:
   python model_config_template.py

The template is organised into ten labelled sections (project metadata, paths, DEM source, mesh, forcing files, simulation settings, output, optional modules, Docker/run, report options). A typical edit for a new project only needs the project name, the GeoTIFF path, the mesh type, and the modset filename.

**What it produces:**

.. list-table::
   :widths: 40 60
   :header-rows: 1

   * - Output
     - Purpose
   * - ``bin/<dem>.asc``
     - Downscaled DEM in ESRI ASCII (or a 10×10 dummy ASC when ``mesh_type = "triangular"``, since elevations then live in the ``.msh`` vertex z-coordinates)
   * - ``bin/<mesh>.msh``
     - Slope-adaptive triangular mesh (triangular path only)
   * - ``bin/<modset>.json``
     - FLUXOS configuration consumed by the C++ binary
   * - ``1_Model_Config/reports/<project>_report.html``
     - Self-contained HTML summary with KPI tiles, domain/mesh statistics, modset table, module status, and a *Next Steps* section with OS-aware Docker / viewer snippets

**Internal layout (don't edit):**

The template delegates to three utility modules inside ``config_support_lib/``:

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Module
     - Responsibility
   * - ``dem_operations.py``
     - GeoTIFF reading (``rasterio``), bilinear downscaling (``scipy.ndimage.zoom``), ESRI ASCII writing
   * - ``mesh_generation.py``
     - Slope-based size field + boundary extraction + Gmsh triangular mesh generation (``pygmsh``/``gmsh``); embeds DEM elevations in vertex z-coordinates
   * - ``driver.py``
     - Orchestrator that calls the preprocessing modules, writes ``modset.json`` (in the format read by ``read_modset()`` in ``src/read_functions.cpp``), captures metadata, and invokes the report generator
   * - ``gen_report.py``
     - Builds the inline-styled HTML report (no Jinja, no external assets beyond CDN fonts)

See :doc:`Digital` for details on the DEM preprocessing pipeline.

External Visualization Tools
-----------------------------

In addition to the built-in scripts, FLUXOS output can be visualized with:

* `ParaView <https://www.paraview.org/>`_ -- Open-source 3D visualization tool. Load ``.pvd`` time series files directly for triangular mesh output.
* `VisIt <https://wci.llnl.gov/simulation/computer-codes/visit>`_ -- Open-source interactive visualization and analysis tool from LLNL.
* `QGIS <https://qgis.org/>`_ -- Open-source GIS for spatial analysis and mapping of simulation results.
* `Google Earth <https://earth.google.com/>`_ -- Overlay inundation maps using KMZ files generated by ``fluxos_viewer.py``.
