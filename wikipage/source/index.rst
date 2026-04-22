.. FLUXOS documentation master file, created by
   sphinx-quickstart on Tue Feb 28 07:27:26 2023.

Welcome to FLUXOS Documentation
=========================================

**FLUXOS** is a high-performance, physics-based hydrodynamic and solute transport model designed for basin-scale environmental simulations. It solves the 2-D shallow water equations on either a regular Cartesian mesh or an unstructured triangular mesh, with optional advection–dispersion solute transport and soil-infiltration loss.

Key Features
------------

* **2D Shallow Water Equations** — Roe's approximate Riemann solver with MUSCL reconstruction
* **Unstructured Triangular Mesh** — edge-based finite volume solver with rotated Roe/HLL flux, MUSCL + Barth–Jespersen limiter, and hydrostatic reconstruction
* **Solute Transport** — advection–dispersion–reaction equation for contaminant tracking
* **Adaptive Time Stepping** — CFL-based automatic time-step control
* **CUDA GPU Acceleration** — full GPU offloading for both regular and triangular mesh solvers
* **Parallel Computing** — hybrid MPI+OpenMP for HPC clusters
* **Wetting/Drying** — robust handling of dynamic wet/dry fronts
* **Template-driven preprocessing** — one editable Python template produces the DEM (``.asc``), Gmsh mesh (``.msh``), and ``modset.json`` in a single run, plus an HTML configuration report (``supporting_scripts/1_Model_Config/``)
* **Template-driven post-processing** — a second template streams simulation output into an HTML flood-statistics report (time-series, max-inundation map, ARR-2019 hazard, depth histogram); the same folder hosts ``fluxos_viewer.py`` for KML / WebGL / MP4 exports (``supporting_scripts/2_Read_Outputs/``)
* **Containers** — ships Docker and Apptainer recipes so newcomers and HPC users can skip system-level dependency wrangling

Recommended workflow
--------------------

The fastest path from clone to a running simulation is the container workflow.
All paths below assume you have cloned the repository:

.. code-block:: bash

   git clone https://github.com/ue-hydro/FLUXOS_cpp.git
   cd FLUXOS_cpp

**1. Build and run inside a container**

.. code-block:: bash

   # (one-time) build the FLUXOS dependency image
   docker compose -f containers/docker-compose.yml build

   # open an interactive shell in the container
   docker compose -f containers/docker-compose.yml run --rm fluxos

   # inside the container shell — /work is the bind-mounted repo root:
   cd /work && mkdir -p build && cd build
   cmake -DMODE_release=ON -DUSE_TRIMESH=ON \
         -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=/work/bin /work
   make -j$(nproc)

   cd /work
   ./bin/fluxos Working_example/modset_trimesh.json

The compiled binary lands at ``bin/fluxos`` on the host (bind-mounted from
``/work/bin``), and simulation outputs are written to ``Results/``. See
:doc:`Containers` for the full Docker / Apptainer guide and
:doc:`Installation` for a native build.

**2. Build your own input from a GeoTIFF DEM**

The ``Working_example/modset_trimesh.json`` above is a canned Rosa Creek
example. For your own domain, install the Python dependencies once:

.. code-block:: bash

   python3 -m venv .venv
   source .venv/bin/activate          # Windows: .venv\Scripts\Activate.ps1
   pip install -r supporting_scripts/requirements.txt

Then edit the ``_config`` dict at the top of
``supporting_scripts/1_Model_Config/model_config_template.py`` (project
name, DEM path or OpenTopography/USGS download bbox, mesh type, forcing
files, optional modules, …) and run it:

.. code-block:: bash

   python supporting_scripts/1_Model_Config/model_config_template.py

The template writes the ``.asc`` DEM, the Gmsh ``.msh`` (if
``mesh_type="triangular"``), and the ``modset.json`` into ``Working_example/``,
then opens an HTML configuration report with copy-paste Docker/Apptainer
commands pre-filled for your project. See :doc:`SupportingScripts` for the
full reference.

**3. Analyse the results**

Edit the ``_config`` dict in
``supporting_scripts/2_Read_Outputs/read_output_template.py`` (point it at
your ``Results*/`` folder and the modset that drove the run), then:

.. code-block:: bash

   python supporting_scripts/2_Read_Outputs/read_output_template.py

The output is an HTML flood-statistics report (volume / flooded-area
time-series, max-inundation map, ARR-2019 hazard classification, depth
histogram, first-inundation map) and, optionally, KML / WebGL / MP4 exports
via ``fluxos_viewer.py``.

Directory layout
----------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Path
     - Purpose
   * - ``src/fluxos/``
     - C++ solver (regular + triangular mesh paths)
   * - ``containers/``
     - Dockerfile + Apptainer recipes
   * - ``Working_example/``
     - Canonical inputs (DEMs, meshes, forcing, modset JSONs)
   * - ``bin/``
     - Compiled binary lands here (git-ignored except ``README.md``)
   * - ``Results/``
     - Simulation outputs (git-ignored)
   * - ``supporting_scripts/1_Model_Config/``
     - Python template for preprocessing + HTML config report
   * - ``supporting_scripts/2_Read_Outputs/``
     - Python template for post-processing + HTML results report, plus ``fluxos_viewer.py``
   * - ``in_a_nutshell/``
     - Interactive HTML slide decks (overview, compile & run, supporting scripts)
   * - ``wikipage/source/``
     - This Sphinx documentation

Contents
--------

.. toctree::
   :maxdepth: 4

   Home
   Documentation
   Installation
   Input
   SupportingScripts
   Calibration
   Developer
   Reference

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
