Output files
==================================

Regular Mesh Output (ASCII Text)
---------------------------------

    The output files are written to the "Results" directory. If non-existant, this folder is created during cmake run. The results are printed at the timestep defined in the Master Configuration file (variable PRINT_STEP, in seconds, e.g., 1h = 3600 seconds). The results are saved in ASCII files; one file for each time step (with the file name being the time in seconds). The basic time unit in the model is seconds. The results are printed in each column are (in order):

+---------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+
|irow [-]             |row number of cell                                                                                                                                         |
+---------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+
|icol[-]              |column number of cell                                                                                                                                      |
+---------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+
|z [m]                |water elevation above sea level (with reference to the DEM)                                                                                                |
+---------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+
|h [m]                |water depth                                                                                                                                                |
+---------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+
|ux [m/s]             |velocity in the x direction                                                                                                                                |
+---------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+
|uy [m/s]             |velocity in the y direction                                                                                                                                |
+---------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+
|qx * dxy [m3/sec]    |flow in the x direction                                                                                                                                    |
+---------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+
|qy * dxy [m3/sec]    |flow in the y direction                                                                                                                                    |
+---------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+
|us [m/s]             |sheer velocity                                                                                                                                             |
+---------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+
|conc_SW [mg/l]       |concentration of solute in water                                                                                                                           |
+---------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+
|fe_1 [m2/s]          |flow in the x direction                                                                                                                                    |
+---------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+
|fn_1 [m2/s]          |flow in the y direction                                                                                                                                    |
+---------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+
|twetimetracer [sec]  |wet time for each cell                                                                                                                                     |
+---------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+

**Visualizing Regular Mesh Results:**

Results can be visualized using:

1. **3D Flood Viewer** (recommended): The built-in ``fluxos_viewer.py`` tool renders results as an animated flood on 3D terrain. See the `3D Flood Viewer (Interactive)`_ section below.
2. **VTK conversion**: The ``fluxos_supporting_scripts/vtk_generator.py`` script converts ``.txt`` output files to VTK format for visualization in ParaView.
3. **Google Earth**: The ``fluxos_supporting_scripts/googleearth_klmgen.py`` script generates KML files for overlay on Google Earth.
4. **Custom scripting**: The CSV format is easily loaded with Python (``numpy``, ``pandas``) or MATLAB for custom analysis and plotting.

Triangular Mesh Output (VTK XML Unstructured Grid)
---------------------------------------------------

When using the triangular mesh (``MESH_TYPE: "triangular"``), output is written in VTK XML Unstructured Grid format (``.vtu``), which can be visualized in `ParaView <https://www.paraview.org/>`_.

**File Structure:**

.. code-block:: text

   fluxos_out/
   ├── 3600.vtu                    # Timestep at 3600 seconds
   ├── 7200.vtu                    # Timestep at 7200 seconds
   ├── 10800.vtu                   # Timestep at 10800 seconds
   ├── ...
   └── fluxos_timeseries.pvd       # ParaView time series collection

**Cell-Centered Data Arrays:**

Each ``.vtu`` file contains the triangular mesh geometry (vertices and cells) along with the following cell-centered data arrays:

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Array Name
     - Units
     - Description
   * - h
     - m
     - Water depth
   * - z
     - m
     - Water surface elevation
   * - zb
     - m
     - Bed elevation
   * - ux
     - m/s
     - Velocity in x direction
   * - uy
     - m/s
     - Velocity in y direction
   * - velocity_magnitude
     - m/s
     - :math:`\sqrt{u_x^2 + u_y^2}`
   * - qx
     - m2/s
     - Unit discharge in x direction
   * - qy
     - m2/s
     - Unit discharge in y direction
   * - us
     - m/s
     - Shear velocity
   * - conc_SW_0
     - mg/l
     - Solute concentration (species 0)
   * - twetimetracer
     - s
     - Wet time tracer
   * - cell_area
     - m2
     - Cell area

**PVD Time Series File:**

The ``.pvd`` file is a ParaView Data collection file that links all timestep ``.vtu`` files into a time series for animation and temporal analysis:

.. code-block:: xml

   <?xml version="1.0"?>
   <VTKFile type="Collection" version="0.1">
     <Collection>
       <DataSet timestep="3600.0" file="3600.vtu"/>
       <DataSet timestep="7200.0" file="7200.vtu"/>
       ...
     </Collection>
   </VTKFile>

Open the ``.pvd`` file in ParaView to load all timesteps as an animation.

**Viewing Results in ParaView:**

1. Open ParaView and load the ``fluxos_timeseries.pvd`` file
2. Select cell data arrays to visualize (e.g., ``h``, ``velocity_magnitude``)
3. Apply color maps and thresholds for analysis
4. Use the animation controls to step through time

3D Flood Viewer (Interactive)
-----------------------------

FLUXOS includes an interactive 3D visualization tool (``fluxos_viewer.py``) that renders simulation results as an animated flood on top of the DEM terrain. The viewer supports both regular and triangular mesh outputs.

**Installation:**

.. code-block:: bash

   cd fluxos_preprocessing
   pip install -r requirements.txt

**Usage:**

.. code-block:: bash

   # Regular mesh results
   python fluxos_viewer.py --results-dir ./Results --dem ./terrain.asc --mesh-type regular

   # Triangular mesh results
   python fluxos_viewer.py --results-dir ./fluxos_out --dem ./terrain.asc --mesh-type triangular

   # Custom settings
   python fluxos_viewer.py --results-dir ./Results --dem ./terrain.asc \
       --variable velocity --exaggeration 3.0 --water-opacity 0.8 --clim 0 1.5

   # View terrain only (no simulation data)
   python fluxos_viewer.py --results-dir ./Results --dem ./terrain.asc --terrain-only

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
     - ``regular`` (Cartesian .txt) or ``triangular`` (VTK .vtu)
   * - ``--variable``
     - h
     - Variable to color: ``h`` (water depth), ``velocity``, or ``conc_SW``
   * - ``--exaggeration``
     - 5.0
     - Vertical exaggeration factor for 3D relief
   * - ``--water-opacity``
     - 0.7
     - Water surface transparency (0 = transparent, 1 = opaque)
   * - ``--clim``
     - auto
     - Color range for the variable (e.g., ``--clim 0 2.0``)
   * - ``--terrain-only``
     - off
     - Show DEM terrain without simulation results

**Controls:**

* **Mouse**: Left-drag to orbit, middle-drag to pan, scroll to zoom
* **Arrow keys**: Left/Right to step through timesteps
* **p key**: Play/pause automatic animation
* **Slider**: Drag the timeline slider to jump to any timestep
* **q key**: Quit
