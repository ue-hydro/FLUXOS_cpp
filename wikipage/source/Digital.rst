Digital Elevation Model (DEM)
==================================

    The DEM provides detailed information about the topography of the terrain. FLUXOS supports two DEM input workflows:

    1. **ESRI ASCII Grid (.asc)** -- The native C++ format. The solver reads ``.asc`` files directly with header keywords (NCOLS, NROWS, XLLCORNER, YLLCORNER, CELLSIZE, NODATA_VALUE).
    2. **GeoTIFF (.tif)** -- Via the Python preprocessing tool ``fluxos_setup.py``, which converts GeoTIFF to ``.asc`` (for regular mesh) or generates an adaptive triangular mesh with DEM elevations embedded in vertex z-coordinates (for triangular mesh).

    Open source spatial data editors such as SAGA (System for Automated Geoscientific Analyses, http://www.saga-gis.org/en/index.html) and QGIS can also be used for conversion between spatial data formats.

Python Preprocessing Tool (GeoTIFF Support)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``fluxos_preprocessing/`` directory contains Python CLI tools for working with GeoTIFF DEMs and exporting simulation results. The preprocessing tool replaces the need for manual DEM format conversion, and the KML exporter enables animated flood visualization in Google Earth (see :doc:`SupportingScripts` for documentation).

**Installation:**

.. code-block:: bash

   cd fluxos_preprocessing
   pip install -r requirements.txt

**DEM Subcommand -- Inspect, Downscale, and Convert:**

.. code-block:: bash

   # Print GeoTIFF metadata
   python fluxos_setup.py dem --input terrain.tif --info

   # Convert GeoTIFF to ASC at native resolution
   python fluxos_setup.py dem --input terrain.tif --output-asc terrain.asc

   # Downscale DEM to 10m resolution and export to ASC (for regular mesh)
   python fluxos_setup.py dem --input terrain.tif --output-asc terrain_10m.asc --resolution 10

**Mesh Subcommand -- Adaptive Triangular Mesh Generation:**

.. code-block:: bash

   # Generate adaptive triangular mesh from DEM (slope-based refinement)
   python fluxos_setup.py mesh --input terrain.tif --output domain.msh \
       --min-size 5.0 --max-size 50.0 --slope-factor 1.0

The mesh generator computes slope from the DEM and creates finer elements in steep areas and coarser elements in flat areas. DEM elevations are embedded as z-coordinates in the ``.msh`` vertex data, so the C++ solver can read elevations directly from the mesh file.

**Config Subcommand -- Generate modset.json:**

.. code-block:: bash

   # Generate FLUXOS JSON config for regular mesh
   python fluxos_setup.py config --dem-file terrain_10m.asc --mesh-type regular \
       --meteo-file forcing.fluxos --roughness 0.005 --output modset.json

   # Generate config for triangular mesh
   python fluxos_setup.py config --dem-file domain_dem.asc --mesh-type triangular \
       --mesh-file domain.msh --meteo-file forcing.fluxos --output modset.json

.. note::

   **Regular Mesh:** The Python tool converts GeoTIFF to ESRI ASCII ``.asc``. The C++ solver reads ``.asc`` files unchanged.

   **Triangular Mesh:** The Python tool generates a Gmsh ``.msh`` file with DEM elevations in vertex z-coordinates. The C++ solver reads vertex z values and computes cell bed elevation as the mean of its three vertices. If the mesh file has z=0 (e.g., manually created mesh), the solver falls back to interpolating from the DEM grid as before. A companion ``.asc`` file is also generated for the ``DEM_FILE`` config entry.

ESRI ASCII Grid Format
^^^^^^^^^^^^^^^^^^^^^^^^

The ESRI ASCII grid format has the following structure. The header keywords can be in upper or lowercase (e.g., NROWS, nrows):

.. code-block:: text

   NCOLS 859
   NROWS 618
   XLLCORNER 425430.000000
   YLLCORNER 5784884.000000
   CELLSIZE 2.000000
   NODATA_VALUE -9999
   <elevation data rows, north to south>

