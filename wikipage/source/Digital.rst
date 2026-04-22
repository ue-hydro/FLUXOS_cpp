Digital Elevation Model (DEM)
==================================

    The DEM provides detailed information about the topography of the terrain. FLUXOS supports two DEM input workflows:

    1. **ESRI ASCII Grid (.asc)** -- The native C++ format. The solver reads ``.asc`` files directly with header keywords (NCOLS, NROWS, XLLCORNER, YLLCORNER, CELLSIZE, NODATA_VALUE).
    2. **GeoTIFF (.tif)** -- Via the Python model configuration template (``supporting_scripts/1_Model_Config/model_config_template.py``), which converts GeoTIFF to ``.asc`` (for regular mesh) or generates an adaptive triangular mesh with DEM elevations embedded in vertex z-coordinates (for triangular mesh).

    Open source spatial data editors such as SAGA (System for Automated Geoscientific Analyses, http://www.saga-gis.org/en/index.html) and QGIS can also be used for conversion between spatial data formats.

Python Preprocessing Template (GeoTIFF Support)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``supporting_scripts/1_Model_Config/`` directory contains a template-driven preprocessing workflow that reads a GeoTIFF DEM and produces ``.asc``, ``.msh`` (triangular), and ``modset.json`` files in one run, plus an HTML report with copy-paste Docker commands for running the simulation. See :doc:`SupportingScripts` for the full reference.

**Installation:**

.. code-block:: bash

   cd supporting_scripts
   pip install -r requirements.txt

**Usage:**

Edit the ``_config = dict(...)`` block at the top of ``model_config_template.py`` (DEM path, target resolution, mesh type, mesh sizing, modset name, forcing files, roughness, modules, ...), then run it:

.. code-block:: bash

   cd supporting_scripts/1_Model_Config
   python model_config_template.py

The template handles all three stages internally:

* **DEM preparation:** reads the GeoTIFF via ``rasterio``, optionally downscales via bilinear interpolation (``scipy.ndimage.zoom``), writes an ESRI ASCII ``.asc`` into ``<repo>/bin/``. For triangular meshes it writes a 10×10 placeholder ``.asc`` (elevations actually live in the ``.msh``), which keeps the ``DEM_FILE`` entry valid for the C++ side without duplicating the DEM on disk.
* **Adaptive triangular mesh (triangular path only):** computes slope from the DEM, builds a size field with finer triangles in steep terrain, extracts the valid-data polygon, invokes ``pygmsh`` / ``gmsh`` to generate the mesh, and interpolates DEM elevations onto the vertices so the C++ solver can read elevations directly from the mesh file.
* **modset.json:** writes the configuration with ``DEM_FILE``, ``MESH_FILE`` / ``MESH_FORMAT`` / ``BOUNDARY_CONDITIONS`` (for trimesh), forcing files, sim start, roughness, output options, and any enabled modules (ADE transport, soil infiltration).

.. note::

   **Regular Mesh:** The template converts GeoTIFF to ESRI ASCII ``.asc``. The C++ solver reads ``.asc`` files unchanged.

   **Triangular Mesh:** The template generates a Gmsh ``.msh`` file with DEM elevations in vertex z-coordinates. The C++ solver reads vertex z values and computes cell bed elevation as the mean of its three vertices. If the mesh file has z=0 (e.g., manually created mesh), the solver falls back to interpolating from the DEM grid as before. A companion placeholder ``.asc`` file is also generated for the ``DEM_FILE`` config entry.

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

