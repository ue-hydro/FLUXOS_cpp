Digital Elevation Model (DEM)
==================================

    The DEM will provide detailed information about the topography of the terrain. FLUXOS-OVERLAND reads DEMs in ASCII-ESRI Arc/info Grid format. Open source spatial data editors such as SAGA (System for Automated Geoscientific Analyses, http://www.saga-gis.org/en/index.html) can be used for convertion between spatial data formats. The header keywords can by in upper or lowercase (e.g., NROWS, nrows)

.. note::

   **Triangular Mesh:** When using unstructured triangular meshes (``MESH_TYPE: "triangular"``), the DEM file is still required. The bed elevation at each triangular cell centroid is interpolated from the DEM grid during initialization. The DEM provides the topographic data, while the mesh file (Gmsh ``.msh`` or Triangle ``.node/.ele``) provides the mesh topology. See the :doc:`TriangularMesh` page for details.

