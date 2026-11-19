Master configuration file
==================================

The master configuration is JSON file that provides FLUXOS with information and instructions to run a simulation. It contains pairs of variable/command names and values. The variables/commands needed to run a simulation are:

.. list-table:: Variables/Keywords in MASTER file (JSON type)

    *   - Variable/Keyword
        - any useful descriptive information about the project
    *   - ``DEM_FILE``
        - full path to the ESRI ASCII DEM file (``.asc``). Can be generated from GeoTIFF using the Python preprocessing tool (see :doc:`Digital`)
    *   - ``SIM_DATETIME_START``
        - Simulation start datetime
    *   - ``RESTART``
        - Restart option. If ``TRUE``, the model will load the last output file and restart from there.
    *   - ``PRINT_STEP``
        - temporal resolution of the output files (this is not the model timestep as that is defined by the Courant–Friedrichs–Lewy Condition of numerical stability)
    *   - ``ROUGNESS_HEIGHT``
        - average roughness of the terrain measured in height (m), it is used in friction model and water will be retained until this height
    *   - ``METEO_FILE`` (optional)
        - full path to the snowmelt/rainfall timeseries
    *   - ``INFLOW_FILE`` (optional)
        - full pa/rainfall timeseries
    *   - ``COURANT_NUMBER`` (optional)
        - Courant number for the adaptive time step. Default ``0.5`` on the regular mesh (the unsplit 2-D scheme needs <= 0.5; the legacy value 1 caused grid-scale oscillations at bed steps such as channel banks) and ``1`` on the triangular mesh. Range (0, 1].
    *   - ``MANNING_N`` (optional, regular mesh)
        - Manning roughness coefficient n (s m\ :sup:`-1/3`) for bed friction, e.g. 0.03 for a natural channel and 0.04-0.08 for floodplains. Absent or ``0`` = no bed friction (legacy behaviour). See below.
    *   - ``EXTERNAL_MODULES`` => ``ADE_TRANSPORT`` => ``STATUS``
        - activate 2D advection-dispersion solver
    *   - ``EXTERNAL_MODULES`` => ``ADE_TRANSPORT`` => ``D_COEF``
        - dispersion coefficient

Bed Friction — regular mesh (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Without ``MANNING_N`` the regular-mesh solver has **no bed friction**:
``ROUGNESS_HEIGHT`` only sets the depth below which a cell counts as dry and a
floor on the depth used to compute velocities. Setting ``"MANNING_N": 0.03``
adds the Manning bed-shear term

.. math::

   \frac{\partial \mathbf{q}}{\partial t} = \dots - g\,n^2\,\frac{|\mathbf{q}|\,\mathbf{q}}{h^{7/3}}

applied point-implicitly every step (``q_new = q / (1 + dt g n² |q| / h^{7/3})``),
which is unconditionally stable and never reverses the flow. It also sets the
shear velocity ``u* = √g n |u| / h^{1/6}``, which drives the existing eddy-viscosity
term. With ``MANNING_N`` the momentum smoothing uses the face discharge rather than
the legacy face volume (``fe1·dt``), whose effect depended on the time step. One
uniform value applies to the whole domain. The CUDA build refuses the key, and
``MESH_TYPE: "triangular"`` rejects it because that mesh uses ``ROUGNESS_HEIGHT``
as its Manning n.

Test case: ``Working_example/manning_friction_test`` (normal depth reproduced
within ±3.5 %).

Water-Level (Stage) Boundary — regular mesh (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Imposes a time-varying water surface elevation (e.g. an estuary tide or a
sea-level forecast) on part or all of a domain edge. When the key is absent
or ``STATUS`` is ``false`` the model behaves exactly as before.

.. code-block:: json

    "BOUNDARY_CONDITIONS": {
        "WATER_LEVEL": {
            "STATUS": true,
            "FILENAME": "sado_mouth_level.csv",
            "EDGE": "west",
            "VELOCITY": "zero_gradient"
        }
    }

.. list-table:: ``BOUNDARY_CONDITIONS`` => ``WATER_LEVEL``

    *   - ``STATUS``
        - ``true`` to enable
    *   - ``FILENAME``
        - CSV with ``t_seconds,level_m`` rows (a header line is allowed). Same time base as ``INFLOW_FILE`` (seconds from the simulation start); linear interpolation in time, first/last value held outside the series. Levels are in the DEM's vertical datum.
    *   - ``EDGE``
        - ``"west"``, ``"east"``, ``"south"``, ``"north"``; a list of edges (``["west", "south"]``); or a list of segments ``{"EDGE": "west", "FROM": 10, "TO": 40}`` (1-based cell index along the edge: rows counted from the south for west/east, columns from the west for south/north) or ``{"EDGE": "west", "FROM_COORD": 4285000, "TO_COORD": 4286200}`` (northing for west/east, easting for south/north)
    *   - ``VELOCITY`` (optional)
        - ``"zero_gradient"`` (default: velocity copied from the first interior cell) or ``"zero"``
    *   - ``COMMENT`` (optional)
        - free text; any other key is rejected

Every time step, after the inflow/meteo forcing, the boundary cells get
``h = max(level(t) - zb, 0)``. They exchange water with the interior through
the normal cell-face fluxes (the legacy edge weir is switched off on them), so
high water flows into the domain and the domain drains into the boundary when
the level is lower. NODATA cells on the edge are skipped. The net volume the
boundary adds (+) or removes (−) is written at every print step to
``<OUTPUT_FOLDER>/water_level_bc_balance.csv``
(``t_s,level_m,vol_in_m3,vol_out_m3,vol_net_m3``). Use ``vol_net_m3`` in budgets;
the gross columns include exchange between neighbouring boundary cells.

Limitations: regular mesh only; the CUDA build and ``MESH_TYPE: "triangular"``
refuse the option with an error; MPI builds are untested. ADE concentration on
boundary cells is copied from the interior (no external solute input).
``INFLOW_FILE`` or ``METEO_FILE`` is still required and sets the simulation
length.

.. note::

   The DEM reader treats an elevation of exactly ``0`` as NODATA and takes the
   absolute value of negative elevations. For estuaries with beds below the
   datum, shift the DEM and the level series by the same positive offset.

Test case: ``Working_example/water_level_bc_test`` (see its ``README.md``).

Triangular Mesh Configuration (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When built with ``USE_TRIMESH``, the following additional variables are available for configuring unstructured triangular mesh simulations:

.. list-table:: Triangular Mesh Variables/Keywords

    *   - ``MESH_TYPE``
        - Mesh type selector: ``"regular"`` (default) or ``"triangular"``
    *   - ``MESH_FILE``
        - Full path to the triangular mesh file (e.g., ``"domain.msh"`` for Gmsh format)
    *   - ``MESH_FORMAT``
        - Mesh file format: ``"gmsh"`` (Gmsh .msh v2.2) or ``"triangle"`` (Triangle .node/.ele)
    *   - ``BOUNDARY_CONDITIONS``
        - JSON object mapping physical group tags to boundary condition types

The JSON file supports C/C++ syntax for comments: single-line comment (``//``) or comment blocks (``/*`` and ``*/``).

Example (Regular Mesh):

.. code-block:: json

    {
        "COMMNET": "Batch_1 - additional tests for paper",
        "DEM_FILE": "Rosa_2m.asc",
        "SIM_DATETIME_START": "2017SEP01 12:15:00",
        "RESTART": false,

        "INFLOW_FILE": {
            "FILENAME": "Flow_forcing.csv",
            "DISCHARGE_LOCATION":{
                "X_COORDINATE": 425894,
                "Y_COORDINATE": 5785553
            }
        },

        "EXTERNAL_MODULES":{
            "ADE_TRANSPORT":{
                "STATUS": true,
                "D_COEF": 0.01
            }
        },
        "METEO_FILE": "Qmelt_synthetic.fluxos",
        "OUTPUT": {
            "OUTPUT_FOLDER": "../fluxos_out",
            "PRINT_STEP": 3600,
            "H_MIN_TO_PRINT": 0.005
        },
        "ROUGNESS_HEIGHT": 0.005

    }

Example (Triangular Mesh):

.. code-block:: json

    {
        "COMMNET": "Triangular mesh simulation",
        "DEM_FILE": "domain_dem.asc",
        "SIM_DATETIME_START": "2017SEP01 12:15:00",
        "RESTART": false,

        "MESH_TYPE": "triangular",
        "MESH_FILE": "domain.msh",
        "MESH_FORMAT": "gmsh",
        "BOUNDARY_CONDITIONS": {
            "1": {"type": "wall"},
            "2": {"type": "outflow"}
        },

        "INFLOW_FILE": {
            "FILENAME": "Flow_forcing.csv",
            "DISCHARGE_LOCATION":{
                "X_COORDINATE": 425894,
                "Y_COORDINATE": 5785553
            }
        },

        "EXTERNAL_MODULES":{
            "ADE_TRANSPORT":{
                "STATUS": true,
                "D_COEF": 0.01
            }
        },
        "METEO_FILE": "Qmelt_synthetic.fluxos",
        "OUTPUT": {
            "OUTPUT_FOLDER": "../fluxos_out",
            "PRINT_STEP": 3600,
            "H_MIN_TO_PRINT": 0.005
        },
        "ROUGNESS_HEIGHT": 0.005
    }

.. note::

   When ``MESH_TYPE`` is set to ``"triangular"``, the DEM file (``DEM_FILE``) is still required. If the ``.msh`` file contains non-zero vertex z-coordinates (generated by the Python preprocessing tool), the solver uses vertex elevations directly -- computing cell bed elevation as the mean of its three vertex z values. Otherwise, the solver falls back to interpolating bed elevation from the DEM grid. The ``MESH_FILE`` provides the mesh topology, while boundary conditions are defined by physical group tags (mapped via ``BOUNDARY_CONDITIONS``).

   **Tip:** Use the Python preprocessing tool to generate both the ``.msh`` file (with embedded DEM elevations) and the ``modset.json`` config. See :doc:`Digital` for details.



