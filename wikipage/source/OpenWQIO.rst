OpenWQ I/O
==================================

.. note::

   This section is under development.

Overview
--------

This section describes the input and output files specific to the OpenWQ biogeochemical module when coupled with FLUXOS.

OpenWQ Input Files
-------------------

When OpenWQ is enabled (``"OPENWQ": 1`` in the master JSON configuration), additional input files are required to define the biogeochemical model:

* **Species definitions**: List of chemical species to be tracked (e.g., nitrogen, phosphorus, dissolved organic carbon)
* **Reaction network**: Kinetic and equilibrium reactions between species
* **Initial concentrations**: Spatially distributed initial conditions for each species
* **Boundary conditions**: Concentration values at inflow boundaries and atmospheric deposition

OpenWQ Output
--------------

OpenWQ results are included in the standard FLUXOS output:

* **Regular mesh**: The ``conc_SW [mg/l]`` column in the ``.txt`` output files contains the concentration of the primary tracked species
* **Triangular mesh**: The ``conc_SW`` cell data array in ``.vtu`` files contains species concentrations

Multiple species concentrations can be visualized using the 3D Flood Viewer:

.. code-block:: bash

   python fluxos_viewer.py --results-dir ./Results --dem ./terrain.asc --variable conc_SW

See :doc:`CouplingOpenWQ` for details on enabling and configuring the OpenWQ coupling.
