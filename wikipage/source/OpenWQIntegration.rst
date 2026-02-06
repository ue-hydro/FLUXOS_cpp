OpenWQ Integration (Optional)
==================================

.. note::

   This section is under development.

Overview
--------

FLUXOS can be coupled with `OpenWQ <https://github.com/ue-hydro/openwq>`_ (Open Water Quality), an open-source biogeochemical cycling framework for integrated water quality modelling. This coupling enables FLUXOS to simulate not only hydrodynamics but also the transport and transformation of multiple chemical species across the landscape.

The coupling is activated via the ``EXTERNAL_MODULES`` section in the master JSON configuration file:

.. code-block:: json

   {
       "EXTERNAL_MODULES": {
           "ADE_TRANSPORT": 1,
           "OPENWQ": 1
       }
   }

Key Capabilities
-----------------

* **Multi-species transport**: Track multiple dissolved and particulate species simultaneously
* **Biogeochemical reactions**: Define reaction networks including decay, sorption, nutrient cycling, and custom kinetics
* **Source/sink terms**: Spatially distributed source and sink terms for pollutant loading
* **Integrated workflow**: OpenWQ reads the hydrodynamic state (water depth, velocity, discharge) from FLUXOS at each timestep and returns updated concentrations

Build Requirements
-------------------

To build FLUXOS with OpenWQ support:

.. code-block:: bash

   cmake -DMODE_release=ON -DUSE_OPENWQ=ON ..
   make -j$(nproc)

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

Multiple species concentrations can be exported for Google Earth visualization:

.. code-block:: bash

   python fluxos_viewer.py --results-dir ./Results --dem ./terrain.asc --variable conc_SW

See :doc:`SupportingScripts` for more visualization tools.
