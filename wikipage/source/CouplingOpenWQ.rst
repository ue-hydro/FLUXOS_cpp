Coupling OpenWQ
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
