.. FLUXOS documentation master file, created by
   sphinx-quickstart on Tue Feb 28 07:27:26 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to FLUXOS-OVERLAND Documentation
=========================================

**FLUXOS-OVERLAND** is a high-performance, physics-based hydrodynamic and solute transport model designed for basin-scale environmental simulations.

Key Features
------------

* **2D Shallow Water Equations** - Roe's approximate Riemann solver with MUSCL reconstruction
* **Solute Transport** - Advection-Dispersion-Reaction equation for contaminant tracking
* **Adaptive Time Stepping** - CFL-based automatic time step control
* **Parallel Computing** - Hybrid MPI+OpenMP for HPC clusters
* **Wetting/Drying** - Robust handling of dynamic wet/dry fronts

Quick Start
-----------

.. code-block:: bash

   # Clone repository
   git clone https://github.com/DiogoCostaPT/FLUXOS_cpp.git
   cd FLUXOS_cpp

   # Build (OpenMP only)
   mkdir build && cd build
   cmake -DMODE_release=ON ..
   make -j$(nproc)

   # Run
   ./bin/fluxos ../input/modset.json

For HPC clusters with MPI:

.. code-block:: bash

   cmake -DMODE_release=ON -DUSE_MPI=ON ..
   make -j$(nproc)
   mpirun -np 4 ./bin/fluxos_mpi ../input/modset.json

Contents
--------

.. toctree::
   :maxdepth: 4

   Home
   Documentation
   Installation
   Input
   Developer
   Reference

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
