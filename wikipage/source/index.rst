.. FLUXOS documentation master file, created by
   sphinx-quickstart on Tue Feb 28 07:27:26 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: uevora_logo.png
   :width: 350px
   :align: center
   :alt: Universidade de Évora

|

Welcome to FLUXOS Documentation
=========================================

**FLUXOS** is a high-performance, physics-based hydrodynamic and solute transport model designed for basin-scale environmental simulations.

Key Features
------------

* **2D Shallow Water Equations** - Roe's approximate Riemann solver with MUSCL reconstruction
* **Unstructured Triangular Mesh** - Edge-based finite volume solver with rotated Roe/HLL flux, MUSCL + Barth-Jespersen limiter, and hydrostatic reconstruction
* **Solute Transport** - Advection-Dispersion-Reaction equation for contaminant tracking
* **Adaptive Time Stepping** - CFL-based automatic time step control
* **CUDA GPU Acceleration** - Full GPU offloading for both regular and triangular mesh solvers
* **Parallel Computing** - Hybrid MPI+OpenMP for HPC clusters
* **Wetting/Drying** - Robust handling of dynamic wet/dry fronts
* **GeoTIFF Preprocessing** - Python CLI tool for DEM import, downscaling, and adaptive mesh generation
* **Google Earth Export** - KMZ exporter for animated flood visualization in Google Earth

Quick Start
-----------

.. code-block:: bash

   # Clone repository
   git clone https://github.com/DiogoCostaPT/FLUXOS_cpp.git
   cd FLUXOS_cpp

   # Build (OpenMP only, regular mesh)
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

For GPU acceleration with CUDA:

.. code-block:: bash

   cmake -DMODE_release=ON -DUSE_CUDA=ON ..
   make -j$(nproc)

For triangular mesh support (with optional GPU):

.. code-block:: bash

   cmake -DMODE_release=ON -DUSE_TRIMESH=ON -DUSE_CUDA=ON ..
   make -j$(nproc)

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
