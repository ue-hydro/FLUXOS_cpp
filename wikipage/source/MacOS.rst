macOS Installation (native)
================================

This guide covers **native** installation on macOS (Intel and Apple Silicon).
For a container-based build — which ships every system dependency
pre-installed and is the easiest path for newcomers — see :doc:`Containers`.

.. note::

   CUDA (``USE_CUDA=ON``) is not supported on macOS because Apple has
   dropped NVIDIA support. For GPU runs, use a Linux workstation or HPC
   cluster with the CUDA Apptainer recipe.

Quick Start (Apple Silicon / Intel, via Homebrew)
--------------------------------------------------

.. code-block:: bash

   # 1. Install Xcode Command Line Tools (one-time)
   xcode-select --install

   # 2. Install Homebrew dependencies
   brew install cmake armadillo libomp nlohmann-json hdf5
   # optional, for MPI runs:
   brew install open-mpi

   # 3. Clone and build
   git clone https://github.com/ue-hydro/FLUXOS_cpp.git
   cd FLUXOS_cpp
   mkdir build && cd build
   cmake -DMODE_release=ON -DUSE_TRIMESH=ON ..
   make -j$(sysctl -n hw.ncpu)

   # 4. Run the bundled Rosa Creek example
   cd ..
   ./build/bin/fluxos Working_example/modset_trimesh.json

Dependency notes
----------------

**Apple Silicon (M1/M2/M3/M4):** use ARM-native Homebrew (installed under
``/opt/homebrew``). Mixing Rosetta / x86_64 Homebrew with an ARM-native
compile breaks linking.

**OpenMP on macOS Clang:** Apple's Clang ships without OpenMP headers.
``brew install libomp`` provides them; CMake's ``FindOpenMP`` picks them up
if Homebrew is on the default search path. If the build reports "OpenMP:
NOT FOUND", prepend the Homebrew OpenMP path explicitly:

.. code-block:: bash

   export OpenMP_ROOT=$(brew --prefix)/opt/libomp
   cmake -DMODE_release=ON -DUSE_TRIMESH=ON ..

**Armadillo:** Homebrew's ``armadillo`` formula is current (12.x). If you
installed an older one, ``brew upgrade armadillo``.

Build variants
--------------

The CMake flag matrix is the same as on Linux — see :doc:`Build` for the
full list. The most common combinations on macOS:

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Goal
     - CMake flags
   * - OpenMP only, regular mesh
     - ``-DMODE_release=ON``
   * - OpenMP + triangular mesh
     - ``-DMODE_release=ON -DUSE_TRIMESH=ON``
   * - MPI + OpenMP
     - ``-DMODE_release=ON -DUSE_MPI=ON``

Running the example
-------------------

After ``make``, the compiled ``fluxos`` binary lives at ``build/bin/fluxos``
inside your checkout. The ``Working_example/`` directory ships with seven
canned modset JSONs that point to the bundled Rosa Creek DEM and forcing
files:

.. code-block:: bash

   ./build/bin/fluxos Working_example/modset_trimesh.json

For your own domain, edit the ``_config`` dict in
``supporting_scripts/1_Model_Config/model_config_template.py`` (see
:doc:`SupportingScripts`) and run it to produce ``.asc`` / ``.msh`` /
``modset.json`` in one shot.

Troubleshooting
---------------

**"clang: error: unsupported option '-fopenmp'"** — Homebrew's ``libomp``
is not installed, or its path isn't on the search path. Install and retry:

.. code-block:: bash

   brew install libomp
   export OpenMP_ROOT=$(brew --prefix)/opt/libomp

**Armadillo link errors on Apple Silicon** — you're likely mixing Rosetta
Homebrew with the ARM-native compiler. Run ``brew config`` and verify the
prefix is ``/opt/homebrew`` (ARM) rather than ``/usr/local`` (x86_64).

**``make`` fails with "operator not supported" in nlohmann/json** — Apple's
stock Clang is too old. Update with ``xcode-select --install`` or install a
newer Clang via ``brew install llvm`` and re-run CMake.
