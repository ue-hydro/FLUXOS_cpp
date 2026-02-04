Building Requirements
==================================

This section describes the software dependencies and build requirements for FLUXOS-OVERLAND.

Required Dependencies
---------------------

.. list-table::
   :widths: 25 25 50
   :header-rows: 1

   * - Dependency
     - Minimum Version
     - Description
   * - C++ Compiler
     - C++11
     - GCC 7+, Clang 8+, or Intel C++ 18+ recommended
   * - CMake
     - 3.10
     - Build system generator
   * - Armadillo
     - 9.9
     - C++ linear algebra library
   * - OpenMP
     - 4.5
     - Shared-memory parallelization (usually bundled with compiler)

Optional Dependencies
---------------------

.. list-table::
   :widths: 25 25 50
   :header-rows: 1

   * - Dependency
     - Version
     - Description
   * - MPI
     - 3.0+
     - Required for distributed computing (OpenMPI, MPICH, or Intel MPI)
   * - HDF5
     - 1.10+
     - Optional parallel I/O support
   * - LAPACK/BLAS
     - Any
     - High-performance linear algebra (Armadillo backend)

Build Modes
-----------

FLUXOS-OVERLAND supports several build configurations:

**Standard Build (OpenMP only)**

.. code-block:: bash

   cmake -DMODE_release=ON ..
   make

**MPI+OpenMP Hybrid Build**

.. code-block:: bash

   cmake -DMODE_release=ON -DUSE_MPI=ON ..
   make

**Debug Build**

.. code-block:: bash

   cmake -DMODE_debug=ON ..
   make

CMake Options
-------------

.. list-table::
   :widths: 30 15 55
   :header-rows: 1

   * - Option
     - Default
     - Description
   * - MODE_release
     - OFF
     - Enable release build with optimizations (-O3)
   * - MODE_debug
     - OFF
     - Enable debug build with symbols (-g)
   * - USE_MPI
     - OFF
     - Enable MPI for distributed computing
   * - USE_OPENWQ
     - OFF
     - Enable OpenWQ water quality coupling

Compiler Optimization Flags
---------------------------

The release build includes aggressive optimization flags for maximum performance:

* ``-O3``: High-level optimization
* ``-march=native``: Optimize for the current CPU architecture
* ``-mtune=native``: Tune for the current CPU
* ``-funroll-loops``: Unroll loops for better performance
* ``-ftree-vectorize``: Enable automatic vectorization
* ``-fno-math-errno``: Disable errno for math functions (faster)
* ``-flto``: Link-time optimization

Platform-Specific Notes
-----------------------

Linux
^^^^^

Most HPC clusters run Linux. Ensure you load the appropriate modules before building:

.. code-block:: bash

   module load gcc/11.2.0
   module load cmake/3.20
   module load armadillo/11.0
   module load openmpi/4.1.1  # if using MPI

macOS
^^^^^

On macOS, install dependencies via Homebrew:

.. code-block:: bash

   brew install cmake armadillo libomp open-mpi

For Apple Silicon (M1/M2), ensure you're using native ARM builds of dependencies.

Windows
^^^^^^^

Windows builds are supported via:

* Visual Studio 2019+ with C++11 support
* MSYS2/MinGW-w64 with GCC
* Windows Subsystem for Linux (WSL) - recommended

Verifying the Build
-------------------

After building, verify the executable:

.. code-block:: bash

   # Check executable exists
   ls -la build/bin/fluxos

   # For MPI build
   ls -la build/bin/fluxos_mpi

   # Run with a test case
   ./build/bin/fluxos ./input/modset.json

   # For MPI (4 processes)
   mpirun -np 4 ./build/bin/fluxos_mpi ./input/modset.json

Troubleshooting
---------------

**Armadillo not found:**

Ensure Armadillo is installed and its include/library paths are accessible:

.. code-block:: bash

   # Check Armadillo installation
   find /usr -name "armadillo" 2>/dev/null

   # Set paths if needed
   cmake -DARMADILLO_INCLUDE_DIR=/path/to/include \
         -DARMADILLO_LIBRARY=/path/to/libarmadillo.so ..

**OpenMP not found:**

For GCC, OpenMP is usually included. For Clang on macOS:

.. code-block:: bash

   brew install libomp
   export OpenMP_ROOT=$(brew --prefix)/opt/libomp

**MPI not found:**

Ensure MPI is in your PATH:

.. code-block:: bash

   which mpicc mpicxx

   # If using environment modules
   module load openmpi

**Link-time optimization (LTO) errors:**

If LTO causes issues, disable it:

.. code-block:: bash

   cmake -DMODE_release=ON -DCMAKE_CXX_FLAGS="-O3 -march=native" ..
