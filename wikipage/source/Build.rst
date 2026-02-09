Building Requirements
==================================

This section describes the software dependencies and build requirements for FLUXOS.

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
   * - CUDA Toolkit
     - 11.0+
     - Required for GPU acceleration (NVIDIA GPUs, Compute Capability 6.0+)
   * - METIS
     - 5.0+
     - Optional graph partitioning for triangular mesh MPI decomposition
   * - HDF5
     - 1.10+
     - Optional parallel I/O support
   * - LAPACK/BLAS
     - Any
     - High-performance linear algebra (Armadillo backend)

Build Modes
-----------

FLUXOS supports several build configurations:

**Standard Build (OpenMP only)**

.. code-block:: bash

   cmake -DMODE_release=ON ..
   make

**MPI+OpenMP Hybrid Build**

.. code-block:: bash

   cmake -DMODE_release=ON -DUSE_MPI=ON ..
   make

**CUDA GPU Build**

.. code-block:: bash

   cmake -DMODE_release=ON -DUSE_CUDA=ON ..
   make

**Triangular Mesh Build**

.. code-block:: bash

   cmake -DMODE_release=ON -DUSE_TRIMESH=ON ..
   make

**Full-Feature Build (Triangular Mesh + GPU + MPI)**

.. code-block:: bash

   cmake -DMODE_release=ON -DUSE_TRIMESH=ON -DUSE_CUDA=ON -DUSE_MPI=ON ..
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
   * - USE_CUDA
     - OFF
     - Enable CUDA GPU acceleration (requires NVIDIA GPU and CUDA Toolkit 11.0+)
   * - USE_TRIMESH
     - OFF
     - Enable unstructured triangular mesh support (Gmsh/Triangle mesh formats)
   * - USE_OPENWQ
     - OFF
     - Enable OpenWQ water quality coupling

.. note::

   ``USE_CUDA`` and ``USE_TRIMESH`` can be combined. When both are enabled, the triangular mesh solver uses GPU acceleration with 7 specialized CUDA kernels.

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

Running the Example
-------------------

FLUXOS includes a test case in the ``bin/`` directory (Rosa Creek watershed, 859x618 cells at 2m resolution).

**Regular Mesh:**

.. code-block:: bash

   mkdir -p Results
   ./build/bin/fluxos bin/modset.json

**Triangular Mesh:**

First generate the mesh from the DEM:

.. code-block:: bash

   python3 fluxos_preprocessing/fluxos_setup.py mesh \
       --input bin/Rosa_2m.asc --output bin/Rosa_trimesh.msh \
       --min-size 10 --max-size 50 --slope-factor 0.5

Then run with the triangular mesh config:

.. code-block:: bash

   mkdir -p Results
   ./build/bin/fluxos bin/modset_trimesh.json

Visualizing Results in Google Earth
-------------------------------------

Export simulation results as KMZ files for animated visualization:

.. code-block:: bash

   # Regular mesh results
   python fluxos_preprocessing/fluxos_viewer.py \
       --results-dir Results --dem bin/Rosa_2m.asc --utm-zone 10

   # Triangular mesh results
   python fluxos_preprocessing/fluxos_viewer.py \
       --results-dir Results --dem bin/Rosa_2m.asc \
       --mesh-type triangular --utm-zone 10

   # Open in Google Earth
   open fluxos_regular.kmz    # macOS
   xdg-open fluxos_regular.kmz  # Linux

Use the time slider in Google Earth to animate through simulation timesteps.

Benchmark Results
-----------------

Tested on the Rosa Creek example (859x618 grid, 5h simulation, 1h output steps) on Apple M-series:

.. list-table::
   :widths: 30 25 25 20
   :header-rows: 1

   * - Configuration
     - Mesh Type
     - Wall Time
     - Output Size
   * - OpenMP (release)
     - Regular (530K cells)
     - 4.73 s
     - 210 MB (6 x 35 MB .txt)
   * - OpenMP (release)
     - Triangular (4528 cells)
     - 0.85 s
     - 9.2 MB (6 x 1.5 MB .vtu)

.. note::

   CUDA acceleration requires an NVIDIA GPU (not available on macOS ARM).
   MPI domain decomposition is available for distributed computing on HPC clusters.

Verifying the Build
-------------------

After building, verify the executable:

.. code-block:: bash

   # Check executable exists
   ls -la build/bin/fluxos

   # Run with the example case
   ./build/bin/fluxos bin/modset.json

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

**CUDA not found:**

Ensure CUDA Toolkit is installed and ``nvcc`` is in your PATH:

.. code-block:: bash

   # Check CUDA installation
   nvcc --version
   nvidia-smi

   # Set CUDA path if needed
   export CUDA_HOME=/usr/local/cuda
   export PATH=$CUDA_HOME/bin:$PATH

**CUDA compute capability mismatch:**

If you get architecture-related errors, specify your GPU's compute capability:

.. code-block:: bash

   cmake -DMODE_release=ON -DUSE_CUDA=ON -DCUDA_ARCH=75 ..  # For RTX 2080
   cmake -DMODE_release=ON -DUSE_CUDA=ON -DCUDA_ARCH=86 ..  # For RTX 3090
