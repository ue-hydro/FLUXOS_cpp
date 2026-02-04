Linux Installation
==================================

This guide provides detailed instructions for installing FLUXOS-OVERLAND on Linux systems, including workstations and HPC clusters.

Quick Start
-----------

For users familiar with Linux development:

.. code-block:: bash

   # Install dependencies (Ubuntu/Debian)
   sudo apt-get install cmake libarmadillo-dev libomp-dev

   # Clone and build
   git clone https://github.com/DiogoCostaPT/FLUXOS_cpp.git
   cd FLUXOS_cpp
   mkdir build && cd build
   cmake -DMODE_release=ON ..
   make -j$(nproc)

   # Run
   ./bin/fluxos ../input/modset.json

Detailed Installation
---------------------

1. Armadillo
^^^^^^^^^^^^

Armadillo is a C++ linear algebra library. The version in package managers may be outdated.

**Option A: From Package Manager (if version >= 9.9)**

.. code-block:: bash

   # Ubuntu/Debian
   sudo apt-get install libarmadillo-dev

   # Fedora/RHEL
   sudo dnf install armadillo-devel

   # Check version
   pkg-config --modversion armadillo

**Option B: From Source (recommended for HPC)**

Download from: http://arma.sourceforge.net/download.html

.. code-block:: bash

   # Download and extract
   wget https://sourceforge.net/projects/arma/files/armadillo-12.6.6.tar.xz
   tar -xf armadillo-12.6.6.tar.xz
   cd armadillo-12.6.6

   # Build and install
   cmake .
   make
   sudo make install

   # Update library cache
   sudo ldconfig

2. OpenMP
^^^^^^^^^

OpenMP is typically included with modern compilers (GCC, Clang, Intel).

.. code-block:: bash

   # Ubuntu/Debian (if not already installed)
   sudo apt-get install libomp-dev

   # Verify OpenMP support
   echo | cpp -fopenmp -dM | grep -i openmp

3. CMake
^^^^^^^^

CMake 3.10 or higher is required.

.. code-block:: bash

   # Ubuntu/Debian
   sudo apt-get install cmake

   # Check version
   cmake --version

   # If version is too old, install from source
   wget https://cmake.org/files/v3.26/cmake-3.26.4-linux-x86_64.tar.gz
   tar -xzf cmake-3.26.4-linux-x86_64.tar.gz
   export PATH=$PWD/cmake-3.26.4-linux-x86_64/bin:$PATH

4. MPI (Optional - for HPC)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

For distributed computing on clusters:

.. code-block:: bash

   # Ubuntu/Debian - OpenMPI
   sudo apt-get install openmpi-bin libopenmpi-dev

   # Ubuntu/Debian - MPICH
   sudo apt-get install mpich libmpich-dev

   # Fedora/RHEL
   sudo dnf install openmpi openmpi-devel

   # Verify MPI
   mpirun --version
   which mpicc mpicxx

5. Build FLUXOS-OVERLAND
^^^^^^^^^^^^^^^^^^^^^^^^

**Standard Build (OpenMP only)**

.. code-block:: bash

   cd FLUXOS_cpp
   mkdir build && cd build

   # Configure
   cmake -DMODE_release=ON ..

   # Build (use all available cores)
   make -j$(nproc)

   # Executable location
   ls -la bin/fluxos

**MPI+OpenMP Hybrid Build**

.. code-block:: bash

   cd FLUXOS_cpp
   mkdir build && cd build

   # Configure with MPI
   cmake -DMODE_release=ON -DUSE_MPI=ON ..

   # Build
   make -j$(nproc)

   # Executable location
   ls -la bin/fluxos_mpi

HPC Cluster Installation
------------------------

On HPC systems with environment modules:

.. code-block:: bash

   # Load required modules (adjust for your system)
   module purge
   module load gcc/11.2.0
   module load cmake/3.20
   module load armadillo/11.0
   module load openmpi/4.1.1

   # Clone repository
   git clone https://github.com/DiogoCostaPT/FLUXOS_cpp.git
   cd FLUXOS_cpp

   # Build with MPI support
   mkdir build && cd build
   cmake -DMODE_release=ON -DUSE_MPI=ON ..
   make -j8

Running FLUXOS
--------------

**Serial/OpenMP Execution**

.. code-block:: bash

   # Set number of OpenMP threads
   export OMP_NUM_THREADS=8

   # Run
   ./build/bin/fluxos ./input/modset.json

**MPI Execution (Workstation)**

.. code-block:: bash

   # Run with 4 MPI processes
   mpirun -np 4 ./build/bin/fluxos_mpi ./input/modset.json

**MPI+OpenMP Hybrid Execution**

.. code-block:: bash

   # 4 MPI processes, 2 OpenMP threads each
   export OMP_NUM_THREADS=2
   mpirun -np 4 ./build/bin/fluxos_mpi ./input/modset.json

**HPC Cluster (SLURM)**

.. code-block:: bash

   # Submit job
   sbatch scripts/run_fluxos_hpc.slurm

   # Check status
   squeue -u $USER

   # View output
   tail -f fluxos_*.out

Environment Variables
---------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Variable
     - Description
   * - OMP_NUM_THREADS
     - Number of OpenMP threads per process
   * - OMP_PROC_BIND
     - Thread-to-core binding (close, spread, master)
   * - OMP_PLACES
     - Thread placement (cores, threads, sockets)

Example optimal settings:

.. code-block:: bash

   export OMP_NUM_THREADS=4
   export OMP_PROC_BIND=close
   export OMP_PLACES=cores

Troubleshooting
---------------

**Armadillo headers not found:**

.. code-block:: bash

   # Locate Armadillo
   find /usr -name "armadillo" 2>/dev/null

   # Set path in CMake
   cmake -DARMADILLO_INCLUDE_DIR=/usr/local/include ..

**OpenMP not working:**

.. code-block:: bash

   # Check compiler supports OpenMP
   g++ -fopenmp -o test test.cpp

   # On macOS with Clang
   brew install libomp
   export OpenMP_ROOT=$(brew --prefix)/opt/libomp

**MPI errors at runtime:**

.. code-block:: bash

   # Check MPI installation
   which mpirun mpiexec

   # Verify library paths
   ldd ./build/bin/fluxos_mpi | grep mpi

   # If using modules, ensure same MPI was used for build and run
   module list

**Permission denied:**

.. code-block:: bash

   # Make executable
   chmod +x ./build/bin/fluxos
   chmod +x ./build/bin/fluxos_mpi

**Library not found at runtime:**

.. code-block:: bash

   # Add library path
   export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

   # Update system library cache
   sudo ldconfig

Performance Tips
----------------

1. **Use native optimization**: Build with ``-march=native`` (default in release mode)

2. **Bind threads to cores**: Set ``OMP_PROC_BIND=close`` and ``OMP_PLACES=cores``

3. **Balance MPI/OpenMP**: For hybrid runs, use 2-4 OpenMP threads per MPI process

4. **Monitor resource usage**: Use ``htop`` or ``top`` to verify all cores are utilized

5. **Check memory usage**: Large domains require significant memory per process
