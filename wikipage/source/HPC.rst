High-Performance Computing (HPC)
==================================

This section describes how to build and run FLUXOS-OVERLAND on HPC clusters using MPI+OpenMP hybrid parallelization.

Overview
--------

FLUXOS-OVERLAND supports three parallelization modes:

1. **OpenMP only**: For workstations and small domains
2. **MPI only**: For distributed memory systems
3. **Hybrid MPI+OpenMP**: For maximum scalability on HPC clusters (recommended)

The hybrid approach uses MPI for communication between nodes and OpenMP for parallelism within each node, providing optimal performance for large-scale simulations.

Building for HPC
----------------

Prerequisites
^^^^^^^^^^^^^

* MPI implementation (OpenMPI, MPICH, or Intel MPI)
* OpenMP-enabled compiler (GCC 7+, Intel C++ 18+, or Clang 8+)
* Armadillo 9.9+
* CMake 3.10+
* HDF5 (optional, for parallel output)

Build Commands
^^^^^^^^^^^^^^

.. code-block:: bash

   # Create build directory
   mkdir build && cd build

   # Configure with MPI support
   cmake -DMODE_release=ON -DUSE_MPI=ON ..

   # Build
   make -j8

   # The executable will be: build/bin/fluxos_mpi

For module-based HPC systems:

.. code-block:: bash

   # Load required modules (adjust for your system)
   module load gcc/11.2.0
   module load openmpi/4.1.1
   module load armadillo/11.0
   module load cmake/3.20

   # Build
   mkdir build && cd build
   cmake -DMODE_release=ON -DUSE_MPI=ON ..
   make -j8

Running on HPC Clusters
-----------------------

SLURM Job Script
^^^^^^^^^^^^^^^^

Example SLURM script for running FLUXOS on a cluster:

.. code-block:: bash

   #!/bin/bash
   #SBATCH --job-name=fluxos_mpi
   #SBATCH --output=fluxos_%j.out
   #SBATCH --error=fluxos_%j.err
   #SBATCH --nodes=4
   #SBATCH --ntasks-per-node=32
   #SBATCH --cpus-per-task=2
   #SBATCH --time=24:00:00
   #SBATCH --partition=compute
   #SBATCH --account=your_account

   # Load modules
   module purge
   module load gcc/11.2.0
   module load openmpi/4.1.1
   module load armadillo/11.0

   # Set OpenMP threads per MPI task
   export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
   export OMP_PROC_BIND=close
   export OMP_PLACES=cores

   # Run FLUXOS
   srun --mpi=pmix ./build/bin/fluxos_mpi ./input/modset.json

PBS Job Script
^^^^^^^^^^^^^^

Example PBS script:

.. code-block:: bash

   #!/bin/bash
   #PBS -N fluxos_mpi
   #PBS -l nodes=4:ppn=32
   #PBS -l walltime=24:00:00
   #PBS -q batch

   cd $PBS_O_WORKDIR

   # Load modules
   module load gcc openmpi armadillo

   # Set OpenMP threads
   export OMP_NUM_THREADS=2
   export OMP_PROC_BIND=close

   # Run
   mpirun -np 128 ./build/bin/fluxos_mpi ./input/modset.json

Scalability Guidelines
----------------------

Use this table to select appropriate parallelization for your domain:

.. list-table::
   :widths: 20 20 15 15 30
   :header-rows: 1

   * - Domain Size
     - Mode
     - MPI Procs
     - OMP Threads
     - Notes
   * - < 500x500
     - OpenMP
     - 1
     - 4-8
     - MPI overhead not justified
   * - 500x500 - 2000x2000
     - Hybrid
     - 4-16
     - 2-4
     - Good balance of communication/computation
   * - 2000x2000 - 5000x5000
     - Hybrid
     - 16-64
     - 2-4
     - Scales well to multi-node
   * - > 5000x5000
     - Hybrid
     - 64-256
     - 2-4
     - Large-scale HPC required

Domain Decomposition
--------------------

FLUXOS uses 2D Cartesian domain decomposition:

* The global domain is automatically divided among MPI processes
* MPI_Cart_create establishes the process topology
* Each process computes a local subdomain with ghost cells
* Ghost cells are exchanged at each time step

.. code-block:: text

   Global Domain (1000 x 1000 cells)
   ┌─────────────┬─────────────┐
   │  Process 0  │  Process 1  │
   │  500x500    │  500x500    │
   ├─────────────┼─────────────┤
   │  Process 2  │  Process 3  │
   │  500x500    │  500x500    │
   └─────────────┴─────────────┘

Ghost Cell Exchange
^^^^^^^^^^^^^^^^^^^

Each subdomain maintains ghost cells (halo regions) from neighboring processes:

.. code-block:: text

   ┌─────────────────────────────┐
   │  Ghost cells (from north)   │
   ├───┬───────────────────┬─────┤
   │ G │                   │ G   │
   │ h │   Local domain    │ h   │
   │ o │                   │ o   │
   │ s │                   │ s   │
   │ t │                   │ t   │
   ├───┴───────────────────┴─────┤
   │  Ghost cells (from south)   │
   └─────────────────────────────┘

Performance Optimization
------------------------

OpenMP Settings
^^^^^^^^^^^^^^^

For optimal OpenMP performance:

.. code-block:: bash

   # Bind threads to cores
   export OMP_PROC_BIND=close
   export OMP_PLACES=cores

   # Set number of threads (typically 2-4 for hybrid)
   export OMP_NUM_THREADS=2

MPI Settings
^^^^^^^^^^^^

For OpenMPI:

.. code-block:: bash

   # Disable InfiniBand if not available
   export OMPI_MCA_btl=^openib

   # Use shared memory for intra-node communication
   export OMPI_MCA_btl=vader,self

For Intel MPI:

.. code-block:: bash

   # Use shared memory fabric
   export I_MPI_FABRICS=shm:ofi

Memory Considerations
^^^^^^^^^^^^^^^^^^^^^

* Each MPI process allocates memory for its local subdomain plus ghost cells
* For very large domains, ensure sufficient memory per node
* Consider using fewer MPI processes with more OpenMP threads to reduce memory overhead

Parallel Output
---------------

FLUXOS supports two parallel output modes:

**Gathered Output (Default)**

* All data gathered to root process for writing
* Simpler file format, single output file
* Suitable for moderate domain sizes

**Parallel Output**

* Each process writes its own portion
* Creates multiple files plus a manifest
* Better scalability for very large domains

Output files are named:

.. code-block:: text

   Results/
   ├── 3600.txt              # Gathered output (single file)
   ├── 3600_rank0.txt        # Parallel output (per-process)
   ├── 3600_rank1.txt
   ├── 3600_manifest.txt     # Manifest listing all files
   └── ...

Troubleshooting
---------------

Common Issues
^^^^^^^^^^^^^

**MPI not found during build:**

.. code-block:: bash

   # Ensure MPI is in PATH
   which mpicc
   which mpicxx

   # Set CC and CXX if needed
   export CC=mpicc
   export CXX=mpicxx

**Poor scaling:**

* Check that ghost cell exchange is not dominating
* Ensure load balance (equal subdomain sizes)
* Verify network bandwidth (InfiniBand recommended)

**Memory errors:**

* Reduce number of MPI processes
* Increase OpenMP threads per process
* Check for memory leaks with valgrind

**Segmentation faults:**

* Ensure consistent MPI library between build and run
* Check Armadillo is compiled with same compiler
* Verify input file format

Profiling
^^^^^^^^^

For performance analysis:

.. code-block:: bash

   # With Intel VTune
   vtune -collect hotspots -- srun ./build/bin/fluxos_mpi input.json

   # With Scalasca
   scalasca -analyze srun ./build/bin/fluxos_mpi input.json

   # With ARM MAP
   map --profile srun ./build/bin/fluxos_mpi input.json

Best Practices
--------------

1. **Start small**: Test with few MPI processes before scaling up
2. **Monitor load balance**: Ensure all processes finish at similar times
3. **Use hybrid mode**: Typically 2-4 OpenMP threads per MPI process works best
4. **Check I/O**: Parallel I/O may become a bottleneck for frequent output
5. **Verify results**: Compare small-domain results with serial version
6. **Use restart files**: For long simulations, implement checkpoint/restart
