CUDA GPU Acceleration
==================================

This section describes how to build and run FLUXOS-OVERLAND with CUDA GPU acceleration for maximum performance on large-scale simulations.

Overview
--------

FLUXOS-OVERLAND supports CUDA GPU acceleration for both regular Cartesian and unstructured triangular mesh solvers. GPU offloading can provide 10-50x speedup over single-core CPU execution for domains with millions of cells.

**Supported GPU operations:**

* Hydrodynamics solver (Roe flux, state update, wet/dry tracking)
* Courant condition (CFL time step computation)
* ADE solute transport (concentration adjustment)

Requirements
------------

* **NVIDIA GPU**: Compute Capability 6.0+ (Pascal or newer)
* **CUDA Toolkit**: Version 11.0 or later
* **GPU Driver**: Compatible with the installed CUDA Toolkit version

.. list-table:: Recommended GPUs
   :widths: 30 20 25 25
   :header-rows: 1

   * - GPU
     - Compute Capability
     - VRAM
     - Typical Domain Size
   * - GTX 1080 Ti
     - 6.1
     - 11 GB
     - Up to 2000x2000
   * - RTX 2080 Ti
     - 7.5
     - 11 GB
     - Up to 3000x3000
   * - RTX 3090
     - 8.6
     - 24 GB
     - Up to 5000x5000
   * - A100
     - 8.0
     - 40/80 GB
     - 10000x10000+
   * - H100
     - 9.0
     - 80 GB
     - 10000x10000+

Building with CUDA
------------------

.. code-block:: bash

   # Standard CUDA build (regular mesh)
   mkdir build && cd build
   cmake -DMODE_release=ON -DUSE_CUDA=ON ..
   make -j$(nproc)

   # CUDA with triangular mesh
   cmake -DMODE_release=ON -DUSE_CUDA=ON -DUSE_TRIMESH=ON ..
   make -j$(nproc)

   # Full-feature build
   cmake -DMODE_release=ON -DUSE_CUDA=ON -DUSE_TRIMESH=ON -DUSE_MPI=ON ..
   make -j$(nproc)

Running with GPU
----------------

No special command-line flags are needed. When built with ``USE_CUDA``, the GPU solver is used automatically:

.. code-block:: bash

   # Single GPU
   ./bin/fluxos ./input/modset.json

   # With MPI (multi-GPU, one GPU per MPI rank)
   mpirun -np 2 ./bin/fluxos_mpi ./input/modset.json

Regular Mesh GPU Architecture
-----------------------------

The regular mesh GPU solver uses a 2D CUDA thread grid matching the domain layout:

* **Thread block size**: 16x16 (256 threads per block)
* **Grid size**: ``(NCOLS/16, NROWS/16)`` blocks
* Each thread processes one cell at ``(irow, icol)``

**Kernels:**

1. ``cuda_courant_condition``: CFL time step with block-level parallel reduction
2. ``cuda_hydrodynamics_calc``: Roe flux computation + state update
3. ``cuda_ade_adjust``: Concentration adjustment for depth changes

**Memory layout:**

* Device memory mirrors the host ``arma::Mat<double>`` layout
* All fields (h, z, qx, qy, etc.) allocated as contiguous 2D arrays on GPU
* Host-device transfers occur at forcing and output steps

Triangular Mesh GPU Architecture
---------------------------------

The triangular mesh GPU solver uses 1D thread indexing with 7 specialized kernels:

.. code-block:: text

   Per timestep execution order:
   ┌─────────────────────────────────────┐
   │ 1. kernel_tri_wetdry     (cells)    │
   │ 2. kernel_tri_gradient   (cells)    │
   │ 3. kernel_tri_limiter    (cells)    │
   │ 4. kernel_tri_edge_flux  (edges)  ← main compute kernel
   │ 5. kernel_tri_accumulate (edges)    │
   │ 6. kernel_tri_update     (cells)    │
   └─────────────────────────────────────┘
   CFL: kernel_tri_courant    (cells)

**Thread indexing:**

* Cell kernels: 1 thread per cell, ``blockDim = 256``, ``gridDim = (ncells + 255) / 256``
* Edge kernels: 1 thread per edge, ``blockDim = 256``, ``gridDim = (nedges + 255) / 256``

**Race condition handling:**

The flux accumulation kernel (step 5) uses ``atomicAdd`` to safely add edge fluxes to the left and right cells concurrently:

.. code-block:: text

   // Each edge thread atomically adds its flux to both cells
   atomicAdd(&d_dh[left_cell],  -flux_mass * dt / cell_area[left_cell]);
   atomicAdd(&d_dh[right_cell], +flux_mass * dt / cell_area[right_cell]);

**Device memory:**

* Mesh topology (read-only): flat arrays for edge connectivity, normals, cell geometry
* Solution data (read-write): flat ``double*`` arrays indexed by cell/edge ID
* Reduction buffer for CFL computation

**CFL Reduction:**

Block-level parallel reduction computes the global minimum time step:

1. Each thread computes its local CFL candidate
2. Shared memory reduction within each block finds the block minimum
3. Block minimums written to a reduction buffer
4. Host reads and reduces across blocks

Performance Guidelines
----------------------

**Maximizing GPU performance:**

1. **Large domains benefit most**: GPU overhead is amortized over more cells
2. **Minimize host-device transfers**: Only transfer at forcing and output steps
3. **Batch chemical species**: Process all species in sequence on GPU before transferring back
4. **Use pinned memory**: For large domains, pinned (page-locked) host memory improves transfer speed

**Expected speedup factors:**

.. list-table::
   :widths: 25 25 25 25
   :header-rows: 1

   * - Domain Size
     - CPU (OpenMP 8T)
     - GPU (RTX 3090)
     - Speedup
   * - 500x500
     - 1x
     - 3-5x
     - Moderate
   * - 2000x2000
     - 1x
     - 10-20x
     - Good
   * - 5000x5000
     - 1x
     - 20-40x
     - Excellent
   * - 10000x10000
     - 1x
     - 30-50x
     - Excellent

.. note::

   Actual speedup depends on GPU model, domain complexity, and wet fraction. Domains with large dry regions may show lower speedup due to early-exit optimizations in the CPU code.

Multi-GPU with MPI
------------------

For multi-GPU execution, combine ``USE_CUDA`` with ``USE_MPI``. Each MPI rank uses one GPU:

.. code-block:: bash

   # Build
   cmake -DMODE_release=ON -DUSE_CUDA=ON -DUSE_MPI=ON ..
   make -j$(nproc)

   # Run on 4 GPUs
   mpirun -np 4 ./bin/fluxos_mpi ./input/modset.json

**SLURM example for multi-GPU:**

.. code-block:: bash

   #!/bin/bash
   #SBATCH --job-name=fluxos_gpu
   #SBATCH --nodes=2
   #SBATCH --ntasks-per-node=4
   #SBATCH --gres=gpu:4
   #SBATCH --time=12:00:00

   module load cuda/11.8
   module load openmpi/4.1.1

   srun ./build/bin/fluxos_mpi ./input/modset.json

Troubleshooting
---------------

**"CUDA error: no CUDA-capable device is detected":**

.. code-block:: bash

   # Check GPU visibility
   nvidia-smi
   echo $CUDA_VISIBLE_DEVICES

**"CUDA error: out of memory":**

* Reduce domain size or use a GPU with more VRAM
* For triangular meshes: reduce mesh cell count
* Check for memory leaks with ``cuda-memcheck``

**Poor GPU performance:**

* Ensure the GPU is not in power-saving mode: ``nvidia-smi -pm 1``
* Check GPU utilization: ``nvidia-smi dmon``
* Verify you're using a release build (debug builds are much slower on GPU)

**Build errors with nvcc:**

* Ensure CUDA Toolkit version matches your GPU driver
* Check compute capability: ``nvidia-smi --query-gpu=compute_cap --format=csv``
