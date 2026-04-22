Containers (Docker and Apptainer)
===================================

FLUXOS ships two container recipes under ``containers/`` that install every
system dependency the solver links against. **Compilation is deliberately
left for the user, inside a shell in the running container** — this keeps
the images small, gives you a reproducible dev loop (edit source on the host,
rebuild inside the container), and drops the compiled binary straight onto
the host via a bind mount.

Recommended paths:

* **Docker** — laptops / workstations (``containers/Dockerfile``,
  ``containers/docker-compose.yml``).
* **Apptainer / Singularity** — HPC clusters where Docker is not available
  (``containers/fluxos_apptainer.def`` for CPU,
  ``containers/fluxos_apptainer_cuda.def`` for GPU).

Both paths give the same four-step flow:

1. Install the container runtime and do a sanity check.
2. Build the image (dependencies only — this is cached after the first build).
3. Open a shell inside a container from that image.
4. Inside the shell, configure with CMake, ``make``, and run FLUXOS.

Mount layout
------------

Both compose / ``.def`` files bind-mount **the entire repository** at a
single path inside the container:

.. list-table::
   :header-rows: 1
   :widths: 25 25 15 35

   * - Host path
     - Container path
     - Access
     - Notes
   * - repo root (``.``)
     - ``/work`` (Docker) or ``/src`` (Apptainer)
     - read-write
     - sources, inputs (``Working_example/``), compiled binary (``bin/``), and every ``Results*/`` folder are the same files on the host and in the container

Mounting the whole repo (rather than picking individual subfolders) means
every file the simulation writes — no matter which ``OUTPUT_FOLDER`` the
modset names — lands on the host automatically. ``Results/``,
``Results_river_30h/``, ``Results_trimesh_soil/``, etc. all work without
changes to the compose / ``.def`` file.

The key CMake flag in the four-step flow is
``-DCMAKE_RUNTIME_OUTPUT_DIRECTORY=/work/bin`` (Docker) or ``/src/bin``
(Apptainer) — it redirects the compiled binary to the bind-mounted ``bin/``
folder so you can run it from the host as ``./bin/fluxos`` after exiting
the shell.

Docker
------

**1. Install Docker**

* macOS / Windows: install Docker Desktop
  (https://www.docker.com/products/docker-desktop).
* Linux: install Docker Engine
  (https://docs.docker.com/engine/install/).

Sanity check the daemon is up:

.. code-block:: bash

   docker --version && docker compose version

**2. Build the FLUXOS dependency image**

From the repo root:

.. code-block:: bash

   docker compose -f containers/docker-compose.yml build

This produces the ``fluxos:latest`` image (~300 MB) containing CMake, GCC,
Armadillo, nlohmann/json, HDF5, OpenMP, and optionally OpenMPI (pass
``--build-arg USE_MPI=ON`` to include it).

**3. Open a shell inside the container**

.. code-block:: bash

   docker compose -f containers/docker-compose.yml run --rm fluxos

You land in ``/work`` with the host repo bind-mounted there. A baked-in
copy of the source at ``/opt/fluxos`` also ships in the image as a
fallback for image-only workflows, but the recommended path is to compile
from ``/work`` directly so host edits are picked up automatically.

**4. Compile and run (inside the container shell)**

.. code-block:: bash

   cd /work && mkdir -p build && cd build
   cmake -DMODE_release=ON -DUSE_TRIMESH=ON \
         -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=/work/bin /work
   make -j$(nproc)

   cd /work
   ./bin/fluxos Working_example/modset_trimesh.json

When you exit the shell (``exit`` / Ctrl-D) the container is removed (thanks
to ``--rm``), but the binary at ``bin/fluxos``, the ``build/`` tree, and
every simulation-output file stay on the host — they were always on the
host through the bind mount.

Re-compile after editing source — no image rebuild needed
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Because the repo is bind-mounted at ``/work``, edits made on the host
(in your editor) are immediately visible inside the container. Just
re-enter the container shell and ``make -j$(nproc)`` from ``build/`` —
no image rebuild required. Run ``docker compose ... build`` only when you
change the Dockerfile or want a fresh dependency layer.

Apptainer / Singularity (HPC)
-----------------------------

Apptainer (the OCI-compatible successor to Singularity) is the container
runtime most HPC clusters support. FLUXOS ships two recipes:

* ``containers/fluxos_apptainer.def`` — CPU image (OpenMP + optional MPI).
* ``containers/fluxos_apptainer_cuda.def`` — NVIDIA CUDA-enabled image for GPU
  runs (based on an ``nvidia/cuda`` image; requires ``--nv`` at runtime).

**1. Activate Apptainer**

Most HPC systems ship it as an environment module:

.. code-block:: bash

   module load apptainer   # or: module load singularity
   apptainer --version

On a local workstation, install from
https://apptainer.org/docs/admin/main/installation.html.

**2. Build the .sif image**

From the repo root (or the login node of your cluster):

.. code-block:: bash

   apptainer build fluxos.sif containers/fluxos_apptainer.def

   # GPU variant (CUDA base image):
   apptainer build fluxos_cuda.sif containers/fluxos_apptainer_cuda.def

.. note::

   Apptainer does not require root to run images, but building ``.sif``
   images with ``--fakeroot`` may require admin configuration on your
   cluster. When that's a blocker, build locally on your workstation and
   ``scp`` the resulting ``.sif`` file to the cluster.

**3. Open a shell inside the container**

Bind-mount the repo at ``/src`` so compilation output lands on the host
filesystem:

.. code-block:: bash

   apptainer shell --bind "$PWD:/src" fluxos.sif

   # GPU variant — expose the host's NVIDIA driver with --nv:
   apptainer shell --nv --bind "$PWD:/src" fluxos_cuda.sif

**4. Compile and run (inside the container shell)**

.. code-block:: bash

   cd /src && mkdir -p build && cd build
   cmake -DMODE_release=ON -DUSE_TRIMESH=ON \
         -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=/src/bin /opt/fluxos
   make -j$(nproc)

   cd /src
   ./bin/fluxos Working_example/modset_trimesh.json

For CUDA GPU builds, add ``-DUSE_CUDA=ON``; for MPI, add ``-DUSE_MPI=ON`` and
invoke with ``mpirun -n <N>``. See :doc:`HPC` for SLURM submission templates.

Container build flags (summary)
-------------------------------

Once inside the shell, the same CMake flags that apply to a native build
apply here. Common combinations:

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Goal
     - CMake flags
   * - OpenMP only, regular mesh
     - ``-DMODE_release=ON``
   * - OpenMP + triangular mesh
     - ``-DMODE_release=ON -DUSE_TRIMESH=ON``
   * - MPI + OpenMP (HPC)
     - ``-DMODE_release=ON -DUSE_MPI=ON``
   * - CUDA GPU
     - ``-DMODE_release=ON -DUSE_CUDA=ON``
   * - Full hybrid (GPU + MPI + triangular)
     - ``-DMODE_release=ON -DUSE_TRIMESH=ON -DUSE_CUDA=ON -DUSE_MPI=ON``

Always include
``-DCMAKE_RUNTIME_OUTPUT_DIRECTORY=/work/bin`` (Docker) or ``/src/bin``
(Apptainer) so the binary lands on the bind-mounted host path.

Troubleshooting
---------------

**"Cannot connect to the Docker daemon"** — Docker Desktop is not running,
or your user is not in the ``docker`` group (Linux native install).

**"apt-get update ... Hash Sum mismatch" during image build** — transient
mirror issue. Rerun ``docker compose ... build``.

**``cmake`` cannot find Armadillo inside the container** — the image already
ships ``libarmadillo-dev``; you most likely ran ``cmake`` from the wrong
directory. Always invoke it from ``/work/build`` (Docker) or
``/src/build`` (Apptainer), pointing the source-tree argument at the
bind-mounted repo (``/work`` or ``/src``).

**Binary ends up in ``build/bin/`` instead of the host ``bin/``** — you
forgot ``-DCMAKE_RUNTIME_OUTPUT_DIRECTORY=/work/bin`` (or ``/src/bin``).
Either re-run CMake with the flag, or copy the binary manually:
``cp /work/build/bin/fluxos /work/bin/``.

**Simulation runs but I cannot find the results** — the output folder is
taken from the modset's ``OUTPUT.OUTPUT_FOLDER`` field (e.g.
``"Results_river_30h/"``), resolved **relative to the current working
directory at run time**. Because the whole repo is bind-mounted at
``/work`` (Docker) or ``/src`` (Apptainer), running from that path means
the results appear alongside ``bin/``, ``Working_example/``, and the rest
of your checkout on the host. If you ``cd`` somewhere else before running,
they land there instead.

**Apptainer ``--nv`` reports "Failed to find NVIDIA driver"** — the host has
no NVIDIA driver (or an incompatible version). Use the CPU ``.sif`` instead.
