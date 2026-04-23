Windows Installation (native)
===============================

Windows support is offered via three routes, in order of preference:

1. **WSL2 (recommended)** — install Ubuntu via the Microsoft Store, then
   follow :doc:`Linux`. This is the smoothest path and matches what the
   project is tested against.
2. **Docker Desktop** — the container workflow documented in
   :doc:`Containers` works natively on Windows once Docker Desktop is
   installed (it runs Linux containers via WSL2 under the hood).
3. **Native MSYS2 / MinGW-w64** — possible but not actively tested.

Route 1 — WSL2 + Ubuntu
-----------------------

.. code-block:: powershell

   wsl --install -d Ubuntu

Open the Ubuntu shell, then follow the native Linux quick-start:

.. code-block:: bash

   sudo apt-get install cmake libarmadillo-dev libomp-dev \
                        nlohmann-json3-dev libhdf5-dev
   git clone https://github.com/ue-hydro/FLUXOS_cpp.git
   cd FLUXOS_cpp
   mkdir build && cd build
   cmake -DMODE_release=ON ..
   make -j$(nproc)
   cd ..
   ./build/bin/fluxos Working_example/modset_trimesh.json

WSL2 sees your Windows filesystem under ``/mnt/c/…`` and its own Linux
filesystem via the ``\\wsl$\Ubuntu`` share in Windows Explorer.

Route 2 — Docker Desktop
------------------------

Install Docker Desktop for Windows
(https://www.docker.com/products/docker-desktop), then use the container
workflow in :doc:`Containers`. All commands in that page work as-is in
PowerShell.

Route 3 — Native MSYS2 / MinGW-w64
-----------------------------------

The dependencies (Armadillo, OpenMP, HDF5, nlohmann/json, CMake) are
installable via pacman under MSYS2. This route is rarely used and not
regularly tested; expect to debug PATH / linker issues on your own. For
most Windows users the WSL2 or Docker routes above are simpler.
