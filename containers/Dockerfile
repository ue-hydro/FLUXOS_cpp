# FLUXOS C++ — Docker image (dependency container)
#
# This image installs everything FLUXOS NEEDS to build (compilers, CMake,
# Armadillo, HDF5, nlohmann/json, OpenMP, optional OpenMPI). It does NOT
# compile FLUXOS during the image build — compilation happens interactively
# inside the container, against the source tree that is bind-mounted from
# the host at /work (see docker-compose.yml).
#
# A copy of the source at /opt/fluxos is also shipped as a fallback (useful
# if the image is used standalone, without the compose bind mount).
#
# The intended workflow — matches the in_a_nutshell deck:
#
#   1. Install / start Docker on the host
#   2. Build the image:
#        docker compose -f containers/docker-compose.yml build
#   3. Open a shell inside the container:
#        docker compose -f containers/docker-compose.yml run --rm fluxos
#   4. Compile FLUXOS inside the shell. The whole repo is bind-mounted at
#      /work, so compiling from there means host edits are immediately
#      visible and every file the binary writes (binary itself, simulation
#      results, etc.) lands on the host filesystem automatically:
#        cd /work && mkdir -p build && cd build
#        cmake -DMODE_release=ON -DUSE_TRIMESH=ON \
#              -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=/work/bin /work
#        make -j$(nproc)
#      Then run — still inside the shell, from /work so modset paths resolve:
#        cd /work
#        ./bin/fluxos Working_example/modset_river_30h.json
#
# Keeping the compile step explicit lets developers iterate on source
# without rebuilding the image, and makes the dependency vs. build
# separation obvious.
#
# CUDA is not supported by this Dockerfile. For GPU runs use
# containers/fluxos_apptainer_cuda.def (NVIDIA base image + --nv at runtime).

ARG UBUNTU_VERSION=24.04

FROM ubuntu:${UBUNTU_VERSION}

ARG DEBIAN_FRONTEND=noninteractive
ARG USE_MPI=OFF

# --- Step 2: dependencies -------------------------------------------------
# Everything FLUXOS links against + the toolchain needed to compile it.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        git \
        ca-certificates \
        libarmadillo-dev \
        nlohmann-json3-dev \
        libhdf5-dev \
        libomp-dev \
    && if [ "$USE_MPI" = "ON" ]; then \
         apt-get install -y --no-install-recommends libopenmpi-dev openmpi-bin; \
       fi \
    && rm -rf /var/lib/apt/lists/*

# --- Source is copied in so it's available for the in-container compile ---
# (Alternative: omit this COPY and bind-mount the repo at run time via
# docker-compose. We copy by default so the image is self-sufficient.)
COPY . /opt/fluxos

# Working dir — docker-compose bind-mounts Working_example/ (ro inputs),
# bin/ (rw compiled binary output), and Results/ (rw simulation outputs) here.
WORKDIR /work

# Default: drop the user into a shell — they'll compile FLUXOS here (step 4).
CMD ["/bin/bash"]
