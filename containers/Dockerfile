# FLUXOS C++ — Docker image
#
# Build (from repo root):
#   docker compose -f containers/docker-compose.yml build
#
# Or directly:
#   docker build -f containers/Dockerfile -t fluxos:latest .
#
# With MPI enabled:
#   docker build -f containers/Dockerfile --build-arg USE_MPI=ON -t fluxos:mpi .
#
# Note: CUDA builds are NOT supported by this Dockerfile. For GPU runs on HPC
# use containers/fluxos_apptainer_cuda.def instead (requires an NVIDIA base image
# and --nv / --gpus at runtime).

ARG UBUNTU_VERSION=24.04

# -----------------------------------------------------------------------------
# Stage 1: builder — compile FLUXOS
# -----------------------------------------------------------------------------
FROM ubuntu:${UBUNTU_VERSION} AS builder

ARG DEBIAN_FRONTEND=noninteractive
ARG USE_TRIMESH=ON
ARG USE_MPI=OFF
ARG USE_OPENWQ=OFF
ARG MODE_release=ON
ARG JOBS

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

WORKDIR /opt/fluxos

# Copy the full source tree (docker-compose build context is repo root)
COPY . /opt/fluxos

# Configure + build
RUN mkdir -p build && cd build \
    && cmake .. \
        -DCMAKE_BUILD_TYPE=Release \
        -DMODE_release=${MODE_release} \
        -DUSE_TRIMESH=${USE_TRIMESH} \
        -DUSE_MPI=${USE_MPI} \
        -DUSE_OPENWQ=${USE_OPENWQ} \
    && make -j${JOBS:-$(nproc)}

# Identify the produced binary name (varies with build flags — see CMakeLists.txt)
# and copy it to a stable path for the runtime stage.
RUN set -e; \
    if [ "$USE_MPI" = "ON" ]; then \
        cp /opt/fluxos/build/bin/fluxos_mpi /opt/fluxos/build/bin/fluxos_final; \
    else \
        cp /opt/fluxos/build/bin/fluxos /opt/fluxos/build/bin/fluxos_final; \
    fi

# -----------------------------------------------------------------------------
# Stage 2: runtime — minimal image with just the binary + shared libraries
# -----------------------------------------------------------------------------
FROM ubuntu:${UBUNTU_VERSION} AS runtime

ARG DEBIAN_FRONTEND=noninteractive
ARG USE_MPI=OFF

# Runtime shared libraries. libgomp1 ships with g++ but is not present in the
# slim runtime base — OpenMP is always linked (see CMakeLists.txt line ~216),
# so libgomp1 is always required.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libarmadillo-dev \
        libhdf5-dev \
        libgomp1 \
    && if [ "$USE_MPI" = "ON" ]; then \
         apt-get install -y --no-install-recommends libopenmpi-dev openmpi-bin; \
       fi \
    && rm -rf /var/lib/apt/lists/*

# Copy the built binary
COPY --from=builder /opt/fluxos/build/bin/fluxos_final /opt/fluxos/bin/fluxos

# Working directory — docker-compose bind-mounts inputs/outputs here.
# Relative paths inside the modset JSON resolve against /work.
WORKDIR /work

ENTRYPOINT ["/opt/fluxos/bin/fluxos"]
CMD ["bin/modset.json"]
