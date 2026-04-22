# `bin/` — compiled binaries

This directory is reserved for the **compiled FLUXOS binary** (`fluxos`,
`fluxos_mpi`, or `fluxos_cuda`, depending on the build flags).

## How binaries end up here

### Container builds (Docker / Apptainer)

The container recipes bind-mount this directory into the container and
direct the CMake build to drop the binary here via the
`-DCMAKE_RUNTIME_OUTPUT_DIRECTORY=/work/bin` (Docker) or `/src/bin`
(Apptainer) flag. Run the simulation from the repo root with:

```bash
./bin/fluxos Working_example/modset.json
```

See [`containers/Dockerfile`](../containers/Dockerfile) +
[`containers/docker-compose.yml`](../containers/docker-compose.yml) and
[`containers/fluxos_apptainer.def`](../containers/fluxos_apptainer.def) for
the full 4-step setup.

### Native build

A plain native build puts the binary under `build/bin/` by default (per
`CMakeLists.txt`). Copy or symlink it here for consistency with the
container workflow:

```bash
cp build/bin/fluxos bin/
# or:  ln -sf ../build/bin/fluxos bin/fluxos
```

## What is NOT here

Example simulation inputs (DEMs, meshes, forcing files, modset JSONs) live
in [`../Working_example/`](../Working_example/). Outputs land in
[`../Results/`](../Results/) (git-ignored).
