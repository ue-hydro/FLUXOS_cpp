<img src="uevora_logo.png" alt="Universidade de Évora" width="350">

# FLUXOS

A 2-D shallow-water flood simulation model (C++), with optional triangular
meshes, MPI + OpenMP + CUDA parallelism, soil infiltration, and
advection–dispersion transport.

📖 **[Documentation on ReadTheDocs](https://fluxos-cpp.readthedocs.io)** · 📽 **[Slide decks in `in_a_nutshell/`](in_a_nutshell/)**

## Table of Contents

* [Quick start](#quick-start)
* [Workflow](#workflow)
* [Repository layout](#repository-layout)
* [Introduction](#introduction)
* [Branches](#branches)
* [Reading material](#reading-material)

## Quick start

Three ways to build FLUXOS, in order of friction:

**1. Docker (zero host dependencies)** — see [`in_a_nutshell/2_How_to_Compile_and_Run_FLUXOS.html`](in_a_nutshell/2_How_to_Compile_and_Run_FLUXOS.html)

```bash
# Build the container image (installs deps, ships source)
docker compose -f containers/docker-compose.yml build
# Open a shell inside the container (the repo is bind-mounted at /work)
docker compose -f containers/docker-compose.yml run --rm fluxos
# Inside the shell — compile, then run
cd /work && mkdir -p build && cd build
cmake -DMODE_release=ON \
      -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=/work/bin /work
make -j$(nproc)
cd /work && ./bin/fluxos Working_example/modset.json
```

**2. Apptainer (for HPC)** — see [`containers/fluxos_apptainer.def`](containers/fluxos_apptainer.def) (CPU+MPI+OpenMP) or [`fluxos_apptainer_cuda.def`](containers/fluxos_apptainer_cuda.def) (GPU+MPI)

```bash
apptainer build --fakeroot fluxos_cpu.sif containers/fluxos_apptainer.def
apptainer shell --bind $PWD:/src fluxos_cpu.sif
# Inside — compile + run:
cd /src && mkdir -p build && cd build
cmake -DMODE_release=ON -DUSE_MPI=ON \
      -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=/src/bin /opt/fluxos
make -j$(nproc)
cd /src && ./bin/fluxos_mpi Working_example/modset.json
```

**3. Native CMake build**

```bash
mkdir build && cd build
cmake -DMODE_release=ON ..   # see the deck for all flags
make -j$(nproc)
cd .. && ./build/bin/fluxos Working_example/modset.json
```

Requires CMake ≥ 3.18, C++17, Armadillo, nlohmann/json (bundled), HDF5, OpenMP.

## Workflow

```
 ┌─────────────────────────┐    ┌───────────────┐    ┌──────────────────────────┐
 │ 1_Model_Config template │ ─► │ FLUXOS solver │ ─► │ 2_Read_Outputs template  │
 │  (DEM + mesh + modset)  │    │  (Docker /    │    │  (stats report + viewer) │
 │                         │    │   Apptainer / │    │                          │
 │  HTML config report     │    │   native)     │    │   HTML results report    │
 └─────────────────────────┘    └───────────────┘    └──────────────────────────┘
```

* **Preprocessing** — edit and run [`supporting_scripts/1_Model_Config/model_config_template.py`](supporting_scripts/1_Model_Config/model_config_template.py).
  Handles GeoTIFF → ESRI-ASCII DEM, slope-adaptive Gmsh triangular mesh,
  `modset.json`, and an interactive HTML config report with copy-paste
  Docker commands and an embedded DEM + mesh preview map. Supports
  auto-downloading DEMs from OpenTopography / USGS 3DEP.
* **Post-processing** — edit and run [`supporting_scripts/2_Read_Outputs/read_output_template.py`](supporting_scripts/2_Read_Outputs/read_output_template.py).
  Produces an HTML results report with time-series (volume, flooded area,
  max depth, velocity), maximum-inundation map, flood-hazard
  classification (ARR-2019), depth histogram, and first-inundation map.
  The same folder also hosts [`fluxos_viewer.py`](supporting_scripts/2_Read_Outputs/output_supporting_lib/fluxos_viewer.py)
  (Google-Earth KML / MP4 / WebGL animation exporter).

## Repository layout

| Path | Purpose |
|---|---|
| `src/fluxos/` | C++ solver — regular and triangular mesh paths |
| `containers/` | Dockerfile + Apptainer recipes (CPU / GPU) |
| `supporting_scripts/1_Model_Config/` | User-editable template that builds DEM / mesh / modset.json + HTML config report |
| `supporting_scripts/2_Read_Outputs/` | User-editable template that builds the HTML results report + KML / WebGL exporters |
| `bin/` | Input examples (DEMs, meshes, forcing files, modset JSONs) |
| `in_a_nutshell/` | Interactive HTML slide decks: overview, compile guide, model setup, supporting scripts |
| `wikipage/source/` | ReadTheDocs (Sphinx) sources |
| `supporting_scripts/2_Read_Outputs/fluxos_web/` | Browser-based WebGL animation viewer (template) |
| `Working_example/` | Legacy self-contained example (pre-JSON modset format) kept for reference |

## Introduction

Source code for the FLUXOS model. The original code (named FLUXOS) was
written in Fortran and consisted of the coupling of 2dmb, +QeS2, MODFLOW
and MT3DMS.

**This C++ port**:
* Uses [Armadillo](https://arma.sourceforge.net/) for linear algebra and [nlohmann/json](https://github.com/nlohmann/json) for configuration
* Removed MODFLOW and MT3DMS (no baseflow component at present)
* Integrates the WINTRA algorithm for runoff–soil interactions and nutrient release ([paper](https://onlinelibrary.wiley.com/doi/full/10.1002/hyp.11346))
* Adds OpenMP + MPI + CUDA parallelism and triangular-mesh support via Gmsh

## Branches

**Active**
* `main` — primary branch with latest verified updates
* `development` — integration branch for verifying feature branches before merging with main
* `supporting_scripts` — feature branch for Python / MATLAB preprocessing and analysis scripts

**Archived**
* `adesolver`, `adesolver_wintra` — earlier feature branches, merged into main

If you have a local clone from before the `master → main` rename:
```
git branch -m master main
git fetch -p origin
git branch -u origin/main main
```

## Reading material

**Theoretical background (original FLUXOS)**
* [EMS paper](https://www.sciencedirect.com/science/article/pii/S1364815216306193?via%3Dihub)
* [PhD thesis](https://scholarbank.nus.edu.sg/handle/10635/124183)

**Applications (original FLUXOS)**
* [STC paper](https://www.sciencedirect.com/science/article/pii/S0169772216300948?via%3Dihub) · [JCH paper](https://www.sciencedirect.com/science/article/pii/S0169772216300948?via%3Dihub) · [JAWRA](https://onlinelibrary.wiley.com/doi/full/10.1111/1752-1688.12316)

**FLUXOS (C++)**
* [Poster — snowmelt flooding in the Canadian Prairies](https://www.researchgate.net/publication/333324452_Hydrodynamic_modelling_of_snowmelt_flooding_events_and_nutrient_transport_in_the_Canadian_Prairies_using_the_FLUXOS_model)
