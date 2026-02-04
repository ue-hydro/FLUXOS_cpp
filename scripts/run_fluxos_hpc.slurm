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

# ============================================================
# FLUXOS MPI+OpenMP Hybrid Job Script for HPC Clusters
# ============================================================
#
# This script runs FLUXOS with hybrid MPI+OpenMP parallelization:
# - MPI for distributed memory (across nodes)
# - OpenMP for shared memory (within each node)
#
# Configuration:
# - 4 nodes x 32 MPI tasks = 128 total MPI processes
# - 2 OpenMP threads per MPI task
# - Total parallelism: 128 x 2 = 256 cores
#
# Adjust --nodes, --ntasks-per-node, and --cpus-per-task
# based on your cluster and domain size.
# ============================================================

# Load required modules (adjust for your HPC system)
module purge
module load gcc/11.2.0
module load openmpi/4.1.1
module load armadillo/11.0
module load hdf5/1.12.0

# Set OpenMP threads per MPI task
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Optimize OpenMP settings
export OMP_PROC_BIND=close
export OMP_PLACES=cores

# MPI settings for hybrid runs
export OMPI_MCA_btl=^openib  # Disable InfiniBand if not available

# Print job information
echo "============================================"
echo "FLUXOS MPI+OpenMP Hybrid Run"
echo "============================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Nodes: ${SLURM_NNODES}"
echo "Tasks per node: ${SLURM_NTASKS_PER_NODE}"
echo "Total MPI tasks: ${SLURM_NTASKS}"
echo "OpenMP threads per task: ${OMP_NUM_THREADS}"
echo "Total cores: $((SLURM_NTASKS * OMP_NUM_THREADS))"
echo "============================================"

# Path to FLUXOS executable and input file
FLUXOS_EXE="./build/bin/fluxos_mpi"
INPUT_FILE="./input/modset.json"

# Check if executable exists
if [ ! -f "${FLUXOS_EXE}" ]; then
    echo "Error: FLUXOS executable not found at ${FLUXOS_EXE}"
    echo "Please build with: cmake -DMODE_release=ON -DUSE_MPI=ON .. && make"
    exit 1
fi

# Check if input file exists
if [ ! -f "${INPUT_FILE}" ]; then
    echo "Error: Input file not found at ${INPUT_FILE}"
    exit 1
fi

# Create output directory
OUTPUT_DIR="./Results_${SLURM_JOB_ID}"
mkdir -p ${OUTPUT_DIR}

# Run FLUXOS with MPI
echo "Starting FLUXOS at $(date)"
srun --mpi=pmix ${FLUXOS_EXE} ${INPUT_FILE}

# Check exit status
if [ $? -eq 0 ]; then
    echo "FLUXOS completed successfully at $(date)"
else
    echo "FLUXOS failed with exit code $?"
fi

echo "============================================"
echo "Job completed"
echo "============================================"
