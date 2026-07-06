#!/bin/bash
# =====================================================================
# submit_task3.sh — Slurm job array for Task 3 Simulation Study
# AUTH Aristotelis HPC Cluster
#
# USAGE:
#   # Submit one config at a time:
#   CONFIG_IDX=1 sbatch submit_task3.sh   # S1_N100_pmax5
#   CONFIG_IDX=2 sbatch submit_task3.sh   # S1_N100_pmax10
#   ...
#
#   # Or pass CONFIG_IDX directly:
#   sbatch --export=ALL,CONFIG_IDX=3 submit_task3.sh
#
# MONITORING:
#   squeue -u $USER
#   tail -f slurm-<JOBID>_<TASKID>.out
#
# AFTER COMPLETION (merge partial files):
#   matlab -nodisplay -nosplash -r "merge_hpc_results('S1_N100_pmax5'); exit"
# =====================================================================

# --- Slurm directives ---
#SBATCH --job-name=task3_svar
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --array=1-20
#SBATCH --output=logs/task3_%A_%a.out
#SBATCH --error=logs/task3_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ddantono@auth.gr

# --- Environment ---
echo "============================================"
echo "Job ID:       $SLURM_JOB_ID"
echo "Array Task:   $SLURM_ARRAY_TASK_ID / $SLURM_ARRAY_TASK_MAX"
echo "Node:         $SLURMD_NODENAME"
echo "Config index: $CONFIG_IDX"
echo "Start time:   $(date)"
echo "============================================"

# Create logs directory if it doesn't exist
mkdir -p logs
mkdir -p results/task3/partial

# --- Load required modules ---
module purge
module load matlab/R2025a
module load gcc/14.2.0
module load r/4.4.1          # adjust to available R version

# Verify R is available
which Rscript || { echo "ERROR: Rscript not found after module load"; exit 1; }
echo "Rscript: $(which Rscript)"
echo "R version: $(Rscript --version 2>&1)"

# --- Export variables for MATLAB ---
export HPC_BATCH_ID=$SLURM_ARRAY_TASK_ID
export HPC_CONFIG_IDX=${CONFIG_IDX:-1}

# --- Framework path ---
FRAMEWORK_ROOT="/home/d/ddantono/svar"

# --- Run MATLAB ---
echo "Running MATLAB batch..."
matlab -nodisplay -nosplash -nojvm -r \
    "addpath('${FRAMEWORK_ROOT}/simulations'); \
     run_sim_task3_hpc; \
     exit" \
    2>&1

EXIT_CODE=$?
echo "MATLAB exit code: $EXIT_CODE"
echo "End time: $(date)"

if [ $EXIT_CODE -ne 0 ]; then
    echo "ERROR: MATLAB job failed with exit code $EXIT_CODE"
    exit $EXIT_CODE
fi

echo "Task $SLURM_ARRAY_TASK_ID complete."
