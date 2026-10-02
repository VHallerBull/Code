#!/bin/bash
# reefmod_B.slurm -- FALLBACK ReefMod route, stage B (HOST MATLAB, not the container): one array task
# per scenario, running ReefModInterface.jl's driver (run_reefmod_scenario.m) on the inputs stage A
# wrote. submit_all.sh submits it as --array=1-N%K (K = REEFMOD_ARRAY_CONCURRENCY, default 4: each
# task holds one MATLAB licence while it runs).
#
# Uses the copies uploaded to ~/CoralModelReview/reefmod_host/ (the host cannot see inside the image):
#   reefmod_host/ReefMod-GBR   the fork with ReefModInterface's hooks (staged by stage_sources.sh)
#   reefmod_host/driver        run_reefmod_scenario.m
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%A_%a.out

set -euo pipefail
source "${CMR_BASE:-$HOME/CoralModelReview}/jobs/common.sh"

TASK="${SLURM_ARRAY_TASK_ID:?run as an array job}"
export REEFMOD_JOB_DIR="$RUNS/$CMR_PROJECT_NAME/reefmod/run"
export REEFMOD_GBR_DIR="$CMR_BASE/reefmod_host/ReefMod-GBR"
DRIVER_DIR="$CMR_BASE/reefmod_host/driver"
[ -f "$REEFMOD_JOB_DIR/inputs/scenario_$TASK.mat" ] || { echo "no input for task $TASK in $REEFMOD_JOB_DIR/inputs" >&2; exit 1; }
[ -f "$DRIVER_DIR/run_reefmod_scenario.m" ] || { echo "upload reefmod_host/ first (RUNNING_ON_HPC.md 6e)" >&2; exit 1; }

module load "$MATLAB_MODULE"
echo "[$(date '+%F %T')] ReefMod task $TASK on $(hostname): $(command -v matlab)"
timed "$RUNS/$CMR_PROJECT_NAME/time_reefmod_B_$TASK.log" \
    matlab -nodisplay -nosplash -batch "addpath('$DRIVER_DIR'); run_reefmod_scenario($TASK)"
