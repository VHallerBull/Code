#!/bin/bash
# pf_chunk.slurm -- one array task = one range of base samples (whole walks) of the pick-and-freeze
# design(s), 03_run_chunk.jl. Task t runs base samples FIRST + (t - 1) WALKS ... FIRST + t WALKS - 1 of
# every job in JOBS ("band:domain" = varied DHW, "band:domain:zero" = zero DHW).
#
# Submitted by pf_submit.sh. By hand, from ~/CoralModelReview/sa/logs (e.g. rerun tasks 7 and 12):
#   sbatch --array=7,12 --export=ALL,MODEL=coralpop,WALKS=50,FIRST=1,JOBS="0.5:SingleReef_SobolDHW 0.5:SingleReef_SobolDHW:zero" ../jobs/pf_chunk.slurm
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=01:30:00
#SBATCH --output=%x_%A_%a.out

set -euo pipefail
source "${SA_BASE:-$HOME/CoralModelReview/sa}/jobs/pf_common.sh"

WALKS="${WALKS:?set WALKS (base samples per task)}"
FIRST="${FIRST:-1}"
JOBS="${JOBS:?set JOBS (e.g. 0.5:SingleReef_SobolDHW 0.5:SingleReef_SobolDHW:zero)}"
TASK="${SLURM_ARRAY_TASK_ID:?run as an array job}"
a=$(( FIRST + (TASK - 1) * WALKS ))
b=$(( a + WALKS - 1 ))
echo "[$(date '+%F %T')] $MODEL base samples $a-$b of: $JOBS, on $(hostname), job ${SLURM_ARRAY_JOB_ID:-?}_$TASK"

load_singularity
band="${JOBS%%:*}"
logdir="$(pf_band_dir "$band")/logs"
mkdir -p "$logdir"
# shellcheck disable=SC2086  # JOBS is a list
timed "$logdir/time_$(printf '%05d' "$a").log" \
    container julia --startup-file=no "$PF_CODE/03_run_chunk.jl" "$MODEL" "--bases=$a:$b" $JOBS

echo "[$(date '+%F %T')] finished"
