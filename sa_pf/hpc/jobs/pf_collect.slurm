#!/bin/bash
# pf_collect.slurm -- merge the chunk files of a band into pf_results.jld2 per domain and set (04_collect.jl).
# It lists missing or failed walks; rerun those tasks with pf_chunk.slurm, then this again.
#
# Submitted by pf_submit.sh after the chunk array. By hand, from ~/CoralModelReview/sa/logs:
#   sbatch --export=ALL,MODEL=coralpop,BAND=0.5 ../jobs/pf_collect.slurm
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.out

set -euo pipefail
source "${SA_BASE:-$HOME/CoralModelReview/sa}/jobs/pf_common.sh"

BAND="${BAND:?set BAND (e.g. 0.5)}"
echo "[$(date '+%F %T')] collect $MODEL band $BAND on $(hostname)"
load_singularity
container julia --startup-file=no "$PF_CODE/04_collect.jl" "$MODEL" "$BAND"
echo "[$(date '+%F %T')] finished: $(pf_band_dir "$BAND")"
