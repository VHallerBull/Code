#!/bin/bash
# sa_chunk.slurm -- one array task = one chunk of a model's Sobol' design on one domain (05_run_chunk.jl).
#
# Submitted by sa_submit.sh. By hand, from ~/CoralModelReview/sa/logs:
#   sbatch --array=1-8 --export=ALL,MODEL=coralpop,DOMAIN=NineReefs_SobolDHW,SET=main,CHUNK_SIZE=1024 ../jobs/sa_chunk.slurm
#
# Rerun single missing chunks with --array=3,7 (06_collect.jl lists the missing rows).
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=00:45:00
#SBATCH --output=%x_%A_%a.out

set -euo pipefail
source "${SA_BASE:-$HOME/CoralModelReview/sa}/jobs/sa_common.sh"

DOMAIN="${DOMAIN:?set DOMAIN}"
SET="${SET:-main}"
CHUNK_SIZE="${CHUNK_SIZE:-1024}"
CHUNK="${SLURM_ARRAY_TASK_ID:?run as an array job}"
echo "[$(date '+%F %T')] $MODEL $DOMAIN $SET chunk $CHUNK (size $CHUNK_SIZE) on $(hostname), job ${SLURM_ARRAY_JOB_ID:-?}_$CHUNK"

load_singularity
mkdir -p "$(run_dir "$DOMAIN")/chunks"
timed "$(run_dir "$DOMAIN")/chunks/time_${SET}_$(printf '%04d' "$CHUNK").log" \
    container julia --startup-file=no "$SA_CODE/05_run_chunk.jl" "$MODEL" "$DOMAIN" "$SET" "$CHUNK" "$CHUNK_SIZE"

echo "[$(date '+%F %T')] finished"
