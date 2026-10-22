#!/bin/bash
# sa_post.slurm -- collect a model's chunks into results tables (06) and draw the figures (07), one domain.
#
# Submitted by sa_submit.sh after the chunk arrays. By hand, from ~/CoralModelReview/sa/logs:
#   sbatch --export=ALL,MODEL=coralpop,DOMAIN=NineReefs_SobolDHW ../jobs/sa_post.slurm
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.out

set -euo pipefail
source "${SA_BASE:-$HOME/CoralModelReview/sa}/jobs/sa_common.sh"

DOMAIN="${DOMAIN:?set DOMAIN}"
echo "[$(date '+%F %T')] collect + analyse $MODEL $DOMAIN on $(hostname)"
load_singularity
container julia --startup-file=no "$SA_CODE/06_collect.jl" "$MODEL" "$DOMAIN"
container julia --startup-file=no "$SA_CODE/07_analyse.jl" "$MODEL" "$DOMAIN"
container julia --startup-file=no "$SA_CODE/08_indicator_outputs.jl" "$MODEL" "$DOMAIN"
echo "[$(date '+%F %T')] finished: $(run_dir "$DOMAIN")/FIGURES.md"
