#!/bin/bash
# analyse.slurm -- collect the Slurm accounting of the engine jobs, then run 09_ANALYSE_SIX_ENGINES.jl.
#
# submit_all.sh submits it with --dependency=afterany on every engine job, so it runs (and compares
# whatever engines succeeded) even if one engine failed. By hand:
#   sbatch --export=ALL,CMR_PROJECT_NAME=TestDomain_N8 ../jobs/analyse.slurm
#SBATCH --job-name=cmr_analyse
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.out

set -euo pipefail
source "${CMR_BASE:-$HOME/CoralModelReview}/jobs/common.sh"
echo "[$(date '+%F %T')] analysing $CMR_PROJECT_NAME on $(hostname)"

bash "$JOBS_DIR/collect_benchmarks.sh" "$CMR_PROJECT_NAME" || echo "WARNING: could not collect Slurm accounting; benchmarks will use in-process numbers only"

load_singularity
container julia --startup-file=no /opt/coral/Code/09_ANALYSE_SIX_ENGINES.jl

echo "[$(date '+%F %T')] done: $RUNS/$CMR_PROJECT_NAME/FIGURES.md"
